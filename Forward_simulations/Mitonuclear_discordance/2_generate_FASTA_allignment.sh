#!/bin/bash
#$ -cwd
#$ -R y
#$ -e slim_fasta.err
#$ -o slim_fasta.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N create_fasta             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

# Script to generate a fasta allignment from the VCF resultuing from the forward mtDNA simulations

# 1. Generate the random fasta sequence

# Define variables
output_file="mito.reference.fasta"
sequence_length=16000

# Generate random sequence
random_sequence=$(cat /dev/urandom | tr -dc 'ACGT' | fold -w $sequence_length | head -n 1)

# Write sequence to output file in FASTA format
echo ">RandomSequence" > "$output_file"
echo "$random_sequence" >> "$output_file"

# 2. Substitute the REF position to obtain the reference fasta

# Define variables
vcf_file="Shearwater_mito.def.vcf"
fasta_file="mito.reference.fasta"
output_file="mito.reference.def.fasta"

# Iterate over each line in the VCF file
while IFS=$'\t' read -r chrom pos id ref alt qual filter info format "${individuals[@]}"; do
    if [[ $chrom == \#* ]]; then
        continue  # Skip header lines
    fi
    
    # Extract position and REF nucleotide
    position=$(echo "$pos" | sed 's/ //g')
    nucleotide=$(echo "$ref" | sed 's/ //g')

    # Replace nucleotide at specified position in the fasta file
    sed -i "s/^\(.\{$position\}\).\{1\}/\1$nucleotide/" "$fasta_file"

done < <(grep -v '^#' "$vcf_file")

# Copy the modified fasta to a new file
cp "$fasta_file" "$output_file"

samtools faidx mito.reference.def.fasta
java -jar /disk3/h-user/jferrer/programari/anaconda3/pkgs/picard-2.20.4-0/share/picard-2.20.4-0/picard.jar CreateSequenceDictionary R=$output_file O=mito.reference.def.dict

# 3. Generate a fasta for each individual and allign them

vcf_file="Shearwater_mito.def.vcf"
input_fasta="mito.reference.def.fasta"
output="mito_MSA.fasta"

bgzip $vcf_file
tabix $vcf_file.gz

samples=($(zcat $vcf_file.gz | grep -v "^##" | grep "^#" | cut -f10- | tr '\t' '\n'))

for sample in ${samples[@]}
        do vcftools --gzvcf $vcf_file.gz --indv $sample --non-ref-ac 1 --recode --out $sample; bgzip $sample.recode.vcf; tabix $sample.recode.vcf.gz; gatk FastaAlternateReferenceMaker -R $input_fasta -O $sample.fasta -V $sample.recode.vcf.gz; echo ">"$sample >> mtDNA_MSA.SLIM.fasta ; tail -n +2 $sample.fasta >> mtDNA_MSA.SLIM.fasta
        done

# With this allignment you cna now run BEAST to infer a callibrated phylogeny from it
