#!/bin/bash
#$ -cwd
#$ -R y
#$ -e roh_plink.err
#$ -o roh_plink.out
#$ -q h13.q
#$ -pe ompi255h13 12
#$ -V                    #export environment var
#$ -N ROH_plink             #name Job

# List of names
inds=("ALT78" "CZA11" "ILA13" "ILA2" "PORQ" "TZE1" "M8" "M1" "M5" "M20" "M4" "G3" "M2" "M3" "M18" "M13" "M19" "M14" "M12" "M11" "M17" "M10" "M21" "G12" "G9" "G10" "G11" "G14" "G15" "M16" "G4" "M6" "COP1" "LT2")
VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.wPP.merged.nomono.masked.vcf.gz

# Iterate through combinations of ind1 and ind2
for ((i=0; i<${#inds[@]}; i++)); do
    for ((j=i+1; j<${#inds[@]}; j++)); do
        ind1="${inds[i]}"
        ind2="${inds[j]}"
        output_file="${ind2}.${ind1}.intersect.bed"

        # Check if the output file already exists, if not, run the command
        if [ ! -f "$output_file" ]; then
            echo "Running bedtools intersect for $ind1 and $ind2"
            bedtools intersect -a $ind1.bed -b $ind2.bed > "$output_file"
            vcftools --gzvcf $VCF --indv $ind1 --indv $ind2 --bed $output_file --window-pi 25000 --out $ind1.$ind2
        else
            echo "Skipping $ind1 and $ind2, output file already exists: $output_file"
        fi
    done
done

echo "Intersect commands completed!"

rm *intersect*

#--------------------------------------------------------------

while IFS= read -r line; do
    fields=($line)  # Split the line into fields

    col10="${fields[9]}"
    col11="${fields[10]}"

    if [[ $col10 == "0/1"* && $col11 == "0/1"* ]]; then
        col10="0/0${col10:3}"
    elif [[ $col10 == "0/1"* && $col11 == "0/0"* ]]; then
        col10="0/0${col10:3}"
    elif [[ $col10 == "0/1"* && $col11 == "1/1"* ]]; then
        :
    elif [[ $col10 == "0/0"* && $col11 == "0/0"* ]]; then
        :
    elif [[ $col10 == "0/0"* && $col11 == "0/1"* ]]; then
        :
    elif [[ $col10 == "0/0"* && $col11 == "1/1"* ]]; then
        col10="0/1${col10:3}"
    elif [[ $col10 == "1/1"* && $col11 == "0/0"* ]]; then
        col10="0/1${col10:3}"
    elif [[ $col10 == "1/1"* && $col11 == "0/1"* ]]; then
        col10="0/1${col10:3}"
    elif [[ $col10 == "1/1"* && $col11 == "1/1"* ]]; then
        :
    fi

    echo "$col10"
done < prova.vcf > prova.2.vcf

