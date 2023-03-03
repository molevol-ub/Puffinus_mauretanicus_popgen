#!/bin/bash
#$ -cwd
#$ -R y
#$ -e covgenome.err
#$ -o covgenome.out
#$ -q h13.q
#$ -pe ompi255h13 6
#$ -V                    #export environment var
#$ -N covgenome             #name Job

# Script to calculate windowed estimates of heterozygosity (taking into account coverage filtered regions)

cd /users-d3/jferrer/gizquierdo/TFM/conservation_genomics/het/PopGenome/DEF_HET
VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP_gatk_filtered.vcf.gz

# First list the individuals and the 1/2x and 2x of their average coverage (here done in a very crappy way)

list=(ALT78 CZA11 ILA13 ILA2 PORQ TZE1 M8 M1 M5 M20 M4 G3 M2 M3 M18 M13 M19 M14 M12 M11 M17 M10 M21 G12 G9 G10 G11 G14 G15 M16 G4 M6 COP1 LT2)
min_cov=(4 4 4 4 4 4 4 3 4 5 4 4 4 4 5 4 4 5 4 4 5 4 4 4 5 4 4 4 4 4 4 4 5 4)
max_cov=(18 15 14 17 15 15 15 14 14 19 15 15 18 17 20 14 15 22 18 17 19 17 14 15 22 17 14 16 16 16 15 15 18 16)

# And the mappability/repeatability masked regions, if you have masked them

bedfile=/users-d3/jferrer/gizquierdo/TFG/genome_scans/masking/final_masked.bed

counter=0

for ind in ${list[*]}
do

# Filter for each one of the inds coverage limits (set at 1/2x and 2x beforehand)

vcftools --gzvcf $VCF --min-meanDP ${min_cov[$counter]} --max-meanDP ${max_cov[$counter]} --min-alleles 2 --max-alleles 2 --minQ 30 --bed $bedfile --indv $ind --remove-filtered-all --recode --out $ind.masked
mv $ind.masked.recode.vcf $ind.masked.vcf

bgzip $ind.masked.vcf
tabix $ind.masked.vcf.gz

# Filter monomorphic sites

bcftools view -Oz -o $ind.masked.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $ind.masked.vcf.gz
tabix $ind.masked.nomono.vcf.gz

# And calculate pi (accurate estimate of heterozygosity) using vcftools - here done for 25k windows

vcftools --gzvcf ./DEF_HET/$ind.masked.nomono.vcf.gz --window-pi 25000 --window-pi-step 25000 --out $ind

# Now let's generate file of regions not masked by coverage (max and min set beforehand)

bedtools genomecov -ibam /users-d3/jferrer/pmau_popgen/SNP_calling/bams/clean_bams/$ind.sorted.dups.bam -bg | awk -v min=${min_cov[$counter]} -v max=${max_cov[$counter]} '$4 > min && $4 < max' > $ind.cov_masked.bed

# Join all overlapping features of the bedfile (if not the file is needlessly long)

cat $ind.cov_masked.bed | cut -f 1,2,3 > $ind.cov_masked.def.bed
bedtools merge -i $ind.cov_masked.def.bed > prova
mv prova $ind.cov_masked.def.bed

# Generate a full masked file that takes into account both the coverage and mappability masks (if needed)

bedtools intersect -a $bedfile -b $ind.cov_masked.def.bed -bed > $ind.full_masked.bed

cat $ind.full_masked.bed | cut -f 1,2,3 > prova
mv prova $ind.full_masked.bed

# Write overlap amount of the unmasked file with the window file (you can obtain it from the vcftools.windowed.pi file with cat your_file | cut -f 1,2,3)

bedtools intersect -a popgenome_windows.bed -b $ind.full_masked.bed -wo -bed > $ind.mask_overlap.bed

counter=$(expr $counter + 1)           # Add to the counter as you iterate through the lists

done

# --------------------------------------------------------------

# If you want 2 inds at a times

ind1=ALT78
ind2=M1

bcftools merge $ind1.masked.vcf.gz $ind2.masked.vcf.gz > $ind1\_$ind2.masked.vcf
bgzip $ind1\_$ind2.masked.vcf
tabix $ind1\_$ind2.masked.vcf.gz

bcftools view -Oz -o $ind1\_$ind2.masked.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' $ind1\_$ind2.masked.vcf.gz
tabix $ind1\_$ind2.masked.nomono.vcf.gz

vcftools --gzvcf $ind1\_$ind2.masked.nomono.vcf.gz --window-pi 25000 --window-pi-step 25000 --out $ind1\_$ind2

bedtools intersect -a $ind1.full_masked.bed -b $ind2.full_masked.bed > $ind1\_$ind2.full_masked.bed
cat $ind1\_$ind2.full_masked.bed | cut -f 1,2,3 > prova
mv prova $ind1\_$ind2.full_masked.bed

# Write overlap of unmasked file with the popgenome window file

bedtools intersect -a popgenome_windows.bed -b $ind1\_$ind2.full_masked.bed -wo -bed > $ind1\_$ind2.mask_overlap.bed

