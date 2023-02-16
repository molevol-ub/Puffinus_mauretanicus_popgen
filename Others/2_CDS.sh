# Script to extract those SNPs found within the coding regions of scfE (scf7180000013763) - after that check whether they are nonSYN/SYN and whether they share the pattern found in local PCAs
# If they share the PCA pattern and are nonSYN changes, they might be driving adaptation

# 1 Extract all the CDSs found within scfE

CDS_file=CDS_indels.scfD.bed
SNP_file=SNPs_gene_indels.scfE.bed
input_vcf=./indels_scfE.recode.vcf

cat /users-d3/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/augustus.hints.gff3 | grep gene | cut -f 1,4,5 | grep scf7180000013763 > $CDS_file

# 2 Extract all SNPs found within scfE (missingness=80%) that are also included in the CDSs using bedtools intersect

bedtools intersect -a $input_vcf -b $CDS_file -bed > $SNP_file

# 3 Select those that conform to the PCA pattern - perhaps using genotype plot? Or manually in any text editor

# 4 Now you must check whether they are synonymous or not: to do so, first extract the complete gff and fasta files for our scaffold

cat /users-d3/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/augustus.hints.gff3 | grep scf7180000013317 > CDS_scfD.full.gff3
samtools faidx ~/pmau_popgen/genome/genome.fa scf7180000013317 > scf7180000013317.fa

# 5 Finally, predict the effects of your selected variants using VEP

vep --cache --offline -i SNPs_CDS.PC1.scfE.vcf --fasta scf7180000013763.fa --gff full_scfE.gff3 -o SNPs_scfE.PC1.vep
