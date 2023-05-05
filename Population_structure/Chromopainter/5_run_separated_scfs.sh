#!/bin/bash
#$ -cwd
#$ -R y
#$ -e divide_vcf.err
#$ -o divide_vcf.out
#$ -q h11.q
#$ -pe ompi64h11 1
#$ -V                    #export environment var
#$ -N divide_vcf             #name Job
#$ -t 1-876
#$ -tc 24

workdir=/users-d3/jferrer/gizquierdo/TFM/chromopainter/scaffolds
#mkdir $workdir
cd $workdir

# 0. Make your list of scaffolds

outvcf_nomono=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.nPP.maxmiss80.nomono.vcf
cat $outvcf_nomono | grep -v "^#" | cut -f1 | sort | uniq > scaffolds.txt
scf=($(cat scaffolds.txt))

# 1. Divide your VCF by scaffolds

vcftools --vcf $outvcf_nomono --chr ${scf[(($SGE_TASK_ID))]} --recode-INFO-all --recode --out ${scf[(($SGE_TASK_ID))]}
mv ${scf[(($SGE_TASK_ID))]}\.recode.vcf ${scf[(($SGE_TASK_ID))]}\.vcf

# Now for every scaffold

while read scf; do

  # 2. Transform VCF to chromopainter input (if it works) using the script provided in: https://github.com/sahwa/vcf_to_chromopainter
  
  vcf_to_chromopainter_main.R -g $scf.vcf -o $scf
  
  # 3. Remove spaces from the "haplotype" lines in the haplotype input file
  
  head $scf.chromopainter.inp -n 3 > prova
  tail -n +4 $scf.chromopainter.inp | sed 's/ //g' >> prova
  mv prova $scf.chromopainter.inp
  
  # 4. Modify the input recombination file so that the last recombination of every gene is substituted by "-9" to accomodate the requirements of chromopainter
  
  python correct_recomb_file.py $scf.recomrates.txt prova.txt
  
  # Copy the last line, that is not included in the script
  
  tail -n 1 $scf.recomrates.txt >> prova.txt
  mv prova.txt $scf.recomrates.txt

done < scaffolds.txt
