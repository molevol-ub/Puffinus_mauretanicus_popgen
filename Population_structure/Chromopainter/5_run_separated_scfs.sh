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
  
  Rscript vcf_to_chromopainter_main.R -g ../vcfs/$scf.vcf -o $scf
  
  # 3. Remove spaces from the "haplotype" lines in the haplotype input file
  
  head $scf.chromopainter.inp -n 3 > prova
  tail -n +4 $scf.chromopainter.inp | sed 's/ //g' >> prova
  mv prova $scf.chromopainter.inp
  
done < ../vcfs/scaffolds.txt

#--------------------------------------------------------------------------------------

# 4. Run chromopainter for the desired populations - remember to use the correct idfile

maindir=/users-d3/jferrer/gizquierdo/TFM/chromopainter

mkdir Menorca_noMen
mkdir Menorca_noMenCab

ChromoPainterv2 -g ${scf[(($SGE_TASK_ID))]}\.chromopainter.inp -r ${scf[(($SGE_TASK_ID))]}\.recomrates.txt -t $maindir/Puffinus.auto.nPP.maxmiss80.idfile.txt \
  -f $maindir/popfiles/Puffinus.auto.nPP.maxmiss80.popfile_men2.txt 0 0 -n 33700.840 -M 0.022739 -o Menorca_noMen/${scf[(($SGE_TASK_ID))]}

ChromoPainterv2 -g ${scf[(($SGE_TASK_ID))]}\.chromopainter.inp -r ${scf[(($SGE_TASK_ID))]}\.recomrates.txt -t $maindir/Puffinus.auto.nPP.maxmiss80.idfile.txt \
  -f $maindir/popfiles/Puffinus.auto.nPP.maxmiss80.popfile_men3.txt 0 0 -n 33700.840 -M 0.022739 -o Menorca_noMenCab/${scf[(($SGE_TASK_ID))]}

