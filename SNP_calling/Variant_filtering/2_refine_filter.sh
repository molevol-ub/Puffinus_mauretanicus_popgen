#!/bin/bash
#$ -R y
#$ -cwd

source $HOME/.profile

### Script to filter and refine vcfs ###
MASTER=/set/your/master/path
reference=/path/to/genome.fa
WORKING_DIR=/set/your/vcf_folder/path
raw_vcf=/path/to/Puffinus_intersect_raw.auto.recode.vcf
DATASET=Puffinus


## Popmaps ## - These have been filtered for low coverage individuals, make sure all are in one directory
SNP_filtered=$WORKING_DIR/${DATASET}_SNP_filter.vcf
gatk_filter_flag=$WORKING_DIR/${DATASET}_SNP_gatk_flagged.vcf
gatk_filtered=$WORKING_DIR/${DATASET}_SNP_gatk_filtered
allele_filtered=$WORKING_DIR/${DATASET}_SNP.minmax2.mindp4maxdp50.filtered
maxmiss_50_filtered=$WORKING_DIR/${DATASET}_SNP.maxmiss50.filtered
maxmiss_80_filtered=$WORKING_DIR/${DATASET}_SNP.maxmiss80.filtered
maxmiss_100_filtered=$WORKING_DIR/${DATASET}_SNP.maxmiss100.filtered
final_filtered=$WORKING_DIR/${DATASET}_pop_SNP.gatk.bi.miss.maf.final.filtered

POPMAP_DIR=$MASTER/SNP_calling/popmaps

## Processing ...

## Select only snps with the "snp_filter"
gatk --java-options "-Xmx20g" SelectVariants -R $reference -V $raw_vcf --select-type-to-include SNP -O $SNP_filtered

### This gatk step does not actually perform any filtering, it just applies the "snp_filter" tag to SNPs that would pass the filtering
gatk --java-options "-Xmx20g" VariantFiltration -R $reference -V $SNP_filtered -O $gatk_filter_flag --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5" --filter-name "snp_filter"

### This stage actually filters out anything that doesn't have the "snp_filter" tag
vcftools --vcf $gatk_filter_flag --recode --remove-filtered-all --out $gatk_filtered

### Use vcftools to filter remaining SNPS ##

## Retain only biallelic SNPs with a min depth of 5 and max depth of 200 and a minimum quality of 30 and remove related individuals
vcftools --vcf $gatk_filtered.recode.vcf --min-alleles 2 --max-alleles 2 --minDP 4 --maxDP 50 --minQ 30 --remove-indv M22 --remove-indv sacella \
--remove-indv M7 --recode --remove-filtered-all --out $allele_filtered

pop_array=(PM PY PP)

for pop in "${pop_array[@]}"
do

## Split vcf file by population and filter by max missing 50 and 75 (for pop gen analyses)
vcftools --vcf $allele_filtered.recode.vcf --keep $POPMAP_DIR/${pop}.popmap --recode --remove-filtered-all --out ${allele_filtered}.${pop}
vcftools --max-missing 0.8 --vcf $allele_filtered.${pop}.recode.vcf --recode --remove-filtered-all --out ${maxmiss_80_filtered}.${pop}
vcftools --max-missing 0.5 --vcf $allele_filtered.${pop}.recode.vcf --recode --remove-filtered-all --out ${maxmiss_50_filtered}.${pop}
bgzip $maxmiss_50_filtered.${pop}.recode.vcf
tabix $maxmiss_50_filtered.${pop}.recode.vcf.gz
bgzip $maxmiss_80_filtered.${pop}.recode.vcf
tabix $maxmiss_80_filtered.${pop}.recode.vcf.gz

done

## Use bcftools isec to generate vcf files for each species that only contain SNPs present in the 3 or 2 species
bcftools isec -p $WORKING_DIR/merge_50_wPP_dir -n=3 $maxmiss_50_filtered.PP.recode.vcf.gz $maxmiss_50_filtered.PY.recode.vcf.gz \
$maxmiss_50_filtered.PM.recode.vcf.gz

bcftools isec -p $WORKING_DIR/merge_50_noPP_dir -n=2 $maxmiss_50_filtered.PY.recode.vcf.gz \
$maxmiss_50_filtered.PM.recode.vcf.gz

bcftools isec -p $WORKING_DIR/merge_80_wPP_dir -n=3 $maxmiss_80_filtered.PP.recode.vcf.gz $maxmiss_80_filtered.PY.recode.vcf.gz \
$maxmiss_80_filtered.PM.recode.vcf.gz

bcftools isec -p $WORKING_DIR/merge_80_noPP_dir -n=2 $maxmiss_80_filtered.PY.recode.vcf.gz \
$maxmiss_80_filtered.PM.recode.vcf.gz

## Use bcftools merge to merge the vcf files
#p50
bgzip $WORKING_DIR/merge_50_wPP_dir/0000.vcf; tabix $WORKING_DIR/merge_50_wPP_dir/0000.vcf.gz
bgzip $WORKING_DIR/merge_50_wPP_dir/0001.vcf; tabix $WORKING_DIR/merge_50_wPP_dir/0001.vcf.gz
bgzip $WORKING_DIR/merge_50_wPP_dir/0002.vcf; tabix $WORKING_DIR/merge_50_wPP_dir/0002.vcf.gz

bgzip $WORKING_DIR/merge_50_noPP_dir/0000.vcf; tabix $WORKING_DIR/merge_50_noPP_dir/0000.vcf.gz
bgzip $WORKING_DIR/merge_50_noPP_dir/0001.vcf; tabix $WORKING_DIR/merge_50_noPP_dir/0001.vcf.gz

bcftools merge -Oz -o $maxmiss_50_filtered.wPP.merged.vcf.gz --threads 16 \
$WORKING_DIR/merge_50_wPP_dir/0000.vcf.gz $WORKING_DIR/merge_50_wPP_dir/0001.vcf.gz $WORKING_DIR/merge_50_wPP_dir/0002.vcf.gz

bcftools merge -Oz -o $maxmiss_50_filtered.noPP.merged.vcf.gz --threads 16 \
$WORKING_DIR/merge_50_noPP_dir/0000.vcf.gz $WORKING_DIR/merge_50_noPP_dir/0001.vcf.gz

tabix $maxmiss_50_filtered.wPP.merged.vcf.gz
tabix $maxmiss_50_filtered.noPP.merged.vcf.gz

#p80
bgzip $WORKING_DIR/merge_80_wPP_dir/0000.vcf; tabix $WORKING_DIR/merge_80_wPP_dir/0000.vcf.gz
bgzip $WORKING_DIR/merge_80_wPP_dir/0001.vcf; tabix $WORKING_DIR/merge_80_wPP_dir/0001.vcf.gz
bgzip $WORKING_DIR/merge_80_wPP_dir/0002.vcf; tabix $WORKING_DIR/merge_80_wPP_dir/0002.vcf.gz

bgzip $WORKING_DIR/merge_80_noPP_dir/0000.vcf; tabix $WORKING_DIR/merge_80_noPP_dir/0000.vcf.gz
bgzip $WORKING_DIR/merge_80_noPP_dir/0001.vcf; tabix $WORKING_DIR/merge_80_noPP_dir/0001.vcf.gz

bcftools merge -Oz -o $maxmiss_80_filtered.wPP.merged.vcf.gz --threads 16 \
$WORKING_DIR/merge_80_wPP_dir/0000.vcf.gz $WORKING_DIR/merge_80_wPP_dir/0001.vcf.gz $WORKING_DIR/merge_80_wPP_dir/0002.vcf.gz

bcftools merge -Oz -o $maxmiss_80_filtered.noPP.merged.vcf.gz --threads 16 \
$WORKING_DIR/merge_80_noPP_dir/0000.vcf.gz $WORKING_DIR/merge_80_noPP_dir/0001.vcf.gz

tabix $maxmiss_80_filtered.wPP.merged.vcf.gz
tabix $maxmiss_80_filtered.noPP.merged.vcf.gz

## Remove monomorphic sites ##
bcftools view -Oz -o $maxmiss_50_filtered.wPP.merged.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
$maxmiss_50_filtered.wPP.merged.vcf.gz

bcftools view -Oz -o $maxmiss_50_filtered.noPP.merged.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
$maxmiss_50_filtered.noPP.merged.vcf.gz

tabix $maxmiss_50_filtered.wPP.merged.nomono.vcf.gz
tabix $maxmiss_50_filtered.noPP.merged.nomono.vcf.gz

bcftools view -Oz -o $maxmiss_80_filtered.wPP.merged.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
$maxmiss_80_filtered.wPP.merged.vcf.gz

bcftools view -Oz -o $maxmiss_80_filtered.noPP.merged.nomono.vcf.gz -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
$maxmiss_80_filtered.noPP.merged.vcf.gz

tabix $maxmiss_80_filtered.wPP.merged.nomono.vcf.gz
tabix $maxmiss_80_filtered.noPP.merged.nomono.vcf.gz

# ## Filter for minor allele frequency
bcftools view -Oz -o $maxmiss_50_filtered.wPP.merged.nomono.minmaf0.03.vcf.gz -q 0.03:minor $maxmiss_50_filtered.wPP.merged.nomono.vcf.gz
bcftools view -Oz -o $maxmiss_50_filtered.noPP.merged.nomono.minmaf0.03.vcf.gz -q 0.03:minor $maxmiss_50_filtered.noPP.merged.nomono.vcf.gz

tabix $maxmiss_50_filtered.wPP.merged.nomono.minmaf0.03.vcf.gz
tabix $maxmiss_50_filtered.noPP.merged.nomono.minmaf0.03.vcf.gz

bcftools view -Oz -o $maxmiss_80_filtered.wPP.merged.nomono.minmaf0.03.vcf.gz -q 0.03:minor $maxmiss_80_filtered.wPP.merged.nomono.vcf.gz
bcftools view -Oz -o $maxmiss_80_filtered.noPP.merged.nomono.minmaf0.03.vcf.gz -q 0.03:minor $maxmiss_80_filtered.noPP.merged.nomono.vcf.gz

tabix $maxmiss_80_filtered.wPP.merged.nomono.minmaf0.03.vcf.gz
tabix $maxmiss_80_filtered.noPP.merged.nomono.minmaf0.03.vcf.gz

## Also filter by minimum allele count = 2 to remove singletons but keep doubletons, etc. as Linck & Battey 2019 show that this is best for pop structure analyses
vcftools --gzvcf $maxmiss_50_filtered.wPP.merged.nomono.vcf.gz --mac 2 --remove-filtered-all --recode --stdout | bgzip > $maxmiss_50_filtered.wPP.merged.nomono.minmac2.vcf.gz
vcftools --gzvcf $maxmiss_50_filtered.noPP.merged.nomono.vcf.gz --mac 2 --remove-filtered-all --recode --stdout | bgzip > $maxmiss_50_filtered.noPP.merged.nomono.minmac2.vcf.gz
vcftools --gzvcf $maxmiss_80_filtered.wPP.merged.nomono.vcf.gz --mac 2 --remove-filtered-all --recode --stdout | bgzip > $maxmiss_80_filtered.wPP.merged.nomono.minmac2.vcf.gz
vcftools --gzvcf $maxmiss_80_filtered.noPP.merged.nomono.vcf.gz --mac 2 --remove-filtered-all --recode --stdout | bgzip > $maxmiss_80_filtered.noPP.merged.nomono.minmac2.vcf.gz

tabix $maxmiss_50_filtered.wPP.merged.nomono.minmac2.vcf.gz
tabix $maxmiss_50_filtered.noPP.merged.nomono.minmac2.vcf.gz
tabix $maxmiss_80_filtered.wPP.merged.nomono.minmac2.vcf.gz
tabix $maxmiss_80_filtered.noPP.merged.nomono.minmac2.vcf.gz

## Generate a dataset with no missing data for some of the analyses
vcftools --gzvcf $maxmiss_50_filtered.wPP.merged.nomono.minmac2.vcf.gz --max-missing 1.0 --remove-filtered-all --recode --stdout | bgzip > $maxmiss_100_filtered.wPP.merged.nomono.minmac2.vcf.gz
vcftools --gzvcf $maxmiss_50_filtered.noPP.merged.nomono.minmac2.vcf.gz --max-missing 1.0 --remove-filtered-all --recode --stdout | bgzip > $maxmiss_100_filtered.noPP.merged.nomono.minmac2.vcf.gz

tabix $maxmiss_100_filtered.wPP.merged.nomono.minmac2.vcf.gz
tabix $maxmiss_100_filtered.noPP.merged.nomono.minmac2.vcf.gz
