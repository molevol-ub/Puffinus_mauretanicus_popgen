# Information about switchError is found in pull requests in: https://github.com/SPG-group/switchError

cd

# 1. Activate the conda enviroment were switchError is installed

source activate switchError

# 2. Run shapeit with all but 1 individual; just for trial we'll iterate through the removal of 5 individuals: CZA11, G15, M19, G11, M17

# 3. Generate input bcf files for switchError (1st shapeit, 2nd whatshap-phased)

vcfdir=/users-d3/jferrer/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes

for ind in (CZA11 G15 M19 G11 M17)
do
bcftools --view --Ou -o Puffinus_shapeit.$ind.bcf $vcfdir/Puffinus_raw.$ind.shapeit4_whatshap_phased.vcf.gz
done

bcftools --view --Ou -o Puffinus_whatshap.bcf $vcfdir/whatshap_phased/Puffinus_raw.auto.whatshap_phased.vcf.gz
