# snpEff

#--------------------------------------------
#------------CREATE--A--DATABASE-------------
#--------------------------------------------

# Steps to create a database for your new genome in snpEff, obtained from: https://pcingola.github.io/SnpEff/se_build_db/

# 1. Configure a new genome

nano snpEff.config

# And write the following lines MANUALLY (e.g. with nano):

# CUSTOM GENOMES

# Balearic shearwater genome, v.1
pmau.genome : Puffinus_mauretanicus

#----------------------------------------------------------------

# 2. Use AGAT to correctly format de gff3 into a gtf (see documentation on the installation process in: https://github.com/NBISweden/AGAT)

agat_convert_sp_gff2gtf.pl --gff genes.gff3 -o genes.gtf

# More info on its use in: https://agat.readthedocs.io/en/latest/gff_to_gtf.html

#----------------------------------------------------------------

# 3. Build a database using gene annotation from a GTF file, a fasta AA alignment and a reference sequence

# Create directory for this new genome (same name as specified in config file (except the ".genome"!))

mkdir ~/Documents/software/snpEff/data/pmau
cd ~/Documents/software/snpEff/data/pmau

# Write the correct names

mv genome.fa sequences.fa
mv genome.fa.fai sequences.fa.fai

mv ../../genes.gtf .						# GTF annotation
mv ../../genes.aa.fasta protein.fa 		# FASTA AA allignment

cd ..

# Create the database - consider the correct names for the genome (name in config file except the ".genome"); gtf (genes.gtf) and AA allignment (do not check CDS, redundant?)

java -jar ~/Documents/software/snpEff/snpEff.jar build -gtf22 -v pmau -noCheckCds

# You have to be in the same directory as snpEff!

#--------------------------------------------
#----------SNP--EFFECT--PREDICTION-----------
#--------------------------------------------

# 4. Finally, if the database was correctly created, you should be able to annotate; more info in: https://www.youtube.com/watch?v=-rmreyRAbkE

# Take into account that VCF files should: 	a) Include monomorphic alternative sites! (Or else homozygous variants will not be detected)
#											b) No missing data! 

mkdir ~/Documents/software/snpEff_output
list=(ALT78 CZA11 ILA13 ILA2 PORQ TZE1 M8 M1 M5 M20 M4 G3 M2 M3 M18 M13 M19 M14 M12 M11 M17 M10 M21 G12 G9 G10 G11 G14 G15 M16 G4 M6 COP1 LT2)

for ind in ${list[*]}
do

#bgzip vcfs/$ind.masked.maxmiss1.vcf
#tabix vcfs/$ind.masked.maxmiss1.vcf.gz

#bcftools view -Oz -o vcfs/$ind.masked.maxmiss1.nomono.vcf.gz -e 'COUNT(GT="RR")=N_SAMPLES' vcfs/$ind.masked.maxmiss1.vcf.gz
#tabix vcfs/$ind.masked.maxmiss1.nomono.vcf.gz

gunzip vcfs/$ind.masked.maxmiss1.vcf.gz

java -jar ~/Documents/software/snpEff/snpEff.jar pmau vcfs/$ind.masked.maxmiss1.vcf > vcfs/$ind.masked.maxmiss1.ann.vcf

mv snpEff_genes.txt snpEff_output/snpEff_genes.$ind.txt 
mv snpEff_summary.html snpEff_output/snpEff_summary.$ind.html

done

#--------------------------------------------
#------------RESULT--MODIFICATION------------
#--------------------------------------------

# 5. Combine all vcfs

cp vcfs/ALT78.masked.maxmiss1.ann.vcf vcfs/all_inds.ann.vcf
bgzip vcfs/all_inds.ann.vcf
tabix vcfs/all_inds.ann.vcf.gz

list=(CZA11 ILA13 ILA2 PORQ TZE1 M8 M1 M5 M20 M4 G3 M2 M3 M18 M13 M19 M14 M12 M11 M17 M10 M21 G12 G9 G10 G11 G14 G15 M16 G4 M6 COP1 LT2)

for ind in ${list[*]}
do

bgzip vcfs/$ind.masked.maxmiss1.ann.vcf
tabix -f vcfs/$ind.masked.maxmiss1.ann.vcf.gz

bcftools merge -Oz -o vcfs/prova.vcf.gz --threads 2 vcfs/all_inds.ann.vcf.gz vcfs/$ind.masked.maxmiss1.ann.vcf.gz
mv vcfs/prova.vcf.gz vcfs/all_inds.ann.vcf.gz
tabix -f vcfs/all_inds.ann.vcf.gz

done

#---------------------------------------------------------------

cd vcfs

# 6. Change polarization of those genotypes that: a) Are monomorphic for the alternate allele in P.puffinus and b) Are polymorphic (mac=2 to exclude singletonces - sequencing errors) in Mediterranean Puffinus

# We have obtained the filtered dataset previously with bcftools.... Now:

cat overlap_noPP_PPmono.bed | cut -f 1,2 > polarise_positions.bed

# 6a. Create vcf with only those SNPs present in the filtered dataset

vcftools --gzvcf all_inds.ann.vcf.gz --positions polarise_positions.bed --recode-INFO-all --out included
mv included.recode.vcf included.vcf

# 6b. Change homozygous genotypes

cat included.vcf | sed 's/0\/0/1\/1/g' > prova
mv prova included.vcf
cat included.vcf | sed 's/0|0/1\/1/g' > prova
mv prova included.vcf
cat included.vcf | sed 's/1\/1/0\/0/g' > prova
mv prova included.vcf
cat included.vcf | sed 's/1|1/0\/0/g' > prova
mv prova included.vcf

# 6c. Create vcf with those SNPs excluded in the filtered dataset

vcftools --gzvcf all_inds.ann.vcf.gz --exclude-positions polarise_positions.bed --recode-INFO-all --out excluded
mv excluded.recode.vcf excluded.vcf

# 6d. Merge excluded vcf + included&modified vcf

bgzip included.vcf
bgzip excluded.vcf
tabix included.vcf.gz
tabix excluded.vcf.gz

bcftools concat -Oz -o all_inds.ann.def.vcf.gz -a --threads 12 included.vcf.gz excluded.vcf.gz
tabix all_inds.ann.def.vcf.gz

#---------------------------------------------------------------

# 7. Due to the difficulty of working with the LOF tags, we extract those positions with LOF tags into another file (first incorporating the header)

zcat all_inds.ann.def.vcf.gz | grep -v "##" | grep '#' > all_inds.ann.LOF.txt
zcat all_inds.ann.def.vcf.gz | grep 'LOF' >> all_inds.ann.LOF.txt

# We'll also create another file without the MODIFIER variants to speed things up

zcat all_inds.ann.def.vcf.gz  |  grep -v "##" | grep -v 'MODIFIER' > all_inds.ann.no_MODIFIER.txt

# 8. "Reduce" the results; for example, if you want to keep the variant, the genotype and their effects:

file_list=(no_MODIFIER LOF)

for file in ${file_list[*]}
do

cat all_inds.ann.$file.txt | grep -v "##" | cut --complement -f 3,4,5,6,7,9 > prova
mv prova all_inds.ann.$file.txt

done

# 9. Now change positions of the columns so you can easily divide the "INFO" column into more columns, of which you'll keep just the 3 interesting ones

file_list=(no_MODIFIER LOF)

for file in ${file_list[*]}
do

# Divide your file to extract the INFO column

cat all_inds.ann.$file.txt | cut -f 1,2 > pos.txt
cat all_inds.ann.$file.txt | cut -f 3 > info.txt
cat all_inds.ann.$file.txt | cut --complement -f 1,2,3 > genotype.txt

# Keep only interesting parts of INFO tag

cat info.txt | sed 's/|/\t/g' | cut -f 2,3,5 > prova
mv prova info.txt

# Keep only the 3 first characters of the genotype tag (previously extracting the header so that it isn't "snipped") - and (previously) change phased genotypes to unphased

cat genotype.txt | sed 's/|/\//g' > prova
mv prova genotype.txt

head -n1 genotype.txt > prova
cat genotype.txt | grep '/' | awk '{for(i=1;i<=NF; i++) {$i=substr($i,1,3) }  OFS=" "; print $0} ' >> prova
mv prova genotype.txt

# Paste all together

paste pos.txt genotype.txt info.txt > outfile.$file.txt

done

#---------------------------------------------------

# And now you are ready to summarise your results, for example with R
# Previously modify the header to better fit your data (remove th "#" from "#CHROM" for example)
# And just in case there are any spaces that should be tabs:

file_list=(no_MODIFIER LOF)

for file in ${file_list[*]}
do

cat outfile.$file.txt | sed 's/\s/\t/g' > prova
mv prova outfile.$file.txt

done
