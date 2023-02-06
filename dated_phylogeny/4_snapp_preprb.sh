# Rscript to run snapp_prep.rb in order to prepare files for time-callibrated SNAPP analyses

ruby snapp_prep-master/snapp_prep.rb -v scfE_downsized.maxmiss100.vcf -c constraints_taxa.normal.txt -t Puffinus_subset.csv -m 1000 -l 500000 -x scfE_downsized.maxmiss100.normal.xml -o scfE_downsized.maxmiss100.normal
