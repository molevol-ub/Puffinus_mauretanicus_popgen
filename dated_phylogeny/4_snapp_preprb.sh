# Bash script to run snapp_prep.rb in order to prepare files for time-callibrated SNAPP analyses

# Options:  -v [string] --> input vcf file
#           -c [string] --> constraints file for the tree callibration --> see an example within this directory
#           -t [string] --> file that relates populations and individuals
#           -m [int] --> nº of iterations between logs
#           -l [int] --> MCMC chain length
#           -x, -o [string] --> name of the ouput

ruby /path/to/snapp_prep.rb -v scfE_downsized.maxmiss100.vcf -c constraints_taxa.normal.txt -t Puffinus_subset.csv \
-m 1000 -l 500000 -x scfE_downsized.maxmiss100.normal.xml -o scfE_downsized.maxmiss100.normal
