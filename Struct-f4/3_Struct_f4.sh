#!/bin/bash

# Rscript to run Struct-f4 according to https://bitbucket.org/plibradosanz/structf4/src/master/
# Beware, this script can only be used if Struct-f4 runs with R-4
# Send to either hercules11 or hercules 12

# Options: -f "CHAR" --> path to input file
#          -T "INT" --> number of threads
#          -K "INT" --> number of ancestral populations
#          -p "CHAR" --> path to plots (pdf) output

Struct-f4 -f /path/file -T 6 -K 2 -p /output-path > /path/output_summary.out

# Additionally you may choose to use: -m "INT" --> number of iterations for 1st MCMC chain
#                                     -n "INT" --> number of iterations for 2nd MCMC chain
# You can reduce them from the default settings (m=1M & n=5M) but first check their efficiency with the output plots with the default settings
