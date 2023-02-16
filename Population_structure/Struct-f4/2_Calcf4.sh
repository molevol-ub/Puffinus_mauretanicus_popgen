#!/bin/bash

# Script to run Calc-f4; this script has to be placed inside the same directory as the input files
# Previously collect all input filenames into a filelist (with the format "path/filename\n")

# Options: -i "STRING" --> filelist name
#          -n "INT" --> nº of individuals included
#          -t "INT" --> nº of threads to run it on

Calc-f4 -i pmau_all.filelist -n 32 -t 6 | gzip -f > pmau.100.minmac2.calcf4.f4.gz
