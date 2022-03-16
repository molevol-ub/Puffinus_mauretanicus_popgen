#!/bin/bash

# Script to run Calc-f4

# Take into account that previously you must have collected all input filenames into a filelist

# Run it in 6 threads and with 32 samples

Calc-f4 -i pmau_all.filelist -n 32 -t 6 | gzip -f > pmau.100.minmac2.calcf4.f4.gz
