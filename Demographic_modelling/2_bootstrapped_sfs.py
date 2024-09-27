#!/usr/bin/python

import pickle
import nlopt
import matplotlib.pyplot as plt
import dadi

# 1. Load VCF data

datafile = './regions/Puffinus.gene_masked.macro_masked.vcf.gz'
dd = dadi.Misc.make_data_dict_vcf(datafile, '../Puffinus.dadi.popfile.txt')

# Define pops
pop_ids, ns = ['Pyel', 'Pmau'], [24, 40]

# 2. Define number of bootstraps and chunk size
Nboot, chunk_size = 100, 1e7

# 3. Define the chunks and generate bootstrapped SFS
chunks = dadi.Misc.fragment_data_dict(dd, chunk_size)
import random
random.seed(12345)
boots = dadi.Misc.bootstraps_from_dd_chunks(chunks, Nboot, pop_ids, ns, polarized=False)
for i in range(len(boots)):
  boots[i].to_file('Pmed.boots_{0}.fs'.format(str(i)))

# 4. Save the boots
import os
pick = open('Puffinus.data_dict','wb')
pickle.dump(dd, pick, 2)

