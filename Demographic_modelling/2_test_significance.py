#!/usr/bin/python

import pickle
import nlopt
import matplotlib.pyplot as plt
import dadi
import random
import os

# 1. Load bootstraped frequency spectrum that has been generated before

import glob
boots_fids = glob.glob('Pmed.boots_*.fs')
boots_syn = [dadi.Spectrum.from_file(fid) for fid in boots_fids]

# 2.Calculate Godambe uncertainties

# Start a file to contain the confidence intervals
popt = [1.2553438154609, 1.11092969299215, 0.356498010726791, 0.116325045782146, 8.50029210425386, 14.2445575192057]         #Best parameter estimates obtained from you models
fi = open('Pmed.demographic_confidence_intervals.txt','w')
fi.write('Optimized parameters: {0}\n\n'.format(popt))

# We want to try a few different step sizes (eps) to see if uncertainties very wildly with changes to step size.

data_fs=data_fs = dadi.Spectrum.from_file('./Puffinus.SFS.dadi.fs')
demo_model = dadi.Demographics2D.split_delay_mig
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
ns = data_fs.sample_sizes
pts_l=[max(ns)+20, max(ns)+30, max(ns)+40]

for eps in [0.01, 0.001, 0.0001]:
    uncerts_adj = dadi.Godambe.GIM_uncert(demo_model_ex, pts_l, boots_syn, popt, data_fs, eps=eps)
    fi.write('Estimated 95% uncerts (with step size '+str(eps)+'): {0}\n'.format(1.96*uncerts_adj[:-1]))
    fi.write('Lower bounds of 95% confidence interval : {0}\n'.format(popt-1.96*uncerts_adj[:-1]))
    fi.write('Upper bounds of 95% confidence interval : {0}\n\n'.format(popt+1.96*uncerts_adj[:-1]))
fi.close()

# 3. Likelihood ratio test with other model

eps=[0.01]
nested_indices=[3]
adj = dadi.Godambe.LRT_adjust(demo_model_ex, pts_l, boots_syn, popt, data_fs, nested_indices, eps)

# With the output, you can obtain your D and weights parameters, and then calculate the p-value

D= # obtained using adj with the formula given in the dadi documentation
weights= # define here
dadi.Godambe.sum_chi2_ppf(D,weights)
