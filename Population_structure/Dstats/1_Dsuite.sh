#!/bin/bash
#$ -cwd
#$ -R y
#$ -e Dsuite.err
#$ -o Dsuite.out
#$ -q h12.q
#$ -pe ompi128h12 4
#$ -V                    #export environment var
#$ -N Dsuite             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

WORKDIR=/users-d3/jferrer/gizquierdo/TFM/Dstats
cd $WORKDIR
VCF=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.wPP.merged.nomono.masked.vcf.gz
VCF_mac=/users-d3/jferrer/gizquierdo/TFG/chr_split/vcfs/auto/Puffinus_SNP.maxmiss80.filtered.noPP.merged.nomono.minmac2.masked.vcf.gz

# 1. Calculate D-stats and f-stats using Dtrios specifying an imput tree where Cabrera is the outgroup to all other Pmau

Dsuite Dtrios -n Puffinus.Cabrera_out -t Puffinus.Cabrera_out.nwk $VCF Puffinus_list.txt

echo "1st run completed"

# Repeat but with a tree where Cabrera is sister to Mallorca

Dsuite Dtrios -n Puffinus.Pit_out -t Puffinus.Pit_out.nwk $VCF Puffinus_list.txt

echo "2nd run completed"

cut -f 2 Puffinus_list.txt | uniq > plot_order.txt

# 2. Plot, limiting the highest value to the maximum D value and f4 value respectively

ruby plot_d.rb Puffinus_list_Puffinus.Cabrera_out_BBAA.txt plot_order.txt 0.0121143 Puffinus.Cab_out.BBAA.D.svg
ruby plot_f4ratio.rb Puffinus_list_Puffinus.Cabrera_out_BBAA.txt plot_order.txt 0.047654 Puffinus.Cab_out.BBAA.f4ratio.svg
