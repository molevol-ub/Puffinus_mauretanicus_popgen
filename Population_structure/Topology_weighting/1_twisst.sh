#!/bin/bash
#$ -cwd
#$ -R y
#$ -e twisst.err
#$ -o twisst.out
#$ -q h14.q
#$ -pe ompi511h14 24
#$ -V                    #export environment var
#$ -N twisst             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

python parseVCF.py -i ~/pmau_popgen/SNP_calling/vcfs/phased_vcfs_def/autosomes/Puffinus_def.auto.shapeit4_whatshap_phased.wPP.maxmiss80.nomono.vcf | gzip > output.geno.gz
python phyml_sliding_windows.py -g output.geno.gz -T 24 --minPerInd 1 --prefix Puffinus.phased.wPP -w 100 --windType sites --model GTR --optimise n
source activate py39
python twisst.py -t Puffinus.phased.wPP.trees.gz -w Puffinus.phased.4pops.weights.csv.gz -g Ppuf -g Pyel -g Pit -g Mal --outgroup Ppuf --groupsFile Puffinus.list
python twisst.py -t Puffinus.phased.wPP.trees.gz -w Puffinus.phased.noCab.weights.csv.gz -g Ppuf -g Pyel -g Men -g Pit -g Mal --outgroup Ppuf --groupsFile Puffinus.list      #including Men
python twisst.py -t Puffinus.phased.wPP.trees.gz -w Puffinus.phased.noMen.weights.csv.gz -g Ppuf -g Pyel -g Cab -g Pit -g Mal --outgroup Ppuf --groupsFile Puffinus.list      #including Cab
