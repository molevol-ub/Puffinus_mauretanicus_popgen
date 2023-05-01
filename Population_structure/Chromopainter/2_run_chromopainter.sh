#!/bin/bash
#$ -cwd
#$ -R y
#$ -e gatk_reference.err
#$ -o gatk_reference.out
#$ -q h12.q
#$ -pe ompi128h12 4
#$ -V                    #export environment var
#$ -N gatk_reference             #name Job
#$ -M 000izquierdoguillem@gmail.com
#$ -m be

cd /users-d3/jferrer/gizquierdo/TFM/chromopainter

# 1. Estimate the required parameters fro chromopainter using a SUBSAMPLE OF INDIVIDUALS AND SITES; use options:
#      · “-s 0” to specify that no paintingsamples are needed in this step
#      · “-i 10 -in -iM” specifies to use 10 iterations of Expectation-Maximisation to infer the switch (“-in”) and emission (“-iM”) parameters. 

ChromoPainterv2 -g Puffinus_subset.chromopainter.inp -r Puffinus_subset.recomrates.txt -a 0 0 -s 0 -i 10 -in -iM -o Puffinus_subset

# 2. Calculate the average estimated values for n and M to use in the final chromopainter run. The script might need some polishing to use the appropiate input (parameters "$infilePREFIX", "$infileSUFFIX", "@chromovec" and "@chromolengths")
# For more info see: https://github.com/hellenthal-group-UCL/fastGLOBETROTTER/blob/master/tutorial/ChromoPainterv2EstimatedNeMutExtractEM.pl

perl scripts/ChromoPainterv2EstimatedNeMutExtractEM.pl 

# In our case, the resulting estimates are: -n=33700.840; -M=0.022739

# 3. Using the output values of the previous step, run chromopainter

ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt -a 0 0 -n 33700.840 -M 0.022739 -o Puffinus.auto.nPP.maxmiss80

#------------------------------------------------------------------------------------------------------------------------------------

# 4. Prepare the input files for GLOBETROTTER

# 4.1 Generate the copy vector input file

ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f Puffinus.auto.nPP.maxmiss80.popfile.txt 0 0 -s 0 -n 33700.840 -M 0.022739 -o Puffinus.auto.nPP.maxmiss80

# Then paint: the different iterations correspond to different "target" (coinciding with the directory name) and source populations

mkdir Pyel
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_pyel.txt 0 0 -n 33700.840 -M 0.022739 -o ./Pyel/Puffinus.auto.nPP.maxmiss80
mkdir Menorca
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_men.txt 0 0 -n 33700.840 -M 0.022739 -o ./Menorca/Puffinus.auto.nPP.maxmiss80
mkdir Menorca_noMen
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_men2.txt 0 0 -n 33700.840 -M 0.022739 -o ./Menorca_noMen/Puffinus.auto.nPP.maxmiss80
mkdir Menorca_noMenCab
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_men3.txt 0 0 -n 33700.840 -M 0.022739 -o ./Menorca_noMenCab/Puffinus.auto.nPP.maxmiss80
mkdir Cabrera
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_cab.txt 0 0 -n 33700.840 -M 0.022739 -o ./Cabrera/Puffinus.auto.nPP.maxmiss80
mkdir Cabrera_noCab
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_cab2.txt 0 0 -n 33700.840 -M 0.022739 -o ./Cabrera_noCab/Puffinus.auto.nPP.maxmiss80
mkdir Cabrera_noMenCab
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_cab3.txt 0 0 -n 33700.840 -M 0.022739 -o ./Cabrera_noMenCab/Puffinus.auto.nPP.maxmiss80
mkdir Mallorca
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_mal.txt 0 0 -n 33700.840 -M 0.022739 -o ./Mallorca/Puffinus.auto.nPP.maxmiss80
mkdir Pitiuses
ChromoPainterv2 -g Puffinus.auto.nPP.maxmiss80.chromopainter.inp -r Puffinus.auto.nPP.maxmiss80.recomrates.txt \
  -t Puffinus.auto.nPP.maxmiss80.idfile.txt -f popfiles/Puffinus.auto.nPP.maxmiss80.popfile_pit.txt 0 0 -n 33700.840 -M 0.022739 -o ./Pitiuses/Puffinus.auto.nPP.maxmiss80
