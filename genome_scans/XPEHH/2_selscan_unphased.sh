#!/bin/bash
#$ -cwd
#$ -R y
#$ -e selscan.err
#$ -o selscan.out
#$ -q h12.q
#$ -pe ompi128h12 12
#$ -V                    #export environment var
#$ -N selscan             #name Job

# Script to create mapfile - later we can run selscan

WORKDIR=/users-d3/jferrer/gizquierdo/TFG/genome_scans/selscan/input
outdir=/users-d3/jferrer/gizquierdo/TFG/genome_scans/selscan/output

# Options:	--unphased (if data unphased)
#		--vcf (input_file)
#		--map (accompanying mapfile)
#		--xpehh (to calculate XP-EHH statistic)
#		--out (outfile prefix)
#		--ehh-win (length of window from side to side of SNP to calculate EHH) - for the moment let's try both default (100kbp) and 12.5kbp
#		--pmap (use phisical map instead of genetical) --> should we?
#		--threads
#		--trunc-ok (should we?)

# Do this for every scaffold and join all ".out" files in one (without header, included later)

#rm $outdir/Puffinus_selscan*			# Remove any previous outfiles just in case

touch $outdir/Puffinus_selscan_Z.out 

while read scaff;
	do
		infile=$WORKDIR/$scaff.mauretanicus.vcf;
		mapfile=$WORKDIR/$scaff.selscan.map;
		reffile=$WORKDIR/$scaff.yelkouan.vcf;
		selscan --xpehh --unphased --vcf $infile --vcf-ref $reffile --map $mapfile --out $outdir/$scaff.win100 --threads 12;
		
		cat $outdir/$scaff.win100.xpehh.out | sed '1d' >> $outdir/Puffinus_selscan_Z.out
		
	done < $WORKDIR/scaff_names

# Falta ficar header

rm $WORKDIR/scf*
rm $outdir/scf*
