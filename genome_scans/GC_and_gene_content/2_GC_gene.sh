#!/bin/bash
#$ -cwd
#$ -R y
#$ -e GC_bedtools.err
#$ -o GC_bedtools.out
#$ -q h12.q
#$ -pe ompi128h12 2
#$ -V                    #export environment var
#$ -N bedtools_nuc             #name Job

# Bash script to calculate GC content using bedtools

window_size=25

WORK_DIR=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/Z/$window_size\_windows
genome_file=$WORK_DIR/../../../masking/genome_masked.def.fa
stat_file=$WORK_DIR/auto_statistics.unmasked.csv
window_file=$WORK_DIR/$window_size\_windows.bed              # Watch out -correct?
#window_file=$WORK_DIR/$window_size\_windows.v2.bed		# Watch out
output_file=$WORK_DIR/gc_output.Z.txt				# Check if Z or autososmes


# 1 Use genome_masked.fa to calculate %GC content for all windows through bedtools nuc; first generate windows bedfile again (to be sure)

cat $stat_file | cut -f 9,10,11 > $window_file

# Remove the "" from the scaffold names and remove 1st line

sed 's|["]||g' $window_file > $window_file.def.bed
mv $window_file.def.bed $window_file

sed '1d' $window_file > $window_file.def.bed
mv $window_file.def.bed $window_file

bedtools nuc -fi $genome_file -bed $window_file > $output_file

# 2 Add it to the autosome_file

cat $output_file | cut -f 5 > $WORK_DIR/gc_prov.csv

sed 's/5_pct_gc/GC_content/' $WORK_DIR/gc_prov.csv > $WORK_DIR/gc_prov_def.csv
mv $WORK_DIR/gc_prov_def.csv $WORK_DIR/gc_prov.csv

paste -d "\t" $stat_file $WORK_DIR/gc_prov.csv > $WORK_DIR/auto_statistics.gc.unmasked.csv
rm $WORK_DIR/gc_prov.csv

# ---------------------------------------------------------------------------------------

# Now calculate the gene density

mask_file=$WORK_DIR/../../../masking/final_masked.redux.bed
annot_file=/users-d3/ccuevas/Functional_Genome_Annotation/scripts/braker3/3pmaureta/augustus.hints.gff3
gene_bedfile=$WORK_DIR/gene_list.Z.bed
overlap_file=$WORK_DIR/gene_overlap.Z.csv
final_output=$WORK_DIR/auto_statistics.gc.unmasked.csv

# 1 Retain only lines with genes and columns with scf_name and pos

grep "gene" $annot_file > $gene_bedfile

cat $gene_bedfile | cut -f 1,4,5 > $gene_bedfile.def.bed
mv $gene_bedfile.def.bed $gene_bedfile

# 2 Run bedtools intersect to mask gene_bedfile with the results of the masking

cat $mask_file | cut -f 1,2,3 >  $WORK_DIR/temp

bedtools intersect -a $gene_bedfile -b $WORK_DIR/temp -bed > $gene_bedfile.def.bed
mv $gene_bedfile.def.bed $gene_bedfile

rm $WORK_DIR/temp

# 3 Run bedtools intersect with the window_file to calculate the amount of overlap (-wao)

bedtools intersect -a $window_file -b $gene_bedfile -wao > $overlap_file

# 4 Run custom python script to calculate overlap percentage

python3 gene_content.py $overlap_file $window_size

# 5 Add to statistics file

sed -i '1s/^/gene_density\n/' $overlap_file

paste -d "\t" $final_output $overlap_file > $WORK_DIR/temp
mv $WORK_DIR/temp $final_output

# 6 Modify the gene_content to take into account the unmasked proportion using custom script, but first calculate the overlapping features again, just in case there are errors

cat $mask_file | cut -f 1,2,3 >  $WORK_DIR/temp

bedtools intersect -a $window_file -b $WORK_DIR/temp -wao > $WORK_DIR/overlap_amount.Z.csv

# Reuse your first script

python3 gene_content.py $WORK_DIR/overlap_amount.Z.csv $window_size

sed -i '1s/^/Unmasked_prop_v2\n/' $WORK_DIR/overlap_amount.Z.csv

paste -d "\t" $final_output $WORK_DIR/overlap_amount.Z.csv > $WORK_DIR/temp
mv $WORK_DIR/temp $final_output
rm $WORK_DIR/overlap_amount.Z.csv

# Now use the new one

python3 correct_denominator.py $final_output $window_size

#------------------------------------------------------------------------------------

# Finally remove the unused header (unmasked_proportion v2)

cat $final_output | cut -f 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 > $final_output.prov
mv $final_output.prov $final_output
