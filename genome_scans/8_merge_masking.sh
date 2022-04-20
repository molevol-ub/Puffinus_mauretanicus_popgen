# Bash script to calculat the masking proportion of each window

window_size=1

auto_file=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/$window_size\_windows/autosome_statistics.csv
mask_file=/users-d3/jferrer/gizquierdo/TFG/genome_scans/masking/final_masked.redux.bed
window_file=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/$window_size\_windows/$window_size\_windows.bed
output_file=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/$window_size\_windows/masking_amount.bed
final_file=/users-d3/jferrer/gizquierdo/TFG/genome_scans/PopGenome/$window_size\_windows/masking_prop.csv

# 1 Generate bed file of the windows used (scf and start+stop)

cat $auto_file | cut -f 9,10,11 > $window_file

# 2 Remove the "" from the scaffold names and remove 1st line

sed 's|["]||g' $window_file > $window_file.def.bed
mv $window_file.def.bed $window_file

sed '1d' $window_file > $window_file.def.bed
mv $window_file.def.bed $window_file

# 3 Use bedtools intersect to calculate the number of overlapping positions between windows

bedtools intersect -a $window_file -b $mask_file -wo > $output_file

# 4 Keep only unique lines (problems arise in "border" areas)

cat -n $output_file | sort -uk2 | sort -n | cut -f2- > $output_file.def.bed
mv $output_file.def.bed $output_file

# 5 Run python script to calculate the percentage

python3 masking_percentage.py $window_size
