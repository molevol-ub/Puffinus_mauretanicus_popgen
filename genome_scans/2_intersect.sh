# Bash script to intersect the mapability map gff3 with the coverage mask

mapability_mask=/users-d3/jferrer/MSMC2_pmau/input/mapp_mask_repeats_subtract.bed
coverage_mask=/users-d3/jferrer/gizquierdo/TFG/genome_scans/masking/scaff_coverage.mask.bed.gz
output=/users-d3/jferrer/gizquierdo/TFG/genome_scans/masking/final_masked.bed

bedtools intersect -a $mapability_mask -b $coverage_mask -bed > $output
