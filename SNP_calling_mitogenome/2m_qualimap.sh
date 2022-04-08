# Script to run qualimap (BAM QC) to evaluate coverage and other statistics of the masked bams

# Qualimap manual: http://qualimap.conesalab.org/doc_html/command_line.html

work_dir=/users-d3/jferrer/pmau_popgen/SNP_calling/mt_bams/clean_bams

ind_list=(ALT78 COP1 CZA11 G10 G11 G12 G14 G15 G3 G4 G9 ILA13 ILA2 LT2 M1 M10 M11 M12 M13 M14 M16 M17 M18 M19 M2 M20 M21 M22 M3 M4 M5 M6 M7 M8 PORQ TZE1 sacella)

for ind in "${ind_list[@]}"
	do
		qualimap bamqc -bam $work_dir/${ind}.sorted.dups.bam -outdir $work_dir/${ind}_qualimap_results -nt 2
	done
# Check effect of "skip duplicates" and "gff feature file"
