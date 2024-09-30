This folder includes the scripts used to perform the SNP calling using an intersect of a) GATK haplotype caller (scripts 1-10) and b) Freebayes (scripts 10-12). After this scripts, the scripts on the following folders were run in succession:

  **· Autosomes_vs_chrZ**: scripts to discriminate between scaffolds belonging to autosomes or chromosome Z using a combination of: a) sex-specific coverage and b) synteny with related species.
  
  **· Variant_filtering**: scripts for filtering the VCFs produced with the previous scripts according to multiple variables (minor allele count, missingness, base quality and coverage, ...)
