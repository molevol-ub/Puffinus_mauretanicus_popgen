# Population genomics in _Puffinus mauretanicus_

# Introduction

This git repository contains those scripts used during the population genomics analyses performed with _Pufifnus mauretanicus_ by the Molecular Evolutionary Genetics Group, University of Barcelona.


# Structure

1. **SNP_calling**: SNP calling pipeline for WGS resequencing data that uses a combination of freebayes and GATK HaplotypeCaller.
2. **SNP_calling_mitogenome**: SNP calling pipeline which uses freebayes - adapted for mitogenomic data.
3. **Population_structure**: includes those analyses performed to visualize population structure within *P. mauretanicus* and *P. yelkouan*. Among these, we include: PCAs, Struct-f4, TreeMix and haplotype genealogy graphs for mitogenomes (using Fitchi).
4. **Dated_phylogeny**: inference of time-callibrated species trees using SNAPP.
5. **Demographic_modelling**: scripts for demographic modelling using dadi, alongside custom models.
6. **Genome_scans**: performance of genome-wide scans for signatures of selection, including detailed analyses of the main outlier windows.
7. **Conservation_genomics**: scripts to calculate heteroziogsity, LROHs, ...
8. **Forward_simulations**: scripts to use SLiM simulations to investigate the plausibility of mito-nuclear discordance in our model and to assess the effects of migration/demography in the species heterozygosity/FROH.
9. **Others**: smaller tasks, such as finding fixed differences between populations, calculating the SFS, ...
