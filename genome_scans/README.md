# Genome-wide scans for divergent selection between *Puffinus* taxa

The 9 **numbered** scripts (except script 5) in this directory are to be run in order to obtain a file that summarizes statistics for the specified window size.

  - Take into account that the 3 first ones aim to refine the masking of the genome; don't use if this has been performed in advance.

After these, 1st you should run the **GC_and_gene** folder and 2nd the **XPEHH** one. Once this is done, you are ready to visualize the results found in the final output file using the scripts found in the **Rplot** directory.

## Characterization of the main outlier window

Aisde from these, you will find 4 folders for different complementary analyses that may be useful to perform once you have identified **candidate scaffolds**:
  - **scf_location**: to detect synteny with chromosomes of closely-related species; useful for identifying neighboring scaffolds to your target ones.
  - **LD**: to generate a linkage disequilibrium heatmap.
  - **GenotypePlot**: to generate genotype plots, which help in visualizing distribution of long haplotypes across our individuals.
  - **PCA_windows**: to generate PCAs across the windows found within the scaffold.
  - **Genes**: to visualize which genes are found within our scaffolds.
  
 ## Scans for balancing selection
 
 Scripts to perform genome-wide scans aimed at detecting regions that might be under **balancing selection** in one of the taxa can be found in:
  - **BetaScan**: which calculate Beta statistics (in order to detect regions with SNPs with correlated MAFS).
  - **BalLerMix+**: which computes a composite likelihood ratio test in order to detects SNPs that might be under balancing selection according to the SFS around it.
 
 It must be stressed that both of these scans are more easily run once the 9 **numbered scripts** (and those in the GC_and_gene & XPEHH directories as well).
