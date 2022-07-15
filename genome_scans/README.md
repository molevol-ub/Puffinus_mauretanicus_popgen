The 9 **numbered** scripts (except script 5) in this directory are to be run in order to obtain a file that summarizes statistics for the specified window size.

After these, 1st you should run the **GC_and_gene** folder and 2nd the **XPEHH** one. Once this is done, you are ready to visualize the results found in the final output file using the scripts found in the **Rplot** directory.

Aisde from these, you will find 4 folders for different complementary analyses that may be useful to perform once you have identified **candidate scaffolds**:
  - **scf_location**: to detect synteny with chromosomes of closely-related species; useful for identifying neighboring scaffolds to your target ones.
  - **LD**: to generate a linkage disequilibrium heatmap.
  - **GenotypePlot**: to generate genotype plots, which help in visualizing distribution of long haplotypes across our individuals.
  - **PCA_windows**: to generate PCAs across the windows found within the scaffold.
  - **Genes**: to visualize which genes are found within our scaffolds.
