#!/bin/bash

# Run locally, hal2fasta does not work on hercules (missing libraries)

#---------------------------------------------------------------------
#---------------PREPARE  THE  CACTUS  363-WAY  ALLIGNMENT------------
#---------------------------------------------------------------------

# 1. Download the 363-way allignment

wget -P /media/guillem/BC90A1CD90A18F08/Guillem/TFM_GIA/Cactus https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020-hub/Gallus_gallus/Gallus_gallus.maf.gz
cd /your/dir

# If needed install hal: instructions in https://github.com/ComparativeGenomicsToolkit/hal
# 2. Export the library paths that are required in the following scripts; if you don't know the path, use "locate"

export PATH=<path to hal>/bin:${PATH}                    #/home/guillem/Documents/software/hal
export PYTHONPATH=<parent of hal>:${PYTHONPATH}          #/home/guillem/Documents/software

# 3. Convert HAL to FASTA using hal2fasta

bird363=/media/guillem/BC90A1CD90A18F08/Guillem/TFM_GIA/Cactus/363-avian-2020.hal
hal2fasta $bird363 $(halStats --root $bird363) --subtree --upper --ucscSequenceNames > /media/guillem/BC90A1CD90A18F08/Guillem/TFM_GIA/Cactus/363_bird.fa &

# 4. Remove ancestral bird genomes

# Remember to run with Singularity instead of Docker

--binariesMode singularity
