#---------------------------------------------------------------------
#---------------PREPARE  THE  CACTUS  363-WAY  ALLIGNMENT------------
#---------------------------------------------------------------------

# 1. Download the 363-way allignment

wget -P /media/guillem/BC90A1CD90A18F08/Guillem/TFM_GIA/Cactus https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020-hub/Gallus_gallus/Gallus_gallus.maf.gz
cd /your/dir

# 2. Export the library paths that are required in the following scripts; if you don't know the path, use "locate"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/disk3/h-user/jferrer/programari/anaconda3/bin/cactus/js/build-temp-027501980/rootfs/usr/lib/x86_64-linux-gnu/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/disk3/h-user/vadim.pisarenco/anaconda3/lib/

# 3. Convert HAL to FASTA using 


# Remember to run with Singularity instead of Docker

--binariesMode singularity
