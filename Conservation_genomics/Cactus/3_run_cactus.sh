# We add our genomes to the 363-way allignment by solving 2 subproblems
# More info here

source activate py39
cd /users-d3/jferrer/gizquierdo/TFM/Cactus/Cactus

mkdir jobstore

cactus ./jobstore Puffinus_run1.txt Puffinus_run1.hal --root birdAnc361 
