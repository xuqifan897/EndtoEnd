#!/bin/sh

resultFolder="/data/qifan/projects/EndtoEnd4/results/InhomoJuly20"
if [ ! -d ${resultFolder} ]
then
    mkdir ${resultFolder}
fi

experimentFolder="${resultFolder}/slab2"
if [ ! -d ${experimentFolder} ]
then
    mkdir ${experimentFolder}
fi

gdb --args ./build/slabAccu \
    --gui false \
    --nParticles 10000000 \
    --resultFolder ${experimentFolder}