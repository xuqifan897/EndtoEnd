#!/bin/bash

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

(time ./build/slabAccu \
    --gui false \
    --nParticles 10000000 \
    --resultFolder ${experimentFolder} \
    ) 2>&1 | tee "${experimentFolder}/myOutput.txt"