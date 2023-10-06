#!/bin/bash

resultFolder="/data/qifan/projects/EndtoEnd/results/slabBench/patient1_g4"
if [ ! -d ${resultFolder} ]; then
    mkdir ${resultFolder}
fi
logFile="${resultFolder}/log.txt"

( time ./build/boxScore \
    --nParticles 10000000 \
    --SegZ 16 \
    --resultFolder ${resultFolder}
) 2>&1 > ${logFile} &