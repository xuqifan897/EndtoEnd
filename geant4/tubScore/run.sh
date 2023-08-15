#!/bin/bash

resultFolder="/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14"

if [ ! -d ${resultFolder} ]; then
    mkdir ${resultFolder}
fi

logFile="${resultFolder}/log.txt"

( time ./build/tubScore \
    --resultFolder /data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14 \
    --nParticles 10000000 ) 2>&1 > ${logFile} &