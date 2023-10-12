#!/bin/bash

if true ; then
    resultFolder="/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_1.0_250_1e8"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi
    for iteration in {0..15} ; do
        logFile="${resultFolder}/log_$( expr ${iteration} + 1 ).txt"
        echo $logFile
        ( time ./build/boxScore \
            --nParticles 100000000 \
            --dimXY 63 \
            --SAD 250 \
            --beamlet-size 0.5 \
            --resultFolder ${resultFolder} \
            --iteration ${iteration} ) 2>&1 > ${logFile}
    done
fi