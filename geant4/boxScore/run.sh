#!/bin/bash

if true ; then
    resultFolder="/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi
    for iteration in {1..15} ; do
        logFile="${resultFolder}/log_$( expr ${iteration} + 1 ).txt"
        echo $logFile
        ( time ./build/boxScore \
            --nParticles 10000000 \
            --resultFolder ${resultFolder} \
            --iteration ${iteration} ) 2>&1 > ${logFile}
    done
fi

# if true; then
#     resultFolder="/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly"
#     if [ ! -d ${resultFolder} ]; then
#         mkdir ${resultFolder}
#     fi
#     logFile="${resultFolder}/log.txt"
#     ( time ./build/boxScore \
#         --gui true \
#         --dummy true \
#         --nParticles 10000000 \
#         --resultFolder ${resultFolder} ) > ${logFile}
# fi