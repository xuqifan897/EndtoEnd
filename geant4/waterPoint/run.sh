#!/bin/bash

resultRoot="/data/qifan/projects/EndtoEnd/results/Sept1Point"
if [ ! -d ${resultRoot} ]; then
    mkdir ${resultRoot}
fi

resultFolder="${resultRoot}/initialDet"
if [ ! -d ${resultFolder} ]; then
    mkdir ${resultFolder}
fi
logFile="${resultFolder}/log.txt"

# ./build/waterPoint \
#     --nParticles 10000000 \
#     --resultFolder ${resultFolder}

( time ./build/waterPoint \
    --nParticles 10000000 \
    --PhantomDimXY 50 \
    --PhantomDimZ 100 \
    --PhantomSZ 5 \
    --logFrequency 10000 \
    --resultFolder ${resultFolder} ) \
    > ${logFile} 2>&1 &