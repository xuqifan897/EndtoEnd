#!/bin/bash

resultRoot="/data/qifan/projects/EndtoEnd/results/Sept1Point"
if [ ! -d ${resultRoot} ]; then
    mkdir ${resultRoot}
fi

# Native dose calculation script
if false; then
    resultFolder="${resultRoot}/largerMat"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi
    logFile="${resultFolder}/log.txt"

    ( time ./build/waterPoint \
        --nParticles 10000000 \
        --PhantomDimXY 50 \
        --PhantomDimZ 100 \
        --PhantomSZ 5 \
        --logFrequency 10000 \
        --resultFolder ${resultFolder} ) \
        > ${logFile} 2>&1 &
fi

# Test trajectory functionality
if true; then
    resultFolder="${resultRoot}/TrjPrint"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi

    logFile="${resultFolder}/log.txt"
    ( time ./build/waterPoint \
        --nParticles 100 \
        --PrintTrajectory true \
        --PhantomDimXY 50 \
        --PhantomDimZ 100 \
        --PhantomSZ 5 \
        --resultFolder ${resultFolder} ) > ${logFile} 2>&1 & 
fi

if false; then
    resultFolder="${resultRoot}/TrjPrint"
    if [ ! -d ${resultFolder} ]; then
        mkdir ${resultFolder}
    fi

    gdb --args ./build/waterPoint \
        --nParticles 100 \
        --PrintTrajectory true \
        --PhantomDimXY 50 \
        --PhantomDimZ 100 \
        --PhantomSZ 5 \
        --resultFolder ${resultFolder}
fi