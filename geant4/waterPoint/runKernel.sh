#!/bin/bash

resultRoot="/data/qifan/projects/EndtoEnd/results/Sept1Point"
if [ ! -d ${resultRoot} ]; then
    mkdir ${resultRoot}
fi

resultFolder="${resultRoot}/pointKernel"
if [ ! -d ${resultFolder} ]; then
    mkdir ${resultFolder}
fi

if true; then
    logFile="${resultFolder}/log.txt"
    (time ./build/waterPoint \
        --nParticles 10000000 \
        --PhantomDimXY 50 \
        --PhantomDimZ 100 \
        --PhantomSZ 5 \
        --PhantomBottom 50 \
        --resultFolder ${resultFolder}) > ${logFile} 2>&1 &
fi

if false; then
    gdb --args ./build/waterPoint \
        --nParticles 10000000 \
        --PhantomDimXY 50 \
        --PhantomDimZ 100 \
        --PhantomSZ 5 \
        --PhantomBottom 50 \
        --resultFolder ${resultFolder}
fi

# break /data/qifan/projects/EndtoEnd/geant4/waterPoint/src/EventAction.cpp:177
# print *(char**)(&(trj->fParticleName))