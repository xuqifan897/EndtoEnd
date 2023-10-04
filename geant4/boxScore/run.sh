#!/bin/bash

if true ; then
./build/boxScore \
    --nParticles 10000000 \
    --scoring false \
    --SegZ 20 \
    --resultFolder /data/qifan/projects/EndtoEnd/results/slabBench/patient1_g4
fi

if false ; then
gdb --args ./build/boxScore \
    --nParticles 10000000 \
    --SegZ 20 \
    --resultFolder /data/qifan/projects/EndtoEnd/results/slabBench/patient1_g4
fi

# break /data/qifan/projects/EndtoEnd/geant4/boxScore/src/Run.cpp:46