#!/bin/bash

(./build/waterKernel \
    --gui false \
    --sizeZ 20.0 \
    --posZ -20.0 \
    --recordEventLog false \
    --nParticles 32 ) 2>&1 | tee myOutput.txt

rm *.rndm