#!/bin/bash

(./build/waterKernel \
    --gui false \
    --sizeZ 20.0 \
    --posZ -20.0 \
    --recordEventLog true \
    --nParticles 32 \
    --resultFolder /home/qifan/projects/EndtoEnd4/IPB6MeV \
    ) 2>&1 | tee myOutput.txt

rm *.rndm