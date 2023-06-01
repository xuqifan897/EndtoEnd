#!/bin/bash

# resultFolder="/data/qifan/projects/EndtoEnd4/results/slab6MeV1e7"
# if [ ! -d ${resultFolder} ]
# then
#     mkdir ${resultFolder}
# fi

# (time ./build/slabIPB \
#     --nParticles 32 \
#     --resultFolder ${resultFolder} \
#     --recordEventLog true \
#     ) 2>&1 | tee "${resultFolder}/myOutput.txt"

resultFolder="/home/qifan/projects/EndtoEnd4/results/slab6MeV1e7"
if [ ! -d ${resultFolder} ]
then
    mkdir ${resultFolder}
fi

(time ./build/slabIPB \
    --gui false \
    --nParticles 32 \
    --resultFolder ${resultFolder} \
    --recordEventLog true \
    ) 2>&1 | tee "${resultFolder}/myOutput.txt"