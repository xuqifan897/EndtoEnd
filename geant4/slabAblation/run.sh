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

# resultFolder="/data/qifan/projects/EndtoEnd4/results/slab6MeVAblationVanilla"
# if [ ! -d ${resultFolder} ]
# then
#     mkdir ${resultFolder}
# fi

# (time ./build/slabIPB \
#     --gui false \
#     --nParticles 10000000 \
#     --resultFolder ${resultFolder} \
#     --recordEventLog false \
#     ) 2>&1 | tee "${resultFolder}/myOutput.txt"


resultFolder="/data/qifan/projects/EndtoEnd4/results/InhomoJuly6"
if [ ! -d ${resultFolder} ]
then
    mkdir ${resultFolder}
fi

experimentFolder="${resultFolder}/boneDose"
if [ ! -d ${experimentFolder} ]
then
    mkdir ${experimentFolder}
fi

(time ./build/slabIPB \
    --gui false \
    --nParticles 10000000 \
    --minRange 0.01 \
    --maxStep 0.01 \
    --resultFolder ${experimentFolder} \
    --recordEventLog false \
    ) 2>&1 | tee "${experimentFolder}/myOutput.txt"