#!/bin/bash

# # For IPB kernel
# (time ./build/waterKernel \
#     --gui false \
#     --sizeZ 20.0 \
#     --posZ -20.0 \
#     --recordEventLog false \
#     --nParticles 10000000 \
#     --resultFolder /data/qifan/projects/EndtoEnd4/results/IPB6MeV1e7 \
#     ) 2>&1 | tee /data/qifan/projects/EndtoEnd4/results/IPB6MeV1e7/myOutput.txt

# rm *.rndm

# # For point kernel
# resultFolder="/data/qifan/projects/EndtoEnd4/results/point6MeV1e7Step"
# if [ ! -d ${resultFolder} ]
# then
#     mkdir ${resultFolder}
# fi
# (time ./build/waterKernel \
#     --gui false \
#     --sizeZ 60.0 \
#     --posZ -60.0 \
#     --kernelType point \
#     --kernelSizeZ 20 \
#     --kernelPosZ -15 \
#     --maxStep 0.01 \
#     --recordEventLog false \
#     --nParticles 10000000 \
#     --resultFolder "${resultFolder}" \
#     ) 2>&1 | tee "${resultFolder}/myOutput.txt"

# rm *.rndm

# for debugging
resultFolder="/data/qifan/projects/EndtoEnd4/results/point6MeVdebug"
if [ ! -d "${resultFolder}" ]
then
    mkdir "${resultFolder}"
fi

(time ./build/waterKernel \
    --gui false \
    --sizeZ 60.0 \
    --posZ -60.0 \
    --kernelType point \
    --kernelSizeZ 20 \
    --maxStep 0.01 \
    --kernelPosZ -15 \
    --nParticles 32 \
    --resultFolder "${resultFolder}" \
    # --recordEventLog true \
    # --debugEventLogTrajectory true \
    ) 2>&1 | tee "${resultFolder}/myOutput.txt"

rm *.rndm