#!/bin/bash

benchFolder="/data/qifan/projects/EndtoEnd/results/CCCSBench"
if [ ! -d ${benchFolder} ]; then
    mkdir ${benchFolder}
fi

logFile="${benchFolder}/benchlog.txt"

nohup ./build/benchmark \
    --debugFolder ${benchFolder} \
    --resultFolder ${benchFolder} \
    --dataFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/tmp \
    --deviceIdx 3 \
    --dicomVolumeDimension 103 103 103 \
    --voxelSize 0.25 0.25 0.25 \
    --doseBoundingBoxStartIndices 1 1 1 \
    --doseBoundingBoxDimensions 101 101 101 \
    --REVConvolutionArrayDimensions 800 800 800 \
    --beamCount 100 \
2>&1 > ${logFile} & 