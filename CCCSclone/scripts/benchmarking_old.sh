#!/bin/bash

./build/CCCSclone \
    --debugLog true \
    --resultFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/results \
    --dataFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/tmp \
    --dicomVolumeDimension 103 103 103 \
    --voxelSize 0.25 0.25 0.25 \
    --doseBoundingBoxStartIndices 1 1 1 \
    --doseBoundingBoxDimensions 101 101 101 \
    --REVConvolutionArrayDimensions 800 800 800 \