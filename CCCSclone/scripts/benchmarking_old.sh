#!/bin/bash

./build/CCCSclone \
    --debugFolder /data/qifan/projects/EndtoEnd/results/CCCSComp/new \
    --resultFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/results \
    --dataFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/tmp \
    --deviceIdx 3 \
    --dicomVolumeDimension 103 103 103 \
    --voxelSize 0.25 0.25 0.25 \
    --doseBoundingBoxStartIndices 1 1 1 \
    --doseBoundingBoxDimensions 101 101 101 \
    --REVConvolutionArrayDimensions 800 800 800 \
    --debugLog false \
    --debugREVTerma false \
    --debugREVDose false \
    --debugBEVDose false \