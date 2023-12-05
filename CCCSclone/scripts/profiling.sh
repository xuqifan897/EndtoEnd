#!/bin/bash
# sudo sh -c 'echo 1 >/proc/sys/kernel/perf_event_paranoid'

debugFolder="/data/qifan/projects/EndtoEnd/results/CCCSComp/new"
resultFolder="/data/qifan/projects/EndtoEnd/results/CCCSclone/results"
dataFolder="/data/qifan/projects/EndtoEnd/results/CCCSclone/tmp"
profileResult="${resultFolder}/prof"

if [ ! -d ${resultFolder} ]; then
    mkdir ${resultFolder}
fi

# nsys profile --stats=true -o ${profileResult} ./build/CCCSclone \
#     --debugFolder /data/qifan/projects/EndtoEnd/results/CCCSComp/new \
#     --resultFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/results \
#     --dataFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/tmp \
#     --deviceIdx 3 \
#     --dicomVolumeDimension 103 103 103 \
#     --voxelSize 0.25 0.25 0.25 \
#     --doseBoundingBoxStartIndices 1 1 1 \
#     --doseBoundingBoxDimensions 101 101 101 \
#     --REVConvolutionArrayDimensions 800 800 800 \

computeResult="${resultFolder}/compute.ncu-rep"
/usr/local/cuda/bin/ncu -f \
    --target-processes all \
    --section SchedulerStats \
    --section WarpStateStats \
    --section SourceCounters \
    --section Occupancy \
    -o ${computeResult} ./build/CCCSclone \
    --debugFolder /data/qifan/projects/EndtoEnd/results/CCCSComp/new \
    --resultFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/results \
    --dataFolder /data/qifan/projects/EndtoEnd/results/CCCSclone/tmp \
    --deviceIdx 3 \
    --dicomVolumeDimension 103 103 103 \
    --voxelSize 0.25 0.25 0.25 \
    --doseBoundingBoxStartIndices 1 1 1 \
    --doseBoundingBoxDimensions 101 101 101 \
    --REVConvolutionArrayDimensions 800 800 800 \