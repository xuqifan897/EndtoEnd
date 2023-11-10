#!/bin/bash

export DOSECALC_DATA="/data/qifan/projects/EndtoEnd/CCCS/data"
sparsity='1e-4'

PreProcess="/data/qifan/projects/EndtoEnd\
/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
DoseCalc="/data/qifan/projects/EndtoEnd\
/CCCS/build/dosecalc-beamlet/dosecalc-beamlet"

expFolder="/data/qifan/projects/EndtoEnd/results/CCCSclone"
dcmFolder="/data/qifan/projects/EndtoEnd/results/slabBench/slab_dicom"
voxelSize=0.25

preprocessLog="${expFolder}/preprocess.log"
dosecalcLog="${expFolder}/dosecalc.log"

# preprocess
cd ${expFolder}
( time ${PreProcess} \
    --dicom=${dcmFolder} \
    --beamlist="${expFolder}/beamlist.txt" \
    --structures="${expFolder}/structures.json" \
    --config="${expFolder}/config.json" \
    --bbox-roi="BODY" \
    --voxsize=${voxelSize} \
    --verbose ) > ${preprocessLog} 2>&1

( time ${DoseCalc} \
    --sparsity-threshold=${sparsity} \
    --ndevices=1 \
    --device=1 \
    --verbose ) > ${dosecalcLog} 2>&1