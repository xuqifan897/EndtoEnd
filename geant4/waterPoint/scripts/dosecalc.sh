#!/bin/bash

PreProcess="/data/qifan/projects/EndtoEnd\
/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
DoseCalc="/data/qifan/projects/EndtoEnd\
/CCCS/build/dosecalc-beamlet/dosecalc-beamlet"

expFolder="/data/qifan/projects/EndtoEnd\
/results/slabBench/patien1_dosecalc"
dcmFolder="/data/qifan/projects/EndtoEnd/results/slabBench/patient1_dicom"
voxelSize=0.1

echo ${PreProcess}

# preprocess
cd ${expFolder}
gdb --args ${PreProcess} \
    --dicom=${dcmFolder} \
    --beamlist="${expFolder}/beamlist.txt" \
    --structures="${expFolder}/structures.json" \
    --config="${expFolder}/config.json" \
    --bbox-roi="BODY" \
    --voxsize=${voxelSize} \
    --verbose