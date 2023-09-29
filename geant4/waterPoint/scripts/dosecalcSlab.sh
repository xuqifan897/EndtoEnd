#!/bin/bash

export DOSECALC_DATA="/data/qifan/projects/EndtoEnd/CCCS/data"
sparsity='1e-4'

PreProcess="/data/qifan/projects/EndtoEnd\
/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
DoseCalc="/data/qifan/projects/EndtoEnd\
/CCCS/build/dosecalc-beamlet/dosecalc-beamlet"

dcmFolder="/data/qifan/projects/EndtoEnd/results/slabBench/slab_dicom"
expTemplate="/data/qifan/projects/EndtoEnd/results/slabBench/slab_dosecalc"
voxelSize=0.1


# customize the data, fluence dimension and beamlet size
fluenceDim=9
beamletSize=1.0

expFolder="/data/qifan/projects/EndtoEnd\
/results/slabBench/slab_dosecalc_${fluenceDim}_${beamletSize}"
preprocessLog="${expFolder}/preprocess.log"
dosecalcLog="${expFolder}/dosecalc.log"
if [ ! -d ${expFolder} ];
then
    mkdir ${expFolder}
fi
cp "${expTemplate}/beamlist.txt" ${expFolder}
sed "s/bs/${beamletSize}/g; s/fd/${fluenceDim}/g" "${expTemplate}/config.json" > "${expFolder}/config.json"
cp "${expTemplate}/structures.json" ${expFolder}

# preprocess
cd ${expFolder}
( time ${PreProcess} \
    --dicom=${dcmFolder} \
    --beamlist="${expFolder}/beamlist.txt" \
    --structures="${expFolder}/structures.json" \
    --config="${expFolder}/config.json" \
    --bbox-roi="water" \
    --voxsize=${voxelSize} \
    --verbose ) > ${preprocessLog}

( time ${DoseCalc} \
    --sparsity-threshold=${sparsity} \
    --ndevices=1 \
    --device=1 \
    --verbose ) > ${dosecalcLog}