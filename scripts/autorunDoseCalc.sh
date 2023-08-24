#!/bin/bash

# please specify the parameters in the block below
numPatients=8
expFolder='/data/datasets/UCLAPatients/experiment'
DataFolder='/data/datasets/UCLAPatients/anonymousDataNew'
projectFolder='/data/qifan/projects/EndtoEnd4'
voxelsize='0.25'
sparsity='1e-4'
CUDAdevice=3


# echo ${expFolder}
prepExec="${projectFolder}/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
doseCalcExec="${projectFolder}/CCCS/build/dosecalc-beamlet/dosecalc-beamlet"
export CUDA_VISIBLE_DEVICES=${CUDAdevice}
export DOSECALC_DATA="${projectFolder}/CCCS/data"

# for ((i=1; i<=${numPatients}; i++))
for i in 5
do
    patientName="patient${i}"
    dicomFolder="${DataFolder}/${patientName}/CTAlignResize"
    patExpFolder="${expFolder}/${patientName}"
    beamlist="${patExpFolder}/beamlist.txt"
    configfile="${patExpFolder}/config.json"
    structuresFile="${patExpFolder}/structures.json"
    BBox="Skin"

    # preprocessing
    cd ${patExpFolder}
    (time ${prepExec} \
    --dicom=${dicomFolder} \
    --beamlist=${beamlist} \
    --structures=${structuresFile} \
    --config=${configfile} \
    --bbox-roi=${BBox} \
    --voxsize=${voxelsize} \
    --verbose ) 2>&1 | tee 'dosecalc-preprocess.log'

    # dose calculation
    (time ${doseCalcExec} \
    --sparsity-threshold=${sparsity} \
    ) 2>&1 | tee 'dosecalc-beamlet.log'
done