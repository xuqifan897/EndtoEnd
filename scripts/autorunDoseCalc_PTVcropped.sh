#!/bin/bash

# please specify the parameters in the block below
numPatients=8
expFolder='/data/datasets/UCLAPatients/experiment'
DataFolder='/data/datasets/UCLAPatients/anonymousDataNew'
projectFolder='/data/qifan/projects/EndtoEnd4'
voxelsize='0.25'
sparsity='1e-4'
CUDAdevice=0


# echo ${expFolder}
prepExec="${projectFolder}/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
doseCalcExec="${projectFolder}/CCCS/build/dosecalc-beamlet/dosecalc-beamlet"
export CUDA_VISIBLE_DEVICES=${CUDAdevice}
export DOSECALC_DATA="${projectFolder}/CCCS/data"

# for ((i=1; i<=${numPatients}; i++))
for ((i=2; i<=${numPatients}; i++))
do
    patientName="patient${i}"
    dicomFolder="${DataFolder}/${patientName}/CT_PTVcrop"
    patExpFolder="${expFolder}/${patientName}"

    optFolderNew="${patExpFolder}/optimizePTVcropped"
    beamlist="${optFolderNew}/beamlist.txt"
    structuresFile="${optFolderNew}/structures.json"
    configfile="${patExpFolder}/config.json"
    BBox="Skin"

    # preprocessing
    cd ${optFolderNew}
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