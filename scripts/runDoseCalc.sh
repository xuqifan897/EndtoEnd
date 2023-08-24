#!/bin/bash

# this file runs dose calculation for all the patients
# please change the settings below
numPatients=8
expFolder="/data/qifan/dataset_qlyu/UCLAPatients/experiment"
DataFolder="/data/qifan/dataset_qlyu/UCLAPatients/anonymousDataNew"
projectFolder="/data/qifan/projects_qlyu/EndtoEnd4"
voxelsize='0.25' # unit: cm, 0.25 is the recommended value
sparsity='1e-4'

export DOSECALC_DATA="${projectFolder}/CCCS/data"

prepExec="${projectFolder}/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
doseCalcExec="${projectFolder}/CCCS/build/dosecalc-beam/dosecalc-beam"

# for ((i=1; i<=${numPatients}; i++))
for i in 1
do
    patientName="patient${i}"
    dicomFolder="${DataFolder}/${patientName}/CTAlignResize"
    patExpFolder="${expFolder}/${patientName}"
    # echo ${dicomFolder}
    # echo ${patExpFolder}

    beamlist="${patExpFolder}/beamlist.txt"
    structures="${patExpFolder}/structures.json"
    configfile="${patExpFolder}/config.json"
    BBox="Skin"

    cd ${patExpFolder}
    (time ${prepExec} \
        --dicom=${dicomFolder} \
        --beamlist=${beamlist} \
        --target-exact="PTV 50 uncropped" \
        --config=${configfile} \
        --bbox-roi=${BBox} \
        --voxsize=${voxelsize} \
        --verbose \
    ) 2>&1 | tee "dosecalc-preprocess.log"
done