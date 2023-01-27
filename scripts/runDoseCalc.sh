#!/bin/bash

# this file runs dose calculation for all the patients
# please change the settings below
numPatients=8
expFolder="/data/qifan/dataset_qlyu/UCLAPatients/experiment"
DataFolder="/data/qifan/dataset_qlyu/UCLAPatients/anonymousDataNew"
projectFolder="/data/qifan/projects_qlyu/EndtoEnd4"


prepExec="${projectFolder}/CCCS/build/dosecalc-preprocess/dosecalc-preprocess"
doseCalcExec="${projectFolder}/CCCS/build/dosecalc-beam/dosecalc-beam"

for ((i=1; i<=${numPatients}; i++))
do
    patientName="patient${i}"
    dicomFolder="${DataFolder}/${patientName}/CTAlignResize"
    patExpFolder="${expFolder}/${patientName}"
    echo ${dicomFolder}
    echo ${patExpFolder}
done