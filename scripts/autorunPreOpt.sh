#!/bin/bash

expFolder="/data/datasets/UCLAPatients/experiment"
numPatients=8
cd ./BOO
for ((i=7; i<=${numPatients}; i++))
do
    patientName="patient${i}"
    logFile="${expFolder}/${patientName}/preOpt.log"
    (time matlab -nodesktop -nodisplay -r "idx=${i};QX_PreOptimize_IMRT") 2>&1 | tee ${logFile}
done