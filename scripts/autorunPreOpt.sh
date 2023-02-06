#!/bin/bash

expFolder="/data/datasets/UCLAPatients/experiment"
nPatients=8
cd ./BOO
for ((i=1; i<=${numPatients}; i++))
do
    patientName="patient${i}"
    logFile="${expFolder}/${patientName}/preOpt.log"
    (time matlab -nodesktop -nodisplay -r "idx=${i};QX_PreOpgimize_IMRT") 2>&1 | tee ${logFile}
done