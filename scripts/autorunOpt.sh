#!/bin/bash
expFolder="/data/datasets/UCLAPatients/experiment"
numPatients=8
trailNO=1
cd ./BOO

for ((i=1; i<=${numPatients}; i++))
do
    patExpFolder="${expFolder}/patient${i}/optimize"
    logFile="${patExpFolder}/optim${trailNO}.log"
    (time matlab -nodesktop -nodisplay -r "patientIdx=${i};trailNO=${trailNO};QX_4piIMRT_cpu;exit") 2>&1 | tee ${logFile}
done