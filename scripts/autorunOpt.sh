#!/bin/bash
# expFolder="/data/datasets/UCLAPatients/experiment"
# numPatients=8
# trailNO=1
# cd ./BOO

# for ((i=5; i<=${numPatients}; i++))
# do
#     patExpFolder="${expFolder}/patient${i}/optimize"
#     logFile="${patExpFolder}/optim${trailNO}.log"
#     (time matlab -nodesktop -nodisplay -r "patientIdx=${i};trailNO=${trailNO};QX_4piIMRT_cpu;exit") 2>&1 | tee ${logFile}
# done


# for patients 3 and 5, there are more beams sellected.
# increased beamWeight to select 20 beams
expFolder="/data/datasets/UCLAPatients/experiment"
numPatients=8
trailNO=2
cd ./BOO

# for ((i=5; i<=${numPatients}; i++))
for i in 3 5
do
    patExpFolder="${expFolder}/patient${i}/optimize"
    logFile="${patExpFolder}/optim${trailNO}.log"
    (time matlab -nodesktop -nodisplay -r "patientIdx=${i};trailNO=${trailNO};QX_4piIMRT_cpu;exit") 2>&1 | tee ${logFile}
done