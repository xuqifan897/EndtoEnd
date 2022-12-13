#!/bin/bash

export DOSECALC_DATA="/home/qlyu/ShengNAS2/SharedProjectData/code_RoboticArm/dose_calculation/v0.8/data"

BBox="Body"
dicomdata="/data/qifan/projects_qlyu/dose-calculation/HN02/dosecalc/Data"
configfile="/data/qifan/projects_qlyu/dose-calculation/HN02/some_folder_beamlet/config.json"
#beamlist="./beam_list_SID50.txt"
structures="/data/qifan/projects_qlyu/dose-calculation/HN02/some_folder_beamlet/structures.json"

voxelsize='0.25' # unit: cm
sparsity='1e-4'

for i in 0 1 2 3 4
do
  beamlist="/data/qifan/projects_qlyu/dose-calculation/HN02/some_folder_beamlet/beam_files/beam_list_"${i}".txt"
#  echo $beamlist
  /data/qifan/projects_qlyu/dose-calculation/build/dosecalc-preprocess/dosecalc-preprocess \
    --dicom=$dicomdata \
    --beamlist=$beamlist \
    --structures=$structures \
    --config=$configfile \
    --bbox-roi=$BBox \
    --voxsize=$voxelsize \
    --verbose \
    2>&1 | tee "dosecalc-preprocess.log"

  echo -e "\n\n=================================================================================\n\n"

  /data/qifan/projects_qlyu/dose-calculation/build/dosecalc-beamlet/dosecalc-beamlet \
    --sparsity-threshold=$sparsity \
    2>&1 | tee "dosecalc-beamlet.log"

  echo -e "\n\n=================================================================================\n\n"

  mv ./Dose_Coefficients.h5 "./Dose_Coefficients_"$i".h5"
#  echo "./Dose_Coefficients_"$i".h5"
done