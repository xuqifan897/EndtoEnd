#!/bin/bas

source_folder="/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation"
num_patients=6;
for((i=1; i<=${num_patients}; i++))
do
    source_file="${source_folder}/patient${i}/optimize/*.csv"
    target_folder="${source_folder}/patient${i}_E2E/beam_angles_VarianIEC.csv"
    cp ${source_file} ${target_folder}
done