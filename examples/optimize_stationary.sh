#!/bin/bash

export CUDA_VISIBLE_DEVICES=3

# cuda-gdb --args ./build/dose_calculation \
./build/optimize_stationary \
    --phantom-dimension 200 200 197 \
    --voxel-size 2. \
    --phantom-isocenter 205.95534 211.23352 162.16011 \
    --phantom-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/CT.dat \
    --PTV-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/PTV_weight.dat \
    --PTV-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/PTV_target.dat \
    --OAR-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/OAR_weight.dat \
    --OAR-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/OAR_target.dat \
    --beam-angle-config-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/beam_angles_E2E.txt \
    --beam-energy 6. \
    --SAD 1000. \
    --fluence-map-dimension 128 128 \
    --fluence-map-convolution-radius 64 64 \
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 640 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --output-folder /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_optimize_stationary \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \
    --iterations 60 \
    --step-size 5e-3 \
    # --fluence-map-init /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/fluence_maps \