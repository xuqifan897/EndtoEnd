#!/bin/bash

export CUDA_VISIBLE_DEVICES=3

./build/EndtoEnd3 \
    --phantom-dimension 200 200 197 \
    --voxel-size 2. \
    --phantom-isocenter 205.95534 211.23352 162.16011 \
    --phantom-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1/CT.dat \
    --PTV-weight-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1/PTVweight.dat \
    --PTV-target-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1/PTVtarget.dat \
    --OAR-weight-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1/OARweight.dat \
    --OAR-target-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1/OARtarget.dat \
    --beam-energy 6. \
    --SAD 1000. \
    --number-of-beams 30 \
    --fluence-map-dimension 128 128 \
    --fluence-map-convolution-radius 64 64 \
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 640 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --output-folder /data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out \
    --zenith-range 30 150 \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \
    --beam-angle-config-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1/beamAngles.txt \
    --iterations 100 \
    --step-size 1e-2