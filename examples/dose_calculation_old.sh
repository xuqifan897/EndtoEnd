#!/bin/bash

export CUDA_VISIBLE_DEVICES=3

./build/dose_calculation \
    --phantom-dimension 200 200 200 \
    --voxel-size 2. \
    --phantom-isocenter 200. 200. 200. \
    --phantom-path /data/qifan/projects_qlyu/EndtoEnd3/data/water/CT.dat \
    --PTV-weight-path /data/qifan/projects_qlyu/EndtoEnd3/data/water/PTVweight.dat \
    --PTV-target-path /data/qifan/projects_qlyu/EndtoEnd3/data/water/PTVtarget.dat \
    --OAR-weight-path /data/qifan/projects_qlyu/EndtoEnd3/data/water/OARweight.dat \
    --OAR-target-path /data/qifan/projects_qlyu/EndtoEnd3/data/water/OARtarget.dat \
    --beam-energy 6. \
    --SAD 1000. \
    --number-of-beams 30 \
    --fluence-map-dimension 128 128 \
    --fluence-map-convolution-radius 64 64 \
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 640 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --output-folder /data/qifan/projects_qlyu/EndtoEnd3/data/water_out \
    --zenith-range 30 150 \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \
    --beam-angle-config-path /data/qifan/projects_qlyu/EndtoEnd3/data/water/beamAngles.txt \
    --iterations 10000 \
    --step-size 1e-2 \
    --azimuth 5.026548245743669