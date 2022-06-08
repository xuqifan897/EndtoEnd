#!/bin/bash

gdb --args /data/qifan/projects_qlyu/EndtoEnd3/build/EndtoEnd3 \
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
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 512 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --fluence-map-output-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/fluence_map.dat \
    --zenith-range 30 150 \
    --dose-path /data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/dose.dat \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv


# tbreak /data/qifan/projects_qlyu/EndtoEnd3/args/argparse.cpp:65
# tbreak /data/qifan/projects_qlyu/EndtoEnd3/args/kernelInit.cpp:81