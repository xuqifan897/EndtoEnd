#!/bin/bash

export CUDA_VISIBLE_DEVICES=2

# cuda-gdb --args ./build/dose_calculation \
./build/optimize_stationary_smoothness \
    --phantom-dimension 170 170 113 \
    --voxel-size 2. \
    --phantom-isocenter 159.22469209972965 149.11895464103335 123.0756984079303 \
    --phantom-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/CT.dat \
    --PTV-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/PTV_weight.dat \
    --PTV-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/PTV_target.dat \
    --OAR-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/OAR_weight.dat \
    --OAR-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/OAR_target.dat \
    --beam-angle-config-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/beam_angles_E2E.txt \
    --beam-energy 6. \
    --SAD 1000. \
    --fluence-map-dimension 128 128 \
    --fluence-map-convolution-radius 64 64 \
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 640 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --output-folder /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_BOO_init \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \
    --iterations 200 \
    --step-size 5e-3 \
    --ideal-dose 0.0 \
    --eta 1e3
    # --fluence-map-init /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/fluence_maps \


./build/optimize_stationary_smoothness \
    --phantom-dimension 200 200 152 \
    --voxel-size 2. \
    --phantom-isocenter 201.29726862901086 190.2460885706709 142.41156191991513 \
    --phantom-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_E2E/CT.dat \
    --PTV-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_E2E/PTV_weight.dat \
    --PTV-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_E2E/PTV_target.dat \
    --OAR-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_E2E/OAR_weight.dat \
    --OAR-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_E2E/OAR_target.dat \
    --beam-angle-config-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_E2E/beam_angles_E2E.txt \
    --beam-energy 6. \
    --SAD 1000. \
    --fluence-map-dimension 128 128 \
    --fluence-map-convolution-radius 64 64 \
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 640 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --output-folder /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient5_BOO_init \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \
    --iterations 200 \
    --step-size 5e-3 \
    --ideal-dose 0.0 \
    --eta 1e3


./build/optimize_stationary_smoothness \
    --phantom-dimension 260 260 128 \
    --voxel-size 2. \
    --phantom-isocenter 258.6639053254438 250.42011834319527 121.75029585798816 \
    --phantom-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_E2E/CT.dat \
    --PTV-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_E2E/PTV_weight.dat \
    --PTV-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_E2E/PTV_target.dat \
    --OAR-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_E2E/OAR_weight.dat \
    --OAR-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_E2E/OAR_target.dat \
    --beam-angle-config-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_E2E/beam_angles_E2E.txt \
    --beam-energy 6. \
    --SAD 1000. \
    --fluence-map-dimension 128 128 \
    --fluence-map-convolution-radius 64 64 \
    --fluence-map-sampling-range 680 1320 \
    --fluence-map-sampling-points 640 \
    --fluence-map-pixel-size 0.7815 0.7815 \
    --output-folder /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient6_BOO_init \
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \
    --iterations 200 \
    --step-size 5e-3 \
    --ideal-dose 0.0 \
    --eta 1e3