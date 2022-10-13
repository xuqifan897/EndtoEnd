import os
import sys


def autorun_annealing():
    template = '#!/bin/bash\n\
export CUDA_VISIBLE_DEVICES=3\n\
./build/optimize_annealing_constrain \\\n\
    --phantom-dimension D1 D2 D3 \\\n\
    --voxel-size 2. \\\n\
    --phantom-isocenter ISO1 ISO2 ISO3 \\\n\
    --phantom-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/CT.dat \\\n\
    --PTV-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/PTV_weight.dat \\\n\
    --PTV-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/PTV_target.dat \\\n\
    --OAR-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/OAR_weight_correct.dat \\\n\
    --OAR-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/OAR_target.dat \\\n\
    --beam-angle-config-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_E2E/beam_angles_uniform.txt \\\n\
    --beam-energy 6. \\\n\
    --SAD 1000. \\\n\
    --fluence-map-dimension 128 128 \\\n\
    --fluence-map-convolution-radius 64 64 \\\n\
    --fluence-map-sampling-range 680 1320 \\\n\
    --fluence-map-sampling-points 640 \\\n\
    --fluence-map-pixel-size 0.7815 0.7815 \\\n\
    --output-folder /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient1_annealing_correct \\\n\
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \\\n\
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \\\n\
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \\\n\
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \\\n\
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \\\n\
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \\\n\
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \\\n\
    --iterations 1000 \\\n\
    --step-size 5e-3 \\\n\
    --ideal-dose 0.0 \\\n\
    --eta 1e3 \\\n\
    --step-size-angular-max 0.6 \\\n\
    --step-size-angular-min 0.06 \\\n\
    --temperature 1e5 \\\n\
    --temperature-max 4e5 \\\n\
    --temperature-min 5e4 \\\n\
    --zenith-range 0.524 2.618'

    for i in range(4, 7):
        input_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient{}_E2E'.format(i)
        input_shape_file = os.path.join(input_folder, 'shape.txt')
        input_iso_file = os.path.join(input_folder, 'isocenter.txt')

        with open(input_shape_file, 'r') as f:
            lineShape = f.readline()
        lineShape = lineShape.split(' ')

        with open(input_iso_file, 'r') as f:
            lineIso = f.readline()
        lineIso = lineIso.split(',')

        command = template.replace('D1', lineShape[0])
        command = command.replace('D2', lineShape[1])
        command = command.replace('D3', lineShape[2])
        command = command.replace('ISO1', lineIso[0])
        command = command.replace('ISO2', lineIso[1])
        command = command.replace('ISO3', lineIso[2])
        command = command.replace('patient1', 'patient{}'.format(i))
        os.system(command)


def autorun_annealing_constrain_init():
    template = 'export CUDA_VISIBLE_DEVICES=2\n\
./build/optimize_stationary_smoothness \\\n\
    --phantom-dimension D1 D2 D3 \\\n\
    --voxel-size 2. \\\n\
    --phantom-isocenter ISO1 ISO2 ISO3 \\\n\
    --phantom-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/CT.dat \\\n\
    --PTV-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/PTV_weight.dat \\\n\
    --PTV-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/PTV_target.dat \\\n\
    --OAR-weight-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/OAR_weight.dat \\\n\
    --OAR-target-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/OAR_target.dat \\\n\
    --beam-angle-config-path /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_E2E/beam_angles_annealing_constrain.txt \\\n\
    --beam-energy 6. \\\n\
    --SAD 1000. \\\n\
    --fluence-map-dimension 128 128 \\\n\
    --fluence-map-convolution-radius 64 64 \\\n\
    --fluence-map-sampling-range 680 1320 \\\n\
    --fluence-map-sampling-points 640 \\\n\
    --fluence-map-pixel-size 0.7815 0.7815 \\\n\
    --output-folder /home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient4_annealing_constrain_init \\\n\
    --spectrum-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/Spectrum.csv \\\n\
    --ATheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperATheta.csv \\\n\
    --atheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerATheta.csv \\\n\
    --BTheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/upperBTheta.csv \\\n\
    --btheta-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/lowerBTheta.csv \\\n\
    --pencil-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/FCBBkernel.csv \\\n\
    --depthDose-path /data/qifan/projects_qlyu/EndtoEnd3/kernels/depthDose.csv \\\n\
    --iterations 200 \\\n\
    --step-size 5e-3 \\\n\
    --ideal-dose 0.0 \\\n\
    --eta 1e3'

    for i in range(1, 7):
        input_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation/patient{}_E2E'.format(i)
        input_shape_file = os.path.join(input_folder, 'shape.txt')
        input_iso_file = os.path.join(input_folder, 'isocenter.txt')

        with open(input_shape_file, 'r') as f:
            lineShape = f.readline()
        lineShape = lineShape.split(' ')

        with open(input_iso_file, 'r') as f:
            lineIso = f.readline()
        lineIso = lineIso.split(',')

        command = template.replace('D1', lineShape[0])
        command = command.replace('D2', lineShape[1])
        command = command.replace('D3', lineShape[2])
        command = command.replace('ISO1', lineIso[0])
        command = command.replace('ISO2', lineIso[1])
        command = command.replace('ISO3', lineIso[2])
        command = command.replace('patient4', 'patient{}'.format(i))
        os.system(command)
        # print(command)
        # break


if __name__ == '__main__':
    # autorun_annealing()
    # autorun_annealing_constrain_init()
    autorun_annealing()