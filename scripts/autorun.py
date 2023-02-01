import os
import json

"""
This function defines the script for dose calculation using Ryan Neph's code
"""

# def autorun():
#     template1 = \
# """#!/bin/bash

# cd {}
# export CUDA_VISIBLE_DEVICES=1
# export DOSECALC_DATA=\"/data/qifan/projects_qlyu/dose-calculation/v0.8/data\"
# BBOX=\"BODY\"
# dicomdata=\"{}\"
# configfile=\"{}\"
# structures=\"{}\"
# beamlist=\"{}\"

# voxelsize='0.25' # unit: cm
# sparsity='1e-4'
# """

#     template2 = \
# """
# # preprocess
# time /data/qifan/projects_qlyu/dose-calculation/build/dosecalc-preprocess/dosecalc-preprocess \\
#     --dicom=${dicomdata} \\
#     --beamlist=${beamlist} \\
#     --structures=${structures} \\
#     --config=${configfile} \\
#     --bbox-roi=${BBOX} \\
#     --voxesize=${voxelsize} \\
#     --verbose \\
#     2>&1 | tee "dosecalc-preprocess.log"

# echo -e \"\\n\\n=================================================================================\n\n\"

# # beamlet dose calculation
# time /data/qifan/projects_qlyu/dose-calculation/build/dosecalc-beamlet/dosecalc-beamlet \\
#     --sparsity-threshold=${sparsity} \\
#     2>&1 | tee "dosecalc-beamlet.log"

# echo -e \"\\n\\n=================================================================================\n\n\"
# """
#     parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation_speedup'
#     num_patients = 6
#     for i in range(num_patients):
#         patient_idx = i + 1
#         target_folder = os.path.join(parent_folder, 'patient{}_BOO'.format(patient_idx))
#         if patient_idx == 1:
#             target_folder = os.path.join(parent_folder, 'patient11_BOO')
#         dicom_folder = os.path.join(parent_folder, 'patient{}_dicom'.format(patient_idx))
#         config_file = os.path.join(target_folder, 'config.json')
#         structures_file = os.path.join(target_folder, 'structures.json')
#         beamlist_file = os.path.join(target_folder, 'beamlist.txt')

#         script = template1.format(target_folder, dicom_folder, 
#             config_file, structures_file, beamlist_file) + template2
        
#         os.system(script)
#         break


# def timing_experiment():
#     command = 'time ls 2>&1 | tee log.txt'
#     os.system(command)

def autorun():
    """
    This function runs dose calculation for all the patients
    please change the variables below
    """
    numPatients = 8
    expFolder = "/data/qifan/dataset_qlyu/UCLAPatients/experiment"
    DataFolder = "/data/qifan/dataset_qlyu/UCLAPatients/anonymousDataNew"
    projectFolder = "/data/qifan/projects_qlyu/EndtoEnd4"
    voxelsize = '0.25' # unit: cm, 0.25 is the recommended value
    sparsity = '1e-4'
    # debug = True
    debug = False

    # All variables below are derived
    DOSECALC_DATA = os.path.join(projectFolder, 'CCCS', 'data')
    command = 'export DOSECALC_DATA={}'.format(DOSECALC_DATA)
    os.system(command)
    prepExec = os.path.join(projectFolder, 'CCCS', 'build', 
        'dosecalc-preprocess/dosecalc-preprocess')
    doseCalcExec = os.path.join(projectFolder, 'CCCS', 'build', 'dosecalc-beamlet', 'dosecalc-beamlet')

    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        # the user may also revise the command below, 'CTAlignResize' specifically
        dicomFolder = os.path.join(DataFolder, patientName, 'CTAlignResize')
        patExpFolder = os.path.join(expFolder, patientName)
        beamlist = os.path.join(patExpFolder, 'beamlist.txt')
        configfile = os.path.join(patExpFolder, 'config.json')
        structuresFile = os.path.join(patExpFolder, 'structures.json')
        BBox='Skin'
        CUDAdevice = 1

        with open(structuresFile, 'r') as f:
            data = json.load(f)
        PTVname = data['ptv']

        # preprocessing
        command = \
        f"""export DOSECALC_DATA="{DOSECALC_DATA}"
    (time {prepExec} \\
    --dicom={dicomFolder} \\
    --beamlist={beamlist} \\
    --structures={structuresFile} \\
    --config={configfile} \\
    --bbox-roi={BBox} \\
    --voxsize={voxelsize} \\
    --verbose \\
    ) 2>&1 | tee "{os.path.join(patExpFolder, 'dosecalc-preprocess.log')}"
"""
        if debug:
            print(command)
        else:
            os.system(command)

        command = 'echo "\n\n{}\n\n"'.format(40*'=')
        if debug:
            print(command)
        else:
            os.system(command)

        # dose calculation
        command = f'cd {patExpFolder}\n' \
            f'export CUDA_VISIBLE_DEVICES={CUDAdevice}\n' \
            f'export DOSECALC_DATA="{DOSECALC_DATA}"\n' \
            f'( time {doseCalcExec} --sparsity-threshold={sparsity} ) ' \
            f'2>&1 | tee {os.path.join(patExpFolder, "dosecalc-beamlet.log")}'
        if debug:
            print(command)
        else:
            os.system(command)


if __name__ == '__main__':
    autorun()