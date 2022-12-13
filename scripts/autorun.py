import os

"""
This function defines the script for dose calculation using Ryan Neph's code
"""

def autorun():
    template1 = \
"""#!/bin/bash

cd {}
export CUDA_VISIBLE_DEVICES=1
export DOSECALC_DATA=\"/data/qifan/projects_qlyu/dose-calculation/v0.8/data\"
BBOX=\"BODY\"
dicomdata=\"{}\"
configfile=\"{}\"
structures=\"{}\"
beamlist=\"{}\"

voxelsize='0.25' # unit: cm
sparsity='1e-4'
"""

    template2 = \
"""
# preprocess
time /data/qifan/projects_qlyu/dose-calculation/build/dosecalc-preprocess/dosecalc-preprocess \\
    --dicom=${dicomdata} \\
    --beamlist=${beamlist} \\
    --structures=${structures} \\
    --config=${configfile} \\
    --bbox-roi=${BBOX} \\
    --voxesize=${voxelsize} \\
    --verbose \\
    2>&1 | tee "dosecalc-preprocess.log"

echo -e \"\\n\\n=================================================================================\n\n\"

# beamlet dose calculation
time /data/qifan/projects_qlyu/dose-calculation/build/dosecalc-beamlet/dosecalc-beamlet \\
    --sparsity-threshold=${sparsity} \\
    2>&1 | tee "dosecalc-beamlet.log"

echo -e \"\\n\\n=================================================================================\n\n\"
"""
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation_speedup'
    num_patients = 6
    for i in range(num_patients):
        patient_idx = i + 1
        target_folder = os.path.join(parent_folder, 'patient{}_BOO'.format(patient_idx))
        if patient_idx == 1:
            target_folder = os.path.join(parent_folder, 'patient11_BOO')
        dicom_folder = os.path.join(parent_folder, 'patient{}_dicom'.format(patient_idx))
        config_file = os.path.join(target_folder, 'config.json')
        structures_file = os.path.join(target_folder, 'structures.json')
        beamlist_file = os.path.join(target_folder, 'beamlist.txt')

        script = template1.format(target_folder, dicom_folder, 
            config_file, structures_file, beamlist_file) + template2
        
        os.system(script)
        break


def timing_experiment():
    command = 'time ls 2>&1 | tee log.txt'
    os.system(command)

if __name__ == '__main__':
    autorun()
    # timing_experiment()