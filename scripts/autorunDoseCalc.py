import os
import json

"""
This function defines the script for dose calculation using Ryan Neph's code
"""

def autorun():
    """
    This function runs dose calculation for all the patients
    please change the variables below
    """
    ## runs on shenggpu4
    numPatients = 8
    expFolder = "/data/datasets/UCLAPatients/experiment"
    DataFolder = "/data/datasets/UCLAPatients/anonymousDataNew"
    projectFolder = "/data/qifan/projects/EndtoEnd4"
    voxelsize = '0.25' # unit: cm, 0.25 is the recommended value
    sparsity = '1e-4'
    # debug = True
    debug = False
    

    ## runs on shenggpu2
    # numPatients = 8
    # expFolder = "/data/qifan/dataset_qlyu/UCLAPatients/experiment"
    # DataFolder = "/data/qifan/dataset_qlyu/UCLAPatients/anonymousDataNew"
    # projectFolder = "/data/qifan/projects_qlyu/EndtoEnd4"
    # voxelsize = '0.25' # unit: cm, 0.25 is the recommended value
    # sparsity = '1e-4'
    # # debug = True
    # debug = False

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
        # break

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