import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
# import mat73


def examineMask():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation_speedup'
    num_patients = 6
    OARweight = 10
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(parent_folder, '{}_BOO'.format(patientName))
        if i == 0:
            patFolder = os.path.join(parent_folder, 'patient11_BOO')
        debugFolder = os.path.join(patFolder, 'debug')
        if not os.path.isdir(debugFolder):
            os.mkdir(debugFolder)
        optFolder = os.path.join(patFolder, 'optimize_OAR{}'.format(OARweight))
        maskPath = os.path.join(patFolder, 'masks.mat')
        masks = scipy.io.loadmat(maskPath)['masks']
        for j in range(len(masks)):
            mask = masks[j]
            name = mask[0][0][0][0][0]
            array = mask[0][0][0][1] # dtype: uint8
            if name == 'BODY':
                BODYmask_path = os.path.join(debugFolder, 'BODYmask')
                if not os.path.isdir(BODYmask_path):
                    os.mkdir(BODYmask_path)
                array = array * 255 # convert 1 to 255
                for j in range(array.shape[2]):
                    filePath = os.path.join(BODYmask_path, '{:03d}.png'.format(j))
                    plt.imsave(filePath, array[:, :, j])


def examineDose():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation_speedup'
    num_patients = 6
    OARweight = 10
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(parent_folder, '{}_BOO'.format(patientName))
        if i == 0:
            patFolder = os.path.join(parent_folder, 'patient11_BOO')
        debugFolder = os.path.join(patFolder, 'debug')
        if not os.path.isdir(debugFolder):
            os.mkdir(debugFolder)
        optFolder = os.path.join(patFolder, 'optimize_OAR{}'.format(OARweight))
        resultFile = os.path.join(optFolder, 
            'result {} Info0 params1 beam20.mat'.format(patientName))
        dose = scipy.io.loadmat(resultFile)['result']['dose'][0][0]
        dosePath = os.path.join(debugFolder, 'dose')
        if not os.path.isdir(dosePath):
            os.mkdir(dosePath)
        doseRoof = np.max(dose)
        for j in range(dose.shape[2]):
            filePath = os.path.join(dosePath, '{:03d}.png'.format(j))
            plt.imsave(filePath, dose[:, :, j], vmin=0, vmax=doseRoof)
        print(patientName)


def draw_DVH():
    parent_folder = '/home/qlyu/ShengNAS2/SharedProjectData/QX_beam_orientation_speedup'
    num_patients = 6
    OARweight = 8
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
        'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
    
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(parent_folder, '{}_BOO'.format(patientName))
        if i == 0:
            patFolder = os.path.join(parent_folder, 'patient11_BOO')
        optFolder = os.path.join(patFolder, 'optimize_OAR{}'.format(OARweight))
        debugFolder = os.path.join(patFolder, 'debug')
        assert os.path.isdir(debugFolder)
        resultFile = os.path.join(optFolder, 
            'result {} Info0 params1 beam20.mat'.format(patientName))
        dose = scipy.io.loadmat(resultFile)['result']['dose'][0][0]

        maskPath = os.path.join(patFolder, 'masks.mat')
        masks = scipy.io.loadmat(maskPath)['masks']
        namelist = []
        for j in range(len(masks)):
            mask = masks[j]
            name = mask[0][0][0][0][0]
            if name == 'BODY':
                continue
            namelist.append(name)
            array = mask[0][0][0][1] # dtype: uint8
            array = array.astype(bool)
            selected = dose[array]
            selected.sort()

            # insert head and tail
            selected = np.insert(selected, 0, 0)
            if selected[-1] < 20:
                selected = np.insert(selected, -1, 20)

            n_voxels = len(selected)
            x_axis = (n_voxels - 1 - np.arange(n_voxels)) / (n_voxels - 1)
            color = colors[j]
            plt.plot(selected, x_axis, color)
        plt.legend(namelist)
        figurePath = os.path.join(debugFolder, 'DVH_OAR{}.png'.format(OARweight))
        plt.savefig(figurePath)
        plt.clf()
        print(patientName)



if __name__ == '__main__':
    # examineMask()
    # examineDose()
    draw_DVH()