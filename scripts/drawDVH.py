import os
import numpy as np
import matplotlib.pyplot as plt
import pydicom
from scipy.io import loadmat
from rt_utils import RTStructBuilder
from PIL import Image
import json
import glob
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


def doseResize():
    """
    As the CCCS dose calculation method reshapes the dose,
    we have to reshape it back to the original shape
    """
    globalFolder = '/data/datasets/UCLAPatients'
    dataFolder = os.path.join(globalFolder, 'anonymousDataNew')
    expFolder = os.path.join(globalFolder, 'experiment')
    numPatients = 8
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimize')
        BOOFile, polishFile = findResultPath(optFolder)

        # get original dose
        polishData = loadmat(polishFile)['polishResult'][0, 0]
        polishDose = polishData['dose']
        print(polishDose.shape)

        # get original dose shape
        dataPatFolder = os.path.join(dataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        MRdose = pydicom.dcmread(MRdoseFile).pixel_array
        MRdoseShape = MRdose.shape
        targetShape = (MRdoseShape[1], MRdoseShape[2], MRdoseShape[0])
        print(targetShape)

        resizedDose = ReSize(polishDose, targetShape)
        # break
        
        # extract trailNO
        trailNO = 1
        BOOFile_ = BOOFile.split('/')[-1]
        while True:
            if str(trailNO) in BOOFile_:
                break
            else:
                trailNO += 1
        
        # save
        outFile = os.path.join(optFolder, 'dose{}.npy'.format(trailNO))
        np.save(outFile, resizedDose)


def ReSize(input, newShape):
    """
    This function resizes a 3D array using trilinear 
    interpoliation. Inputs:
    input: the input 3D array
    newShape: target shape
    output: resized input of the target shape
    """

    newIndices = np.indices(newShape)
    newIndices = newIndices.astype(np.float32)
    oldShape = input.shape

    oldIndices = np.zeros_like(newIndices)
    oldIndices[0, :, :, :] = (newIndices[0, :, :, :] + 0.5) \
        / newShape[0] * oldShape[0] - 0.5
    oldIndices[1, :, :, :] = (newIndices[1, :, :, :] + 0.5) \
        / newShape[1] * oldShape[1] - 0.5
    oldIndices[2, :, :, :] = (newIndices[2, :, :, :] + 0.5) \
        / newShape[2] * oldShape[2] - 0.5
    
    # # for debug purposes
    # coeffs = []
    
    base = np.floor(oldIndices).astype(int)
    diff = oldIndices - base
    result = np.zeros(newShape, dtype=np.float32)
    for i in range(2):
        coeffx = np.abs(i - (1 - diff[0, :, :]))
        coordx = base[0, :, :, :] + i
        flag = np.logical_and(coordx >= 0, coordx < oldShape[0])
        coeffx = coeffx * flag
        coordx[coordx<0] = 0
        coordx[coordx>=oldShape[0]] = oldShape[0] - 1
        for j in range(2):
            coeffy = np.abs(j - (1 - diff[1, :, :]))
            coordy = base[1, :, :, :] + j
            flag = np.logical_and(coordy >= 0, coordy < oldShape[1])
            coeffy = coeffy * flag
            coordy[coordy<0] = 0
            coordy[coordy>=oldShape[1]] = oldShape[1] - 1
            for k in range(2):
                coeffz = np.abs(k - (1 - diff[2, :, :]))
                coordz = base[2, :, :, :] + k
                flag = np.logical_and(coordz >= 0, coordz < oldShape[2])
                coeffz = coeffz * flag
                coordz[coordz<0] = 0
                coordz[coordz>=oldShape[2]] = oldShape[2] - 1
                coeff = coeffx * coeffy * coeffz
                result += input[coordx, coordy, coordz] * coeff
                # # for debug purposes
                # coeffs.append(coeff)
    # # for debug purposes
    # coeffs = [np.expand_dims(a, axis=0) for a in coeffs]
    # coeffs = np.concatenate(coeffs, axis=0)
    # outFile = '/data/datasets/UCLAPatients/visNew/temp/coeffs.npy'
    # np.save(outFile, coeffs)

    return result


def viewDose():
    """
    This function compares the resized dose and original dose
    """

    # # before the goal above, we firstly take a look at the coefficients matrix
    # file = '/data/datasets/UCLAPatients/visNew/temp/coeffs.npy'
    # coeffs = np.load(file)
    # reduce = np.sum(coeffs, axis=0)
    # return

    globalFolder = '/data/datasets/UCLAPatients'
    dataFolder = os.path.join(globalFolder, 'anonymousDataNew')
    expFolder = os.path.join(globalFolder, 'experiment')
    numPatients = 8

    outFolder = os.path.join(globalFolder, 'visNew', 'temp')
    if not os.path.isdir(outFolder):
        os.mkdir(outFolder)

    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimize')

        # find trailNO
        BOOfile, polishFile = findResultPath(optFolder)
        trailNO = 1
        while True:
            if str(trailNO) in BOOfile:
                break
            else:
                trailNO += 1

        polishData = loadmat(polishFile)['polishResult'][0, 0]
        polishDose = polishData['dose']
        reshapeFile = os.path.join(optFolder, 'dose{}.npy'.format(trailNO))
        reshapeData = np.load(reshapeFile)
        print(polishDose.shape, reshapeData.shape)

        # polishProfile = np.mean(polishDose, axis=(0, 1))
        # reshapeProfile = np.mean(reshapeData, axis=(0, 1))
        polishProfile = np.mean(polishDose, axis=(0, 2))
        reshapeProfile = np.mean(reshapeData, axis=(0, 2))
        plt.plot(polishProfile)
        plt.plot(reshapeProfile)
        plt.legend(['polish', 'reshape'])
        outFile = os.path.join(outFolder, '{}.png'.format(patientName))
        plt.savefig(outFile)
        plt.clf()
        break


def DVHCompClinicalOptimize():
    """
    This function compares the DVHs for clinical plans and MATLAB optimized plans
    """
    # this function runs on shenggpu4

    globalFolder = '/data/datasets/UCLAPatients'
    dataFolder = os.path.join(globalFolder, 'anonymousDataNew')
    expFolder = os.path.join(globalFolder, 'experiment')
    visFolder = os.path.join(globalFolder, 'visNew')
    numPatients = 8

    for idx in range(numPatients):
    # for idx in range(6, numPatients):
        patientName = 'patient{}'.format(idx+1)

        # load clinical dose
        # We state that the array shape should be of (height, width, slice)
        dataPatFolder = os.path.join(dataFolder, patientName)
        clinicalDoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        clinicalDose = pydicom.dcmread(clinicalDoseFile)
        clinicalDose = clinicalDose.pixel_array
        clinicalDose = np.transpose(clinicalDose, (1, 2, 0))
        # print(clinicalDose.shape)

        # load BOO dose. Firstly, find the correct trailNO
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimize')
        doseFiles = glob.glob(os.path.join(optFolder, 'dose*.npy'))
        doseFiles.sort()
        doseFile = doseFiles[-1]
        BOOdose = np.load(doseFile)

        last = doseFile.split('/')[-1]
        digit = ""
        for c in last:
            if c.isdigit():
                digit = digit + c
        trailNO = int(digit)
        # print(trailNO)
        # print(BOOdose.shape, '\n')
        # continue

        # load structures
        structuresFile = os.path.join(expPatFolder, 'structures.json')
        with open(structuresFile, 'r') as f:
            structures = json.load(f)
        PTVname = structures['ptv']
        OARnames = structures['oar']
        ROIs = [PTVname, ] + OARnames[1:]  # to exclude skin/body
        # print(ROIs)

        # load masks
        MRFolder = os.path.join(dataPatFolder, 'MR')
        rtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=MRFolder, rt_struct_path=rtFile)
        masks = {name: rtstruct.get_roi_mask_by_name(name) for name in ROIs}
        # print(masks[PTVname].shape)

        visPatFolder = os.path.join(visFolder, patientName)
        outFile = os.path.join(visPatFolder, 'DVHcomp{}.png'.format(trailNO))
        DVHgroup(BOOdose, masks)
        plt.legend(ROIs)
        plt.savefig(outFile)
        plt.clf()
        print(patientName)


def DVHgroup(dose, masks):
    for mask, array in masks.items():
        doseMasked = dose[array]
        doseMasked = np.sort(doseMasked)
        doseMasked = np.insert(doseMasked, 0, 0)
        yAxis = 1 - np.arange(len(doseMasked)) / len(doseMasked)
        plt.plot(doseMasked, yAxis)


def findResultPath(path):
    """
    This function finds the path to the most updated results.
    file names are in the form 'BOOresult${trailNO}.mat' and
    'polishResult${trailNO}.mat'. We find the files with the
    largest ${trailNO}
    """
    BOOFileTemplate = os.path.join(path, 'BOOresult{}.mat')
    BOOFile = None
    polishFileTemplate = os.path.join(path, 'polishResult{}.mat')
    polishFile = None
    count = 1
    while True:
        if os.path.isfile(BOOFileTemplate.format(count)) \
            and os.path.isfile(polishFileTemplate.format(count)):
            BOOFile = BOOFileTemplate.format(count)
            polishFile = polishFileTemplate.format(count)
            count += 1
        else:
            break
    return BOOFile, polishFile


def getPTVname(expPatFolder):
    structuresFile = os.path.join(expPatFolder, 'structures.json')
    with open(structuresFile, 'r') as f:
        data = json.load(f)
        PTVname = data['ptv']
    return PTVname


if __name__ == '__main__':
    # examineMask()
    # examineDose()
    # draw_DVH()
    # doseResize()
    # viewDose()
    DVHCompClinicalOptimize()