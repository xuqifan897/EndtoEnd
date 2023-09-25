import os
import numpy as np
import json
import matplotlib.pyplot as plt
import pydicom
from rt_utils import RTStructBuilder

rootFolder = "/data/qifan/projects/EndtoEnd/results/slabBench"

def viewDcm():
    """
    In this function, we view the dicom images and rt structures
    """
    patDicomFolder = os.path.join(rootFolder, 'patient1_dicom')
    patImgFolder = os.path.join(rootFolder, 'patient1_image')
    if not os.path.isdir(patImgFolder):
        os.mkdir(patImgFolder)

    if False:
        files = os.listdir(patDicomFolder)
        for file in files:
            if 'RTstructure' in file:
                continue
            fileFull = os.path.join(patDicomFolder, file)
            dcmData = pydicom.dcmread(fileFull)
            dcmArray = dcmData.pixel_array

            fileName = file.split('.')[0]
            figureName = fileName + '.png'
            figureNameFull = os.path.join(patImgFolder, figureName)
            plt.imsave(figureNameFull, dcmArray, cmap='gray')
    
    rtFile = os.path.join(patDicomFolder, 'RTstructure.dcm')
    rtstruct = RTStructBuilder.create_from(
        dicom_series_path=patDicomFolder, 
        rt_struct_path=rtFile)
    
    # get PTV mask
    maskName = 'BODY'
    RTstructNames = rtstruct.get_roi_names()
    assert maskName in RTstructNames, "PTV is not in rtstruct"
    mask = rtstruct.get_roi_mask_by_name(maskName)
    mask = np.flip(mask, axis=2)
    maskFolder = os.path.join(rootFolder, 'bodyMask')
    if not os.path.isdir(maskFolder):
        os.mkdir(maskFolder)
    scale = 255
    for i in range(mask.shape[2]):
        slice = mask[:, :, i]
        slice = np.uint8(slice * scale)
        imgFile = os.path.join(maskFolder, 'slice{}.png'.format(i+1))
        plt.imsave(imgFile, slice)

    if False:
        rtData = pydicom.read_file(rtFile)
        ctrs = rtData.ROIContourSequence
        print('great!')


def structPrepare():
    """Prepare the RTstructure json file"""
    bodyName = 'BODY'
    struct = {"prescription": 20,
        "ptv": bodyName,
        'oar': []}
    
    patDicomFolder = os.path.join(rootFolder, 'patient1_dicom')
    rtFile = os.path.join(patDicomFolder, 'RTstructure.dcm')
    rtstruct = RTStructBuilder.create_from(
        dicom_series_path=patDicomFolder, 
        rt_struct_path=rtFile)
    RTstructNames = rtstruct.get_roi_names()
    struct['oar'] = [a for a in RTstructNames if a != bodyName]
    structPath = '/data/qifan/projects/EndtoEnd/results' \
        '/slabBench/patien1_dosecalc'
    jsonFile = os.path.join(structPath, 'structures.json')
    with open(jsonFile, 'w') as f:
        json.dump(struct, f)


if __name__ == '__main__':
    # viewDcm()
    structPrepare()