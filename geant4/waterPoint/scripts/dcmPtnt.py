import os
import numpy as np
import json
import matplotlib.pyplot as plt
import pydicom
from rt_utils import RTStructBuilder
import h5py
from PIL import Image
import pdb
import pickle

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


def doseRead():
    """
    I ran the beamlet dose calculation code on shenggpu6 bug-free. 
    This function is to take a look at the dose
    """
    expFolder = os.path.join(rootFolder, 'patient1_dosecalc')
    coeffFile = os.path.join(expFolder, 'Dose_Coefficients.h5')
    # doseShape = (160, 160, 158)
    doseShape = (158, 160, 160)

    Dataset = h5py.File(coeffFile, 'r')
    # print_hdf5_structure(Dataset)
    Dataset_Beams = Dataset['beams']
    Dataset_Beams_Data = Dataset_Beams['data']
    beam_data = Dataset_Beams_Data['beam_00000']
    
    nBeamlets = 400
    globalMat = np.zeros(doseShape, dtype=np.float32)
    globalMat = globalMat.flatten()
    for i in range(nBeamlets):
        beamletName = 'beamlet_{:05d}'.format(i)
        beamlet = beam_data[beamletName]
        coeffs = beamlet['coeffs'][()]
        lindex = beamlet['lindex'][()]
        # I guess the coeffs and lindex represent sparse data. Let's take a look
        localMat = np.zeros_like(globalMat)
        localMat = localMat.flatten()
        sparseToDense(coeffs, lindex, localMat)
        globalMat += localMat
        print('beamlet {}'.format(i))
    # reshape it back
    globalMat = np.reshape(globalMat, doseShape)
    file = os.path.join(expFolder, 'Dose_Coefficients.npy')
    np.save(file, globalMat)


def sparseToDense(coeffs, lindex, localMat):
    """
    This function converts the sparse matrix into dense matrix
    """
    assert len(coeffs) == len(lindex), "The sparse matrix does not align in size"
    for idx, value in zip(lindex, coeffs):
        localMat[idx] = value


def print_hdf5_structure(group, level=0):
    """Recursively print the hierarchy of groups and datasets."""
    for key in group.keys():
        if isinstance(group[key], h5py.Group):
            print("  " * level + f"Group: {key}/")
            print_hdf5_structure(group[key], level + 1)
        elif isinstance(group[key], h5py.Dataset):
            print("  " * level + f"Dataset: {key}")


def visDose():
    """
    This function shows the dose 
    """
    expFolder = os.path.join(rootFolder, 'patient1_dosecalc')
    dataFile = os.path.join(expFolder, 'Dose_Coefficients.npy')
    array = np.load(dataFile)

    DoseVisFolder = os.path.join(rootFolder, 'DoseVis')
    if not os.path.isdir(DoseVisFolder):
        os.mkdir(DoseVisFolder)
    
    # Normalize
    roof = np.max(array)
    array /= roof

    for i in range(array.shape[0]):
        slice = array[i, :, :]
        FileName = os.path.join(DoseVisFolder, 'slice{}.png'.format(i+1))
        plt.imsave(FileName, slice, vmin=0, vmax=1)


def doseOverlap():
    """
    In the function above, we showed that 
    """
    doseSliceNo = 98
    # convert dose slice number to CT slice number
    CTSliceNo = 198 - int(doseSliceNo * 197 / 158)
    
    dicomFolder = os.path.join(rootFolder, 'patient1_dicom')
    dicomFile = os.path.join(dicomFolder, '{:03d}.dcm'.format(CTSliceNo))
    DoseVisFolder = os.path.join(rootFolder, 'DoseVis')
    doseSliceImg = os.path.join(DoseVisFolder, 'slice{}.png'.format(doseSliceNo))

    # then we overlay the two images
    dcmImg = pydicom.dcmread(dicomFile).pixel_array
    dcmImgRoof = np.max(dcmImg)
    uint8Max = 255
    dcmImg = np.uint8(dcmImg / dcmImgRoof * uint8Max)

    doseImg = Image.open(doseSliceImg)
    doseImg = doseImg.resize(dcmImg.shape)
    doseImg = np.array(doseImg)

    print(dcmImg.shape, dcmImg.dtype, doseImg.shape, doseImg.dtype)
    result = np.zeros_like(doseImg)
    weight = 0.3
    result[:, :, 0] = weight * doseImg[:, :, 0] + (1 - weight) * dcmImg
    result[:, :, 1] = weight * doseImg[:, :, 1] + (1 - weight) * dcmImg
    result[:, :, 2] = weight * doseImg[:, :, 2] + (1 - weight) * dcmImg
    result[:, :, 3] = uint8Max
    resultPIL = Image.fromarray(result)
    file = os.path.join(rootFolder, 'patient1_dosecalc', 'patDose.png')
    resultPIL.save(file)


def slabPhantomReconstruction():
    """
    This function creates a slab phantom the same 
    as the one shown in Dr. Neph's CCCS paper
    """
    if False:
        # firstly, have a dicom to take a look
        dcmFile80 = os.path.join(rootFolder, 'patient1_dicom', '080.dcm')
        data80 = pydicom.dcmread(dcmFile80)
        dcmFile81 = os.path.join(rootFolder, 'patient1_dicom', '081.dcm')
        data81 = pydicom.dcmread(dcmFile81)

        data80_attributes = dir(data80)
        data81_attributes = dir(data81)
        assert data80_attributes == data81_attributes, "The two object has different attributes"

    # get a template and modify
    templateFile = os.path.join(rootFolder, 'patient1_dicom', '001.dcm')
    templateData = pydicom.dcmread(templateFile)
    attributes = dir(templateData)
    
    Rows = 256
    Columns = 256
    Slices = 256

    # edit slice
    slice = np.ones((Rows, Columns), dtype=np.uint16) * 1000
    slice[:16, :] = 920  # 1.6 cm of adipose, 0.92 g/cm^3
    slice[16:32, :] = 1040  # 1.6 cm of muscle, 1.04 g/cm^3
    slice[32:48, :] = 1850  # 1.6 cm of bone, 1.85 g/cm^3
    slice[48:64, :] = 1040  # 1.6 cm of muscle, 1.04 g/cm^3
    slice[64:160, :] = 250  # 9.6 cm of lung, 0.25 g/cm^3
    slice[160:176, :] = 1040  # 1.6 cm of muscle, 1.04 g/cm^3
    slice[176:192, :] = 1850  # 1.6 cm of bone, 1.85 g/cm^3
    slice[192:208, :] = 920  # 1.6 cm of adipose, 0.92 g/cm^3
    slice[208:224, :] = 1850  # 1.6 cm of bone, 1.85 g/cm^3
    slice[224:240, :] = 1040  # 1.6 cm of muscle, 1.04 g/cm^3
    slice[240:256, :] = 920  # 1.6 cm of adipose, 0.92 g/cm^3

    SOPInstanceUIDtemplate = getattr(templateData, 'SOPInstanceUID')
    SOPInstanceUIDtemplate = SOPInstanceUIDtemplate.split('.')
    SOPInstanceUIDtemplate[-1] = '{}'
    SOPInstanceUIDtemplate = '.'.join(SOPInstanceUIDtemplate)
    templateData.Rows = Rows
    templateData.Columns = Columns
    templateData.PixelSpacing = [1.0, 1.0]
    templateData.SliceThickness = '1.0'
    templateData.PixelData = slice.tobytes()

    # print(templateData)
    # return

    if False:
        image = templateData.pixel_array
        file = os.path.join(rootFolder, 'patient1_dosecalc', 'water.png')
        plt.imsave(file, image, cmap='gray')
    
    slabPhantomFolder = os.path.join(rootFolder, 'slab_dicom')
    if not os.path.isdir(slabPhantomFolder):
        os.mkdir(slabPhantomFolder)
    
    for i in range(Slices):
        idx = i + 1
        SOPInstanceUID = SOPInstanceUIDtemplate.format(idx)
        InstanceNumber = str(idx)
        offset = i + 0.5
        ImagePosition = [0., 0., -offset]
        SliceLocation = str(offset)

        templateData.SOPInstanceUID = SOPInstanceUID
        templateData.InstanceNumber = InstanceNumber
        templateData.ImagePositionPatient = ImagePosition
        templateData.SliceLocation = SliceLocation

        dcmFile = os.path.join(slabPhantomFolder, '{:03d}.dcm'.format(idx))
        templateData.save_as(dcmFile)
        # for debug purposes
        if (i == 2):
            print(templateData)


def RTcreate():
    """
    This function is to create an RTstruct file for the 
    slab phantom created in the function above
    """
    slabPhantomFolder = os.path.join(rootFolder, 'slab_dicom')
    rtstruct = RTStructBuilder.create_new(slabPhantomFolder)
    maskShape = (256, 256, 256)
    mask = np.ones(maskShape, dtype=bool)
    rtstruct.add_roi(mask=mask, color=[255, 0, 255], name="water")
    rtFile = os.path.join(slabPhantomFolder, 'RTstructure.dcm')
    rtstruct.save(rtFile)


def RTexamine():
    """
    A simple function to examine the rtstructures
    """
    slabPhantomFolder = os.path.join(rootFolder, 'slab_dicom')
    rtFile = os.path.join(slabPhantomFolder, 'RTstructure.dcm')
    rtstruct = RTStructBuilder.create_from(slabPhantomFolder, rtFile)
    names = rtstruct.get_roi_names()
    mask = rtstruct.get_roi_mask_by_name(names[0])
    breakpoint()


def CCCSexamine():
    """
    We ran the dose calculation code bug-free. 
    Now we take a look at the dose
    """
    slabExpFolder = os.path.join(rootFolder, 'slab_dosecalc')
    fluenceDim = 20

    shape = (256, 256, 256)
    if True:
        # read the data
        outputFile = os.path.join(slabExpFolder, 'Dose_Coefficients.h5')
        dataset = h5load(outputFile, shape)

    # for visualization
    if False:
        slabVisFolder = os.path.join(rootFolder, 'SlabDoseVis')
        if not os.path.isdir(slabVisFolder):
            os.mkdir(slabVisFolder)

        beamKey = list(dataset.keys())[0]
        beamlets = dataset[beamKey]
        beamlets_keys = list(beamlets.keys())
        beamlets_keys.sort()
        for i, key in enumerate(beamlets_keys):
            beamletMat = beamlets[key]
            slice = beamletMat[:, 128, :]
            file = os.path.join(slabVisFolder, '{:03d}.png'.format(i+1))
            plt.imsave(file, slice)
            print(i)
    
    if True:
        # we got to know that the beamlets are arrange in a certaim manner. 
        # So the central beam index should be 61 (60 starting from 0)
        beamKey = list(dataset.keys())[0]
        beamlets = dataset[beamKey]

        centralIdx = 40
        centralBeamletKey = 'beamlet_{:05d}'.format(centralIdx)
        centralBeamlet = beamlets[centralBeamletKey]
        centerlineFile = os.path.join(slabExpFolder, 'centerline.npy')
        np.save(centerlineFile, centralBeamlet)
    
    if True:
        centerlineFile = os.path.join(slabExpFolder, 'centerline.npy')
        array = np.load(centerlineFile)

        # centerline cross-section
        crossSection = array[128]
        image = os.path.join(slabExpFolder, 'centerline.png')
        plt.imsave(image, crossSection)

        if False:
            # depth-dose cerve
            examineRange = 2
            depth = np.arange(256) * 0.1  # cm
            labels = []
            for i in range(128-examineRange,128+examineRange):
                for j in range(128-examineRange, 128+examineRange):
                    line = array[i, :, j]
                    plt.plot(depth, line)
                    labels.append((i, j))
            plt.legend(labels)
            file = os.path.join(slabExpFolder, 'depthDoseCurve.png')
            plt.savefig(file)
            plt.clf()
        
        depth = np.arange(256) * 0.1  # cm
        line = array[129, :, 126]
        plt.plot(depth, line)
        file = os.path.join(slabExpFolder, 'depthDoseCurveTentative.png')
        plt.savefig(file)
        plt.clf()


def h5load(file, shape):
    """
    This function reads the matrices from file
    """
    matSize = 1
    for a in shape:
        matSize *= a

    result = {}
    dataset = h5py.File(file, 'r')
    beams = dataset['beams']['data']
    beams_keys = list(beams.keys())
    beams_keys.sort()
    for beamKey in beams_keys:
        beamDataset = beams[beamKey]
        beamlets_keys = list(beamDataset.keys())
        beamlets_keys.sort()
        result[beamKey] = {}
        for beamlet_key in beamlets_keys:
            coeffs = beamDataset[beamlet_key]['coeffs'][()]
            lindex = beamDataset[beamlet_key]['lindex'][()]
            element = np.zeros(matSize, dtype=coeffs.dtype)
            for a, b in zip(lindex, coeffs):
                element[a] = b
            element = np.reshape(element, shape)
            result[beamKey][beamlet_key] = element
            print(beamlet_key)
    return result


def CCCSBatchStudy():
    """
    This function carries a batch study for the CCCS dose calculation,
    by varing the fluence dimension and the beamlet size
    """
    FluenceDim = 9
    BeamletSize = 2.0
    shape = (256, 256, 256)
    size = shape[0] * shape[1] * shape[2]

    expFolder = os.path.join(rootFolder, 'slab_dosecalc_{}_{}'.format(FluenceDim, BeamletSize))
    Dose_Coefficients_file = os.path.join(expFolder, 'Dose_Coefficients.h5')
    dataset = h5py.File(Dose_Coefficients_file, 'r')
    dataset = dataset['beams']['data']
    beam = dataset['beam_00000']
    centralBeamletIdx = int((FluenceDim - 1) / 2 * FluenceDim + (FluenceDim - 1) / 2)
    BeamletKey = 'beamlet_{:05d}'.format(centralBeamletIdx)
    beamlet = beam[BeamletKey]
    coeffs = beamlet['coeffs'][()]
    lindex = beamlet['lindex'][()]

    DoseMat = np.zeros(size, dtype=coeffs.dtype)
    for key, value in zip(lindex, coeffs):
        DoseMat[key] = value
    DoseMat = np.reshape(DoseMat, shape)

    # find the location for the beamlet by calculating the mass center
    slice = DoseMat[:, 128, :]
    result1 = 0.
    for i in range(slice.shape[0]):
        result1 += i * np.sum(slice[i, :])
    result1 /= np.sum(slice)
    result2 = 0.
    for i in range(slice.shape[1]):
        result2 += i * np.sum(slice[:, i])
    result2 /= np.sum(slice)
    coords = (result1, result2)
    print(coords)

    if False:
        sliceLowerIdx = int(np.floor(result1))
        sliceHigherIdx = int(np.ceil(result1))

        sliceLower = DoseMat[sliceLowerIdx, :, :]
        sliceLowerFile = os.path.join(expFolder, 'sliceLower.png')
        plt.imsave(sliceLowerFile, sliceLower)

        sliceHigher = DoseMat[sliceHigherIdx, :, :]
        sliceHigherFile = os.path.join(expFolder, 'sliceHigher.png')
        plt.imsave(sliceHigherFile, sliceHigher)
    
    coords00 = (int(np.floor(coords[0])), int(np.floor(coords[1])))
    coords01 = (int(np.floor(coords[0])), int(np.ceil(coords[1])))
    coords10 = (int(np.ceil(coords[0])), int(np.floor(coords[1])))
    coords11 = (int(np.ceil(coords[0])), int(np.ceil(coords[1])))
    depth = np.arange(shape[1]) * 0.1  # cm
    line00 = DoseMat[coords00[0], :, coords00[1]]
    line01 = DoseMat[coords01[0], :, coords01[1]]
    line10 = DoseMat[coords10[0], :, coords10[1]]
    line11 = DoseMat[coords11[0], :, coords11[1]]
    plt.plot(depth, line00)
    plt.plot(depth, line01)
    plt.plot(depth, line10)
    plt.plot(depth, line11)
    plt.legend(['line00', 'line01', 'line10', 'line11'])
    plt.xlabel('depth (cm)')
    plt.ylabel('dose a.u.')
    plt.title('beamlet ddp with fdim {}, size {}'.format(FluenceDim, BeamletSize))
    figureFile = os.path.join(expFolder, 'depthDose_{}_{}.png'.format(FluenceDim, BeamletSize))
    plt.savefig(figureFile)
    plt.clf()


if __name__ == '__main__':
    # viewDcm()
    # structPrepare()
    # doseRead()
    # visDose()
    # doseOverlap()
    # slabPhantomReconstruction()
    # RTcreate()
    # RTexamine()
    # CCCSexamine()
    CCCSBatchStudy()