import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pydicom
import cv2
from rt_utils import RTStructBuilder

# this script runs on shenggpu4. Here we specify the path to the data
globalFolder = '/data/datasets/UCLAPatients'
patients = ['0530793', '0678774', '0999272', '2415674', '4135918', '5758498', '6865205']
# patients = ['0530793', '0678774', '2415674', '4135918', '5758498', '6865205']
num_patients = 7
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 
    'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

def examineMR():
    """
    this function anonymizes the patient data
    """
    # create the folder for the anonymized data
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    if not os.path.exists(anonymousDataPath):
        os.makedirs(anonymousDataPath)

    # then we examine the individual dicom files
    # the *MRdose* directory has 6 subdirectories:
    # *MR*, *RTDOSE*, *RTPLAN*, *RTst*00000, *RTst*00001, *SEG*
    # the *CT* directory has 2 subdirectories:
    # *CT*, *RTst*
    rawDataPath = os.path.join(globalFolder, 'rawData')
    for patient in patients:
        patientFolder = os.path.join(rawDataPath, patient)
        MRdoseTemplate = os.path.join(patientFolder, '*MRdose')
        MRdoseFolder = folderEvaluate(MRdoseTemplate)

        MRTemplate = os.path.join(MRdoseFolder, '*MR*')
        RTDOSETemplate = os.path.join(MRdoseFolder, '*RTDOSE*')
        RTPLANTemplate = os.path.join(MRdoseFolder, '*RTPLAN*')
        RTst0Template = os.path.join(MRdoseFolder, '*RTst*00000')
        RTst1Template = os.path.join(MRdoseFolder, '*RTst*00001')
        SEGTemplate = os.path.join(MRdoseFolder, '*SEG*')

        MRfolder = folderEvaluate(MRTemplate)
        RTDOSEfolder = folderEvaluate(RTDOSETemplate)
        RTPLANfolder = folderEvaluate(RTPLANTemplate)
        RTst0folder = folderEvaluate(RTst0Template)
        RTst1folder = folderEvaluate(RTst1Template)
        SEGfolder = folderEvaluate(SEGTemplate)

        # examine MR dicoms
        # MR files: 
        # InstitutionName, Manufacturer, ManufacturerModelName, 
        # PatientBirthDate, PatientID, PatientName, PatientSex
        #
        # RTDOSE is the same as MR files. It has an attribute, 
        # pixel_data, which contains dose information
        #
        # RTst:
        # InstitutionName, Manufacturer, ManufacturerModelName, 
        # OperatorName, PatientBirthDate, PatientID, PatientName, 
        # PatientSex, PerformingPhysicianName

        files = os.listdir(MRfolder)
        print("There are {} dicom files in the MR folder of patient".\
            format(len(files), patient))
        files.sort()
        file = os.path.join(MRfolder, files[0])
        dicomData = pydicom.dcmread(file)

        # examine RTDOSE
        files = os.listdir(RTDOSEfolder)
        file = os.path.join(RTDOSEfolder, files[0])
        dicomData = pydicom.dcmread(file)

        # examine RTPLAN
        files = os.listdir(RTPLANfolder)
        file = os.path.join(RTPLANfolder, files[0])
        dicomData = pydicom.dcmread(file)

        # examine RTst*00000
        files = os.listdir(RTst0folder)
        file = os.path.join(RTst0folder, files[0])
        dicomData = pydicom.dcmread(file)

        # examine RTst*00001
        files = os.listdir(RTst1folder)
        file = os.path.join(RTst1folder, files[0])
        dicomData = pydicom.dcmread(file)

        # examine SEG
        files = os.listdir(SEGfolder)
        file = os.path.join(SEGfolder, file)
        dicomData = pydicom.dcmread(file)

        break


def anonymizeMR():
    rawDataPath = os.path.join(globalFolder, 'rawData')

    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    if not os.path.exists(anonymousDataPath):
        os.makedirs(anonymousDataPath)
    
    # attributes = ['InstitutionName', 'Manufacturer', 'ManufacturerModelName', \
    #     'OperatorName', 'PatientBirthDate', 'PatientID', 'PatientName', \
    #     'PatientSex', 'PerformingPhysicianName']

    for idx, patient in enumerate(patients):
        patientNameNew = 'patient{}'.format(idx+1)
        patFolderOld = os.path.join(rawDataPath, patient)
        patFolderNew = os.path.join(anonymousDataPath, patientNameNew)

        MRdoseOld = os.path.join(patFolderOld, '*MRdose')
        MRdoseOld = folderEvaluate(MRdoseOld)
        MRdoseNew = os.path.join(patFolderNew, 'MRdose')
        assert os.path.isdir(MRdoseOld)
        if not os.path.isdir(MRdoseNew):
            os.makedirs(MRdoseNew)

        # anomymizing MR dicoms
        MRFolderOld = os.path.join(MRdoseOld, '*MR*')
        MRFolderOld = folderEvaluate(MRFolderOld)
        MRFolderNew = os.path.join(MRdoseNew, 'MR')
        if not os.path.isdir(MRFolderNew):
            os.mkdir(MRFolderNew)
        files = os.listdir(MRFolderOld)
        files.sort()
        for file in files:
            fileOld = os.path.join(MRFolderOld, file)
            dicomData = pydicom.dcmread(fileOld)
            wash(dicomData)
            fileNew = os.path.join(MRFolderNew, file)
            dicomData.save_as(fileNew)
        
        # anonymizing RTDOSE dicom
        RTDOSEFolderOld = os.path.join(MRdoseOld, '*RTDOSE*')
        RTDOSEFolderOld = folderEvaluate(RTDOSEFolderOld)
        files = os.listdir(RTDOSEFolderOld)
        assert len(files) == 1
        file = files[0]
        fileOld = os.path.join(RTDOSEFolderOld, file)
        dicomData = pydicom.dcmread(fileOld)
        wash(dicomData)
        fileNew = os.path.join(MRdoseNew, 'RTdose.dcm')
        dicomData.save_as(fileNew)

        # anonymizing RTPLAN dicom
        RTPLANFolderOld = os.path.join(MRdoseOld, '*RTPLAN*')
        RTPLANFolderOld = folderEvaluate(RTPLANFolderOld)
        files = os.listdir(RTPLANFolderOld)
        assert len(files) == 1
        file = files[0]
        fileOld = os.path.join(RTPLANFolderOld, file)
        dicomData = pydicom.dcmread(fileOld)
        wash(dicomData)
        fileNew = os.path.join(MRdoseNew, 'RTPLAN.dcm')
        dicomData.save_as(fileNew)

        # anonymizing RTstruct
        for i in range(2):
            RTstFolderOld = os.path.join(MRdoseOld, '*RTst*0000{}'.format(i))
            RTstFolderOld = folderEvaluate(RTstFolderOld)
            files = os.listdir(RTstFolderOld)
            assert len(files) == 1
            file = files[0]
            fileOld = os.path.join(RTstFolderOld, file)
            dicomData = pydicom.dcmread(fileOld)
            wash(dicomData)
            fileNew = os.path.join(MRdoseNew, 'RTst{}.dcm'.format(i))
            dicomData.save_as(fileNew)
        
        # anonymizing SEG
        SEGFolderOld = os.path.join(MRdoseOld, '*SEG*')
        SEGFolderOld = folderEvaluate(SEGFolderOld)
        files = os.listdir(SEGFolderOld)
        file = files[0]
        fileOld = os.path.join(SEGFolderOld, file)
        dicomData = pydicom.dcmread(fileOld)
        wash(dicomData)
        fileNew = os.path.join(MRdoseNew, 'SEG.dcm')
        dicomData.save_as(fileNew)

        print('{}: {}'.format(patientNameNew, patient))

        

def folderEvaluate(template):
    results = glob.glob(template)
    assert len(results) == 1, \
        "{} matches 0 or more than 1 file/folders.".format(template)
    return results[0]


def wash(dicomData):
    attributes = ['InstitutionName', 'Manufacturer', 'ManufacturerModelName', \
        'OperatorsName', 'PatientBirthDate', 'PatientID', 'PatientName', \
        'PatientSex', 'PerformingPhysicianName']
    for attribute in attributes:
        if hasattr(dicomData, attribute):
            setattr(dicomData, attribute, '')
    if hasattr(dicomData, 'file_meta'):
        if hasattr(dicomData.file_meta, 'ImplementationVersionName'):
            dicomData.file_meta.ImplementationVersionName = ''
        if hasattr(dicomData.file_meta, 'SourceApplicationEntityTitle'):
            dicomData.file_meta.SourceApplicationEntityTitle = ''


def visualizeMR():
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    if not os.path.isdir(visFolder):
        os.mkdir(visFolder)

    
    roof = 255
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(anonymousDataPath, patientName)
        MRFolder = os.path.join(patFolder, 'MRdose', 'MR')
        files = dicomSort(MRFolder, 'MR')
        
        visPatFolder = os.path.join(visFolder, patientName)
        visMRFolder = os.path.join(visPatFolder, 'MR')
        if not os.path.isdir(visMRFolder):
            os.makedirs(visMRFolder)
        
        for j, file in enumerate(files):
            oldFile = os.path.join(MRFolder, file)
            newFile = os.path.join(visMRFolder, '{:03d}.png'.format(j+1))
            dicomData = pydicom.dcmread(oldFile)
            pixel_array = dicomData.pixel_array
            roof = np.max(pixel_array)
            pixel_array = np.uint8(pixel_array / roof * 255)
            plt.imsave(newFile, pixel_array, cmap='gray')
        print(patientName)


def dicomSort(folder, modality):
    files = os.listdir(folder)
    result = []
    for file in files:
        file_ = os.path.join(folder, file)
        dicomData = pydicom.dcmread(file_)
        assert dicomData.Modality == modality
        result.append((file, int(dicomData.InstanceNumber)))
    result.sort(key = lambda x: x[1])
    result = [a[0] for a in result]
    
    return result


def visualizeRTPLAN():
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
        file = os.path.join(MRdoseFolder, 'RTPLAN.dcm')
        RTplan = pydicom.dcmread(file)
        print('number of beams: {}'.format(len(RTplan.BeamSequence)))


def visualizeRTdose():
    """
    This function visualizes RT dose.
    The data matrix of RT dose is of the shape (slice, hight, width),
    the slice dimension is reversed from InstanceNumber
    """
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
        RTdosefile = os.path.join(MRdoseFolder, 'RTdose.dcm')
        dose = pydicom.dcmread(RTdosefile).pixel_array
        roof = np.max(dose)

        visPatFolder = os.path.join(visFolder, patientName)
        visDoseFolder = os.path.join(visPatFolder, 'dose')
        if not os.path.isdir(visDoseFolder):
            os.mkdir(visDoseFolder)
        
        for i in range(dose.shape[0]):
            outputFile = os.path.join(visDoseFolder, '{:03d}.png'.format(i+1))
            plt.imsave(outputFile, dose[dose.shape[0]-1-i, :, :], vmax=roof)
        print(patientName)


def visualizeRTstruct():
    """
    This function visualizes the RT structures of different patients,
    each patient has three files that are of RTSTRUCT modality, they
    are: RTst0.dcm, RTst1.dcm, and SEG.dcm

    Result:
    we could not open SEG.dcm using RTStructBuilder package
    RTst0 and RTst1 are the same (at least visually the same)
    the shape of RTSTRUCT mask is (hight, width, slice)
    the slice dimension is reversed from InstanceNumber
    """
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    
    scale = 255
    RTfiles = ['RTst0.dcm', 'RTst1.dcm']
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
        MRFolder = os.path.join(MRdoseFolder, 'MR')
        for file in RTfiles:
            fullFile = os.path.join(MRdoseFolder, file)
            rtstruct = RTStructBuilder.create_from(
                dicom_series_path=MRFolder, rt_struct_path=fullFile)
            ROInames = rtstruct.get_roi_names()
            print('{} {}\n{}\n\n'.format(patientName, file, ROInames))
            masks = {}
            for name in ROInames:
                try:
                    masks[name] = rtstruct.get_roi_mask_by_name(name)
                except:
                    print('fail to extract the mask for {}'.format(name))
            
            # create folder
            visRTstructFolder = os.path.join(visFolder, patientName, file.split('.')[0])
            if not os.path.isdir(visRTstructFolder):
                os.makedirs(visRTstructFolder)
            
            # draw contours
            MRfiles = dicomSort(MRFolder, 'MR')
            for j, file in enumerate(MRfiles):
                MRfile = os.path.join(MRFolder, file)
                result = draw_mask(MRfile, masks, len(MRfiles) - 1 - j)
                outputPath = os.path.join(visRTstructFolder, '{:03d}.png'.format(j+1))
                cv2.imwrite(outputPath, result)
                print('{} {} {}'.format(patientName, file, j+1))


def draw_mask(MRfile, masks, frameID):
    pixel_array = pydicom.dcmread(MRfile).pixel_array
    pixel_array = gray2rgb(pixel_array)
    for name, mask in masks.items():
        maskFrame = mask[:, :, frameID]
        maskFrame = np.uint8(255 * maskFrame)
        contours, hierarchy = cv2.findContours(maskFrame, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        pixel_array = cv2.drawContours(pixel_array, contours, -1, (0, 255, 0), 1)
    # cv2.imshow('MR contour', pixel_array)
    # cv2.waitKey()
    return pixel_array



def gray2rgb(input):
    scale = 255
    roof = np.max(input)
    input = np.uint8(input / roof * scale)
    input = np.expand_dims(input, axis=2)
    input = np.repeat(input, 3, axis=2)
    return input


def test_cv2_drawContour():
    imagePath = '/data/datasets/UCLAPatients/visualize/patient1/MR/041.png'
    image = cv2.imread(imagePath)
    imgray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    ret, thresh = cv2.threshold(imgray, 128, 128, 128)
    contours, hierachy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    
    canvas = np.zeros_like(image)
    result = cv2.drawContours(canvas, contours, -1, (0, 255, 0), 3)
    cv2.imshow('result', result)
    cv2.waitKey(0)


def visualizeDVH():
    """
    This function creates the DVH histograms, and write them to files
    """
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    DVHFolder = os.path.join(visFolder, 'DVH')
    if not os.path.isdir(DVHFolder):
        os.mkdir(DVHFolder)
    

    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
        doseFile = os.path.join(MRdoseFolder, 'RTdose.dcm')
        dose = pydicom.dcmread(doseFile).pixel_array

        # flip to normal geometry, shape: (slice, height, width)
        dose = np.flip(dose, axis=0)

        RTSTRUCTFile = os.path.join(MRdoseFolder, 'RTst0.dcm')
        MRFolder = os.path.join(MRdoseFolder, 'MR')
        rtstruct = RTStructBuilder.create_from(MRFolder, RTSTRUCTFile)
        ROInames_unfiltered = rtstruct.get_roi_names()
        PTVuncropped = [a for a in ROInames_unfiltered if 'PTV' in a and 'uncrop' in a]
        ROInames = PTVuncropped + [a for a in ROInames_unfiltered if a[:2] == 'O_']
        masks = {}
        for name in ROInames:
            try:
                masks[name] = rtstruct.get_roi_mask_by_name(name)
            except:
                print('fail to extract the mask for {}'.format(name))
        # flip the masks
        for name, mask in masks.items():
            masks[name] = np.transpose(np.flip(mask, axis=2), (2, 0, 1))
        # after this, both dose and RTstruct mask are in shape (slice, height, width)

        # start to draw DVH
        count = 0
        for name, mask in masks.items():
            maskedDose = dose[mask]
            maskedDose = np.sort(maskedDose)
            size = maskedDose.size
            xAxis = np.zeros(size+1, dtype=maskedDose.dtype)
            xAxis[1:] = maskedDose
            yAxis = np.zeros(size+1, dtype=np.float32)
            yAxis[0] = 1
            yAxis[1:] = 1 - np.arange(size) / size
            plt.plot(xAxis, yAxis)
            count += 1
        plt.legend(masks.keys())
        plt.xlabel('dose (a.u.)')
        plt.ylabel('relative volume')
        plt.title(patientName)
        fileName = os.path.join(DVHFolder, patientName + '.png')
        plt.savefig(fileName)
        plt.clf()
        # plt.show()
        # break


def printAnatomy():
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    numPatients = 7
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
        MRFolder = os.path.join(MRdoseFolder, 'MR')
        RTSTRUCTFile = os.path.join(MRdoseFolder, 'RTst0.dcm')
        rtstruct = RTStructBuilder.create_from(MRFolder, RTSTRUCTFile)
        ROINames = rtstruct.get_roi_names()
        print(patientName)
        print(ROINames, '\n\n')


def DVHdraft():
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    
    # for i in range(num_patients):
    #     patientName = 'patient{}'.format(i+1)
    #     MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
    #     doseFile = os.path.join(MRdoseFolder, 'RTdose.dcm')
    #     dose = pydicom.dcmread(doseFile).pixel_array
    #     dose = np.reshape(dose, -1)
    #     dose = np.sort(dose)
    #     size = dose.size
    #     yAxis = 1 - np.arange(size) / size
    #     plt.plot(dose, yAxis)
    #     plt.show()
    #     # print(dose[-100:])
    #     break

    # # generate mask images
    # for i in range(num_patients):
    #     patientName = 'patient{}'.format(i+1)
    #     MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')

    #     RTSTRUCTFile = os.path.join(MRdoseFolder, 'RTst0.dcm')
    #     MRFolder = os.path.join(MRdoseFolder, 'MR')
    #     rtstruct = RTStructBuilder.create_from(MRFolder, RTSTRUCTFile)
    #     RprintAnatomy   try:
    #             masks[name] = rtstruct.get_roi_mask_by_name(name)
    #         except:
    #             print('fail to extract the mask for {}'.format(name))
        
    #     # if 'N_AIR' in masks:
    #     #     print('number of N_AIR voxels in {}: {}'.format(patientName, np.sum(masks['N_AIR'])))
        
    #     visRTFolder = os.path.join(visFolder, patientName, 'anatomy')
    #     for key, mask in masks.items():
    #         keyFolder = os.path.join(visRTFolder, key)
    #         if not os.path.isdir(keyFolder):
    #             os.makedirs(keyFolder)
            
    #         mask = np.uint8(255 * mask)
    #         # mask shape: (height, width, slice), slice is reversed
    #         nslices = mask.shape[2]
    #         for j in range(nslices):
    #             slice = mask[:, :, nslices - 1 - j]
    #             fileName = os.path.join(keyFolder, '{:03d}.png'.format(j+1))
    #             plt.imsave(fileName, slice)
    #         print('{} {}'.format(patientName, key))

    # show DVH at PTV uncropped
    target = 'PTV uncropped'
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')

        RTSTRUCTFile = os.path.join(MRdoseFolder, 'RTst0.dcm')
        MRFolder = os.path.join(MRdoseFolder, 'MR')
        rtstruct = RTStructBuilder.create_from(MRFolder, RTSTRUCTFile)
        ROInames = rtstruct.get_roi_names()
        # print('PTV uncropped in {}: {}'.format(patientName, target in ROInames))
        assert target in ROInames

        doseFile = os.path.join(MRdoseFolder, 'RTdose.dcm')
        dose = pydicom.dcmread(doseFile).pixel_array
        dose = np.flip(dose, axis=0)

        targetMask = rtstruct.get_roi_mask_by_name(target)
        targetMask = np.transpose(np.flip(targetMask, axis=2), (2, 0, 1))

        targetDose = dose[targetMask]
        targetDose = np.sort(targetDose)
        size = targetDose.size
        xAxis = np.zeros(size+1, dtype=targetDose.dtype)
        xAxis[1:] = targetDose
        yAxis = np.zeros(size+1, dtype=np.float32)
        yAxis[0] = 1
        yAxis[1:] = 1.0 - np.arange(size) / size

        plt.plot(xAxis, yAxis)
        plt.show()
        break


def visRawMR():
    """
    It seems that the MR images for patient 2 and patient 3 are the same.
    Here we examine it in the raw data

    the results show that, though the filenames of the dicom files of 
    patient 2 and patient 3 are different, they are the same
    """
    rawDataPath = os.path.join(globalFolder, 'rawData')
    visFolder = os.path.join(globalFolder, 'visualize')
    idxs = [1, 2]
    scale = 255
    for idx in idxs:
        patientID = patients[idx]
        rawMRFolder = os.path.join(rawDataPath, patientID, '*MRdose', '*MR*')
        rawMRFolder = folderEvaluate(rawMRFolder)
        files = dicomSort(rawMRFolder, 'MR')

        patientName = 'patient{}'.format(idx + 1)
        imageFolder = os.path.join(visFolder, patientName, 'rawMR')
        if not os.path.isdir(imageFolder):
            os.mkdir(imageFolder)
        for j, file in enumerate(files):
            inputPath = os.path.join(rawMRFolder, file)
            outputPath = os.path.join(imageFolder, '{:03d}.png'.format(j+1))
            pixel_array = pydicom.dcmread(inputPath).pixel_array
            roof = np.max(pixel_array)
            pixel_array = np.uint8(pixel_array / roof * scale)
            cv2.imwrite(outputPath, pixel_array)
        print('{} {}'.format(patientName, patientID))


def PTVexamination():
    """
    through visual exam, we found that seems PTV uncropped is divided 
    into PTV_HIGH and PTV_LOW. PTV_CROPPED is smaller than PTV uncropped, 
    but I haven't figured out other relationship
    """

    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'MRdose')
        MRFolder = os.path.join(MRdoseFolder, 'MR')
        RTSTRUCTFile = os.path.join(MRdoseFolder, 'RTst0.dcm')
        rtstruct = RTStructBuilder.create_from(MRFolder, RTSTRUCTFile)
        ROInames = rtstruct.get_roi_names()

        PTVuncroppedName = [a for a in ROInames if 'PTV' in a and 'uncrop' in a]
        PTVhighName = [a for a in ROInames if 'PTV' in a and 'HIGH']


def anonymizeCT():
    rawDataPath = os.path.join(globalFolder, 'rawData')
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    for i in range(num_patients):
        patientID = patients[i]
        rawPatFolder = os.path.join(rawDataPath, patientID)
        rawCTct = os.path.join(rawPatFolder, '*CT', '*_CT_*')
        rawCTct = folderEvaluate(rawCTct)
        rawCTRT = os.path.join(rawPatFolder, '*CT', '*_RTst_*')
        rawCTRT = folderEvaluate(rawCTRT)

        patientName = 'patient{}'.format(i+1)
        newPatFolder = os.path.join(anonymousDataPath, patientName)
        newCTFolder = os.path.join(newPatFolder, 'CT')
        newDicomFolder = os.path.join(newCTFolder, 'dicom')
        if not os.path.isdir(newDicomFolder):
            os.makedirs(newDicomFolder)

        # copy CT dicoms
        files = os.listdir(rawCTct)
        for file in files:
            fullFile = os.path.join(rawCTct, file)
            dicomData = pydicom.dcmread(fullFile)
            wash(dicomData)
            outputPath = os.path.join(newDicomFolder, file)
            dicomData.save_as(outputPath)
        
        # copy RTstruct
        rawRTfile = os.listdir(rawCTRT)[0]
        rawRTfile = os.path.join(rawCTRT, rawRTfile)
        dicomData = pydicom.dcmread(rawRTfile)
        wash(dicomData)
        outputPath = os.path.join(newCTFolder, 'RTst.dcm')
        dicomData.save_as(outputPath)
        print('{} {}'.format(patientName, patientID))


def visualizeCT():
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    scale = 255
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(anonymousDataPath, patientName)
        dicomFolder = os.path.join(patFolder, 'CT', 'dicom')
        files = dicomSort(dicomFolder, 'CT')

        visPatFolder = os.path.join(visFolder, patientName)
        visPatCTFolder = os.path.join(visPatFolder, 'CT')
        if not os.path.isdir(visPatCTFolder):
            os.mkdir(visPatCTFolder)
        
        for j, file in enumerate(files):
            inputFile = os.path.join(dicomFolder, file)
            pixel_array = pydicom.dcmread(inputFile).pixel_array
            roof = np.max(pixel_array)
            pixel_array = np.uint8(pixel_array / roof * scale)
            outputFile = os.path.join(visPatCTFolder, '{:03d}.png'.format(j+1))
            cv2.imwrite(outputFile, pixel_array)
        print(patientName)


def visualizeRTstructCT():
    """
    This function visualizes the RT structures of different patients,
    each patient has three files that are of RTSTRUCT modality, they
    are: RTst0.dcm, RTst1.dcm, and SEG.dcm

    Result:
    we could not open SEG.dcm using RTStructBuilder package
    RTst0 and RTst1 are the same (at least visually the same)
    the shape of RTSTRUCT mask is (hight, width, slice)
    the slice dimension is reversed from InstanceNumber
    """
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    visFolder = os.path.join(globalFolder, 'visualize')
    
    for i in range(1, num_patients):
        patientName = 'patient{}'.format(i+1)
        MRdoseFolder = os.path.join(anonymousDataPath, patientName, 'CT')
        MRFolder = os.path.join(MRdoseFolder, 'dicom')
        file = 'RTst.dcm'
        fullFile = os.path.join(MRdoseFolder, file)
        try:
            rtstruct = RTStructBuilder.create_from(
                dicom_series_path=MRFolder, rt_struct_path=fullFile)
        except:
            print('{} error'.format(patientName))
            continue
        ROInames = rtstruct.get_roi_names()
        print('{} {}\n{}\n\n'.format(patientName, file, ROInames))
        masks = {}
        for name in ROInames:
            try:
                masks[name] = rtstruct.get_roi_mask_by_name(name)
            except:
                print('fail to extract the mask for {}'.format(name))
        
        # create folder
        visRTstructFolder = os.path.join(visFolder, patientName, 'CTrt')
        if not os.path.isdir(visRTstructFolder):
            os.makedirs(visRTstructFolder)
        
        # draw contours
        MRfiles = dicomSort(MRFolder, 'CT')
        for j, file in enumerate(MRfiles):
            MRfile = os.path.join(MRFolder, file)
            result = draw_mask(MRfile, masks, len(MRfiles) - 1 - j)
            outputPath = os.path.join(visRTstructFolder, '{:03d}.png'.format(j+1))
            cv2.imwrite(outputPath, result)
            print('{} {} {}'.format(patientName, file, j+1))


def CTrtExamine():
    """
    The function visualizeRTstructCT has some errors, specifically, it 
    reports that the frame of reference UID of the RTstruct file is not 
    consistent with the CT dicom files
    """
    anonymousDataPath = os.path.join(globalFolder, 'anonymousData')
    for i in range(num_patients):
        patientName = 'patient{}'.format(i+1)
        patCTFolder = os.path.join(anonymousDataPath, patientName, 'CT')
        RTstFile = os.path.join(patCTFolder, 'RTst.dcm')
        dicomFolder = os.path.join(patCTFolder, 'dicom')
        CTfiles = dicomSort(dicomFolder, 'CT')
        CTfile = os.path.join(dicomFolder, CTfiles[0])
        CTdata = pydicom.dcmread(CTfile)
        RTdata = pydicom.dcmread(RTstFile)
        continue



if __name__ == '__main__':
    # examineMR()
    # anonymizeMR()
    # visualizeMR()
    # visualizeRTPLAN()
    # visualizeRTdose()
    # visualizeRTstruct()
    # test_cv2_drawContour()
    # visualizeDVH()
    # DVHdraft()
    # printAnatomy()
    # visRawMR()
    # anonymizeCT()
    # visualizeCT()
    # visualizeRTstructCT()
    CTrtExamine()