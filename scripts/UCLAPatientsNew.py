import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import pydicom
import cv2
from rt_utils import RTStructBuilder

"""
Last time we processed some patients. However, we noticed that the RTstruct 
of CT is rigidly registered from that of corresponding MR. And the geometry
of CT and MR differred a lot. So this time, we selected another set of patients.
"""

globalFolder = '/data/datasets/UCLAPatients'
rawFolder = '/data/datasets/UCLAPatients/rawDataNew'
anonymousFolder = os.path.join(globalFolder, 'anonymousDataNew')
visuallize = os.path.join(globalFolder, 'visualizeNew')
patients = os.listdir(rawFolder)
patientTree = []
numPatients = 8

# list all attributes to wipe out in CT
AttrWipeCT = ['AccessionNumber', 'AcquisitionDate', 'AcquisitionNumber', 'AcquisitionTime', 'ContentDate',
    'ContentTime', 'InstanceCreationDate', 'InstanceCreationTime', 'InstitutionName', 'Manufacturer', 
    'ManufacturerModelName', 'OperatorsName', 'PatientAge', 'PatientBirthDate', 'PatientName', 'PatientSex', 
    'PatientSize', 'PatientWeight', 'PerformingPhysicianName', 'ReferringPhysicianName', 'SeriesDate', 
    'StationName', 'StudyDate']

AttrWipeCTrt = ['AccessionNumber', 'ApprovalStatus', 'InstanceCreationDate', 'InstanceCreationTime', 'InstitutionName',
    'Manufacturer', 'ManufacturerModelName', 'OperatorsName', 'PatientAge', 'PatientBirthDate', 'PatientID', 'PatientName',
    'PatientSex', 'PatientSize', 'PatientWeight', 'PerformingPhysicianName', 'ReferringPhysicianName', 'SeriesDate', 
    'SeriesDescription', 'SoftwareVersions', 'StationName', 'StructureSetDate', 'StructureSetName', 'StructureSetTime', 
    'StudyDate', 'parent']

AttrWipeMR = ['AccessionNumber', 'AcquisitionDate', 'AcquisitionNumber', 'AcquisitionTime', 'ContentDate', 'ContentTime',
    'InstanceCreationDate', 'InstanceCreationTime', 'InstitutionName', 'Laterality', 'MRAcquisitionType', 'Manufacturer', 
    'ManufacturerModelName', 'OperatorsName', 'PatientAge', 'PatientBirthDate', 'PatientID', 'PatientName', 'PatientSex', 
    'ReferringPhysicianName', 'SeriesDate', 'StationName', 'StudyDate', 'StudyDescription', 'StudyID', 'StudyTime']

AttrWipeMRrt = ['AccessionNumber', 'ApprovalStatus', 'InstanceCreationDate', 'InstanceCreationTime', 'InstitutionName', 
    'Manufacturer', 'ManufacturerModelName', 'OperatorsName', 'PatientAge', 'PatientBirthDate', 'PatientID', 'PatientName', 
    'PatientSex', 'PerformingPhysicianName', 'ReferringPhysicianName', 'SeriesDate', 'SeriesDescription', 'SeriesNumber', 
    'SeriesTime', 'SoftwareVersions', 'StationName', 'StructureSetDate', 'StructureSetName', 'StructureSetTime', 'StudyDate',
    'StudyDescription', 'StudyID', 'StudyTime', 'parent']

AttrWipeDose = ['AccessionNumber', 'AcquisitionDate', 'AcquisitionTime', 'ContentDate', 'ContentTime', 'InstanceCreationDate', 
    'InstanceCreationTime', 'InstitutionName', 'Laterality', 'Manufacturer', 'ManufacturerModelName', 'OperatorsName', 
    'PatientAge', 'PatientBirthDate', 'PatientID', 'PatientName', 'PatientSex', 'ReferringPhysicianName', 'SeriesDate', 
    'SeriesDescription', 'SeriesTime', 'StudyDate', 'StudyDescription', 'StudyTime', 'parent']


def folderEvaluate(template):
    results = glob.glob(template)
    assert len(results) == 1, \
        "{} matches 0 or more than 1 file/folders.".format(template)
    return results[0]


def findNumber(chunks, modality='CT'):
    for i, chunk in enumerate(chunks):
        if chunk == modality:
            return chunks[i+2]


def FullPatientTree():
    """
    This function finds the corresponding files of individual patients.
    """
    global patientTree
    for patient in patients:
        patFolder = os.path.join(rawFolder, patient)
        subFolder = os.path.join(patFolder, '*Studies')
        subFolder = folderEvaluate(subFolder)

        # find CT folder
        CTFolder = folderEvaluate(os.path.join(subFolder, '*_CT_*'))
        chunks = CTFolder.split('/')[-1]
        chunks = chunks.split('_')
        number = findNumber(chunks)
        
        # find CT RTst folder
        CTRTstFolder = folderEvaluate(os.path.join(subFolder, '*_RTst_*_{}_*00000'.format(number)))
        
        # find MR folder
        MRFolder = folderEvaluate(os.path.join(subFolder, '*_MR_*'))
        chunks = MRFolder.split('/')[-1]
        chunks = chunks.split('_')
        number = findNumber(chunks, modality='MR')

        # find MR RTst folder
        MRRTstFolder = folderEvaluate(os.path.join(subFolder, '*_RTst_*_{}_*00000'.format(number)))

        # find MR dose folder
        RTDOSEFolder = folderEvaluate(os.path.join(subFolder, '*_RTDOSE_*_{}_*'.format(number)))
        
        patientTree.append({'CT': CTFolder, 'CTrt': CTRTstFolder, 'MR': MRFolder, 
            'MRrt': MRRTstFolder, 'dose': RTDOSEFolder})

    # # take a look
    # attrs = patientTree[0].keys()
    # for pat in patientTree:
    #     for key in attrs:
    #         print('{} {}'.format(key, pat[key]))
    #     print('\n')


FullPatientTree()


def anonymize():
    if not os.path.isdir(anonymousFolder):
        os.mkdir(anonymousFolder)
    modality = 'dose'

    # anonymize CT files
    if modality == 'CT':
        # attrs = []
        for i in range(len(patientTree)):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            if not os.path.isdir(patFolder):
                os.mkdir(patFolder)

            patInfo = patientTree[i]
            CTFolderOld = patInfo['CT']
            CTfiles = dicomSort(CTFolderOld, modality)

            # # take a look
            # CTfile = os.path.join(CTFolderOld, CTfiles[5])
            # CTdata = pydicom.dcmread(CTfile)
            # attr = CTdata.__dir__()

            # for a in AttrWipeCT:
            #     if a not in attr:
            #         print(a)
            # break

            CTFolderNew = os.path.join(patFolder, 'CT')
            if not os.path.isdir(CTFolderNew):
                os.mkdir(CTFolderNew)
            for j, CTfile in enumerate(CTfiles):
                oldFile = os.path.join(CTFolderOld, CTfile)
                newFile = os.path.join(CTFolderNew, '{:03d}.dcm'.format(j+1))
                CTdata = pydicom.dcmread(oldFile)
                wipeCT(CTdata)
                CTdata.save_as(newFile)
            print(patientName)

            # attrs.append(attr)
        # checkAttrs(attrs)
    
    elif modality == 'CTrt':
        # attrs = []
        for i in range(len(patientTree)):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)

            patInfo = patientTree[i]
            CTrtFolder = patInfo[modality]
            file = os.listdir(CTrtFolder)[0]
            file = os.path.join(CTrtFolder, file)
            data = pydicom.dcmread(file)
            # attr = data.__dir__()
            # for a in AttrWipeCTrt:
            #     if a not in attr:
            #         print(a)
            # break
            # attrs.append(attr)
            wipeCTrt(data)

            newPath = os.path.join(patFolder, 'CTrt.dcm')
            data.save_as(newPath)
            print(patientName)
        
        # # take a look
        # for attr in attrs:
        #     print(attr, '\n')
    
    elif modality == 'MR':
        # attrs = []
        for i in range(len(patientTree)):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            newMRFolder = os.path.join(patFolder, 'MR')
            if not os.path.isdir(newMRFolder):
                os.mkdir(newMRFolder)

            patInfo = patientTree[i]
            MRFolder = patInfo[modality]
            files = dicomSort(MRFolder, modality)
            
            # # have a look
            # file = os.path.join(MRFolder, files[0])
            # data = pydicom.dcmread(file)
            # attr = data.__dir__()
            # # attrs.append(attr)
            # for a in AttrWipeMR:
            #     if a not in attr:
            #         print(a)

            for j, file in enumerate(files):
                file_ = os.path.join(MRFolder, file)
                data = pydicom.dcmread(file_)
                wipeMR(data)
                newPath = os.path.join(newMRFolder, '{:03d}.dcm'.format(j+1))
                data.save_as(newPath)
            print(patientName)

        # for attr in attrs:
        #     print(attr, '\n')
    
    elif modality == 'MRrt':
        # attrs = []
        for i in range(len(patientTree)):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            newFile = os.path.join(patFolder, 'MRrt.dcm')

            patInfo = patientTree[i]
            MRrtFolder = patInfo[modality]
            rtFile = os.listdir(MRrtFolder)[0]
            rtFile = os.path.join(MRrtFolder, rtFile)
            data = pydicom.dcmread(rtFile)

            # attr = data.__dir__()
            # # attrs.append(attr)
            # for a in AttrWipeMRrt:
            #     if a not in attr:
            #         print(a)
            # break

            wipeMRrt(data)
            data.save_as(newFile)
            print(patientName)
        
        # checkAttrs(attrs)
    
    elif modality == 'dose':
        for i in range(len(patientTree)):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            newFile = os.path.join(patFolder, 'MRdose.dcm')

            patInfo = patientTree[i]
            doseFolder = patInfo[modality]
            file = os.listdir(doseFolder)[0]
            file = os.path.join(doseFolder, file)
            data = pydicom.dcmread(file)

            wipeMRdose(data)
            data.save_as(newFile)
            print(patientName)


def wipeCT(CTdata):
    for attr in AttrWipeCT:
        setattr(CTdata, attr, '')


def wipeCTrt(CTrtdata):
    for attr in AttrWipeCTrt:
        setattr(CTrtdata, attr, '')


def wipeMR(MRdata):
    for attr in AttrWipeMR:
        setattr(MRdata, attr, '')

def wipeMRrt(data):
    for attr in AttrWipeMRrt:
        setattr(data, attr, '')


def wipeMRdose(data):
    for attr in AttrWipeDose:
        setattr(data, attr, '')


def checkAttrs(attrs):
    """
    Some sttributes are longer, others are shorter
    """
    template = attrs[0]
    for i in range(1, len(attrs)):
        comp = attrs[i]
        for attr in template:
            if attr not in comp:
                print(attr)
        print('\n')


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


def visCTMR():
    """
    Firstly, convert the DICOM images of CT and MR to PNG images for 
    visualization
    """



if __name__ == '__main__':
    anonymize()