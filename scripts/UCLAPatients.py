import os
import glob
import pydicom

# this script runs on shenggpu4. Here we specify the path to the data
globalFolder = '/data/datasets/UCLAPatients'
patients = ['0530793', '0678774', '0999272', '0999272', '0999272', '0999272', '0999272']

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
        'OperatorName', 'PatientBirthDate', 'PatientID', 'PatientName', \
        'PatientSex', 'PerformingPhysicianName']
    for attribute in attributes:
        if hasattr(dicomData, attribute):
            setattr(dicomData, attribute, '')
    if hasattr(dicomData, 'file_meta'):
        if hasattr(dicomData.file_meta, 'ImplementationVersionName'):
            dicomData.file_meta.ImplementationVersionName = ''
        if hasattr(dicomData.file_meta, 'SourceApplicationEntityTitle'):
            dicomData.file_meta.SourceApplicationEntityTitle = ''



if __name__ == '__main__':
    # examineMR()
    anonymizeMR()