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
visFolder = os.path.join(globalFolder, 'visNew')
patients = os.listdir(rawFolder)
patients.sort()
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


anatomies = [
    ['Skin', 'GTV', 'PTV 50 uncropped', 'O_Cord', 'O_Esgs', 'O_Hrt', 'O_Livr', 'O_Stmc', 'O_Kdny_Rt', 
    'O_Kdny_Lt', 'O_Bwel_Sm', 'O_Bwel_Lg', 'O_Duodenum', 'PTV_CROP', 'PTV_LOW', 'PTV_HIGH'],

    ['Skin', 'O_Cord', 'O_Esgs', 'O_Hrt', 'O_Livr', 'O_Stmc', 'O_Kdny_Rt', 'O_Kdny_Lt', 'O_Bwel_Sm', 
    'O_Bwel_Lg', 'GTV ACR', 'PTV ACR uncroppe', 'PTV_CROPPED', 'PTV_LOW', 'PTV_HIGH'],

    ['Skin', 'O_Cord', 'O_Esgs', 'O_Hrt', 'O_Livr', 'O_Stmc', 'O_Kdny_Rt', 'O_Kdny_Lt', 'O_Bwel', 
    'O_Duod', 'GTV', 'PTV uncropped', 'PTV_CROPPED', 'PTV_LOW', 'PTV_HIGH'],

    ['Skin', 'O_Cord', 'O_Hrt', 'O_Livr', 'O_Stmc', 'O_Kdny_Rt', 'O_Bwel', 'O_Vssl', 'GTV', 
    'PTV uncropped', 'O_DUOD', 'PTV_CROPPED', 'PTV_HIGH', 'PTV_LOW'],

    ['Skin', 'O_Cord', 'O_Lung_Lt', 'O_Lung_Rt', 'O_Hrt', 'O_Livr', 'G_Site', 'O_Lung_Tt', 
    'O_Kdny_Rt', 'O_Kdny_Lt', 'O_Vssl', 'O_Bwl', 'O_Duod', 'O_Spln', 'PTV 3mm uncroppe', 
    'O_Stmch', 'PTV_HIGH', 'PTV_LOW', 'PTV_CROPPED'],

    ['O_Bldr', 'O_Femr_Lt', 'O_Femr_Rt', 'O_Rctm', 'Skin', 'O_Cord', 'O_Spln', 'O_Livr', 'O_Stmc', 
    'O_Kdny_Rt', 'O_Kdny_Lt', 'O_Bwel', 'O_Duod', 'O_Vssl', 'GTV', 'PTV 50 uncropped', 'PTV_CROPPED', 
    'PTV_HIGH', 'PTV_LOW'],

    ['Skin', 'GTV', 'PTV 50_3mm', 'PTV_High', 'O_Duod', 'O_Stmc', 'O_Bwel_Sm', 'O_Bwel_Lg', 'O_Cord', 
    'O_Esgs', 'O_Vessle', 'O_Kdny_Rt', 'O_Kdny_Lt', 'O_Livr', 'PTV_crop', 'O_Nrlms'],

    ['Skin', 'GTV', 'PTV_3MM', 'PTV_crop', 'PTV_High', 'O_Duod', 'O_Stmc', 'O_Bwel_Lg', 'O_Esgs', 
    'O_Vessle', 'O_Cord', 'O_Kdny_Rt', 'O_Kdny_Lt', 'O_Livr', 'O_Bwel_Sm', 'O_Nrlms', 'O_Panc']
]


colors = [(231, 203, 107), (152, 153, 139), (97, 111, 175), (60, 42, 162), (75, 172, 12), 
    (35, 69, 171), (86, 201, 174), (16, 103, 52), (204, 19, 136), (121, 104, 80), 
    (96, 32, 200), (57, 131, 240), (17, 97, 12), (255, 218, 243), (245, 242, 167), 
    (94, 100, 206), (238, 227, 161), (202, 254, 156), (76, 41, 217), (158, 239, 248)]


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


def allAnonymize():
    modalities = ['CT', 'CTrt', 'MR', 'MRrt', 'dose']
    for modality in modalities:
        anonymize(modality)


def anonymize(modality):
    if not os.path.isdir(anonymousFolder):
        os.mkdir(anonymousFolder)
    # modality = 'CT'
    # modality = 'CTrt'
    # modality = 'MR'
    # modality = 'MRrt'
    # modality = 'dose'

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
    if not os.path.isdir(visFolder):
        os.mkdir(visFolder)
    scale = 255
    modality = 'CT'

    # save MR images
    if modality == 'MR':
        for i in range(numPatients):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            MRFolder = os.path.join(patFolder, 'MR')
            MRoutput = os.path.join(visFolder, patientName, 'MR')
            if not os.path.isdir(MRoutput):
                os.makedirs(MRoutput)
            files = os.listdir(MRFolder)
            files.sort()
            for file in files:
                file_ = os.path.join(MRFolder, file)
                pixel_array = pydicom.dcmread(file_).pixel_array
                roof = np.max(pixel_array)
                pixel_array = np.uint8(pixel_array / roof * scale)
                outFile = file.split('.')[0]
                outFile = os.path.join(MRoutput, outFile + '.png')
                plt.imsave(outFile, pixel_array, cmap='gray')
            print(patientName)
    elif modality == 'CT':
        for i in range(numPatients):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            CTFolder = os.path.join(patFolder, 'CT')
            CToutput = os.path.join(visFolder, patientName, 'CT')
            if not os.path.isdir(CToutput):
                os.makedirs(CToutput)
            files = os.listdir(CTFolder)
            files.sort()
            for file in files:
                file_ = os.path.join(CTFolder, file)
                pixel_array = pydicom.dcmread(file_).pixel_array
                roof = np.max(pixel_array)
                pixel_array = np.uint8(pixel_array / roof * scale)
                outFile = file.split('.')[0]
                outFile = os.path.join(CToutput, outFile + '.png')
                plt.imsave(outFile, pixel_array, cmap='gray')
            print(patientName)


def showAnatomy():
    # task = 'sanity'
    task = 'MR'
    # task = 'CT'

    # Firstly, check that the selected anatomies exist
    if task == 'sanity':
        for i in range(numPatients):
            idx = i + 1
            patientName = 'patient{}'.format(idx)
            patFolder = os.path.join(anonymousFolder, patientName)
            MRrtPath = os.path.join(patFolder, 'MRrt.dcm')
            MRFolder = os.path.join(patFolder, 'MR')
            MRstruct = RTStructBuilder.create_from(dicom_series_path=MRFolder, 
                rt_struct_path=MRrtPath)
            MRROInames = MRstruct.get_roi_names()

            CTrtPath = os.path.join(patFolder, 'CTrt.dcm')
            CTFolder = os.path.join(patFolder, 'CT')
            CTstruct = RTStructBuilder.create_from(dicom_series_path=CTFolder, 
                rt_struct_path=CTrtPath)
            CTROInames = CTstruct.get_roi_names()

            ana = anatomies[i]
            flag = True
            for name in ana:
                if name not in MRROInames:
                    flag = False
                    break
            print('{} MR {}'.format(patientName, 'correct' if flag else 'incorrect'))
            flag = True
            for name in ana:
                if name not in CTROInames:
                    flag = False
                    break
            print('{} CT {}\n'.format(patientName, 'correct' if flag else 'incorrect'))

            # print('{}\n{}\n{}\n'.format(patientName, MRROInames, CTROInames))
    
    # Secondly, draw contours
    else:
        for i in range(numPatients):
            # i = 4
            patientName = 'patient{}'.format(i+1)
            patFolder = os.path.join(anonymousFolder, patientName)
            dicomFolder = os.path.join(patFolder, task)
            RTfile = os.path.join(patFolder, task + 'rt.dcm')
            RT = RTStructBuilder.create_from(dicom_series_path=dicomFolder, rt_struct_path=RTfile)
            anatomyList = anatomies[i]
            result = drawContours(dicomFolder, RT, anatomyList)
            outFolder = os.path.join(visFolder, patientName, task+'contour')
            writeResult(outFolder, result)
            print(patientName)


def writeResult(outFolder, result):
    if not os.path.isdir(outFolder):
        os.makedirs(outFolder)
    nSlices = result.shape[3]
    for i in range(nSlices):
        fileName = '{:03d}.png'.format(i+1)
        fileName = os.path.join(outFolder, fileName)
        plt.imsave(fileName, result[:, :, :, i])
        

def drawContours(dicomFolder, RT, anatomyList):
    """
    This function draws contour on the dicom images
    """
    # firstly, get the shape of the images
    exampleMask = RT.get_roi_mask_by_name(anatomyList[0])
    shape = exampleMask.shape
    height, width, slices = shape
    channels = 3
    result = np.zeros((height, width, channels, slices), dtype=np.uint8)

    # # take a look
    # print(np.max(exampleMask))
    # plt.imshow(exampleMask[:, :, 40])
    # plt.show()
    # exit()

    # get slices first
    scale = 255
    for i in range(slices):
        slicePath = os.path.join(dicomFolder, '{:03d}.dcm'.format(i+1))
        pixel_array = pydicom.dcmread(slicePath).pixel_array
        roof = np.max(pixel_array)
        pixel_array = np.uint8(pixel_array / roof * scale)
        pixel_array = np.stack((pixel_array, ) * channels, axis=2)
        result[:, :, :, i] = pixel_array
    
    image = np.ones((height, width, channels), dtype=np.uint8) * scale
    for j, anatomy in enumerate(anatomyList):
        mask = RT.get_roi_mask_by_name(anatomy)
        mask = np.flip(mask, axis=2)
        for i in range(slices):
            maskSlice = mask[:, :, i]
            maskSlice = np.uint8(maskSlice * scale)
            contours, _ = cv2.findContours(maskSlice, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            image[:] = result[:, :, :, i]
            result[:, :, :, i] = cv2.drawContours(image, contours, -1, colors[j], 1)
    return result


def colorGen():
    """
    this is to generate 20 different colors
    """
    colors = []
    Ncolors = 20
    norm = 256
    for i in range(Ncolors):
        b = int(np.random.rand() * norm)
        g = int(np.random.rand() * norm)
        r = int(np.random.rand() * norm)
        color = (b, g, r)
        colors.append(color)
    
    for anatomy in anatomies:
        print(len(anatomy))


def visDose():
    """
    This function is to visualize the dose of MR images
    """
    roof = 65535
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(anonymousFolder, patientName)
        # MRFolder = os.path.join(patFolder, 'MR')
        doseFile = os.path.join(patFolder, 'MRdose.dcm')
        dose = pydicom.dcmread(doseFile)
        dose = dose.pixel_array
        dose = np.flip(dose, axis=0)

        visFolderDose = os.path.join(visFolder, patientName, 'dose')
        if not os.path.isdir(visFolderDose):
            os.mkdir(visFolderDose)
        for j in range(dose.shape[0]):
            outFile = '{:03d}.png'.format(j+1)
            outFile = os.path.join(visFolderDose, outFile)
            plt.imsave(outFile, dose[j, :, :], vmin=0, vmax=roof)
        print(patientName)


def visDVH():
    """
    This function draws DVH
    """
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(anonymousFolder, patientName)
        MRFolder = os.path.join(patFolder, 'MR')
        MRrtFile = os.path.join(patFolder, 'MRrt.dcm')
        MRdoseFile = os.path.join(patFolder, 'MRdose.dcm')

        # get mask
        MRrt = RTStructBuilder.create_from(dicom_series_path=MRFolder, rt_struct_path=MRrtFile)
        anatomy = anatomies[i]
        anatomy.remove('Skin')
        dose = pydicom.dcmread(MRdoseFile).pixel_array
        dose = np.flip(dose, axis=0)  # (slice, height, width)

        for name in anatomy:
            mask = MRrt.get_roi_mask_by_name(name)
            mask = np.flip(mask, axis=2)  # (height, width, slice)
            mask = np.transpose(mask, (2, 0, 1))
            maskedDose = dose[mask]
            maskedDose = np.sort(maskedDose)
            size = maskedDose.size
            xAxis = np.zeros(size+1, dtype=maskedDose.dtype)
            xAxis[1:] = maskedDose
            yAxis = np.zeros(size+1, dtype=np.float32)
            yAxis[0] = 1
            yAxis[1:] = 1 - np.arange(size) / size
            plt.plot(xAxis, yAxis)
        plt.legend(anatomy)
        plt.title('{} DVH'.format(patientName))
        plt.xlabel('dose (a.u.)')
        plt.ylabel('relative volume')
        outputFile = os.path.join(visFolder, patientName, 'DVH.png')
        plt.savefig(outputFile)
        plt.clf()
        print(patientName)


def checkCTrt():
    """
    It seems that the RTstruct for CT is empty. This is to check the RTstructures.
    It turns out that CT anatomy is empty
    """
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        patFolder = os.path.join(anonymousFolder, patientName)
        CTFolder = os.path.join(patFolder, 'CT')
        CTrtFile = os.path.join(patFolder, 'CTrt.dcm')
        # anatomy = anatomies[i]
        CTrt = RTStructBuilder.create_from(dicom_series_path=CTFolder, rt_struct_path=CTrtFile)

        # skin is in anatomy
        skinMask = CTrt.get_roi_mask_by_name('Skin')
        print(np.max(skinMask))


def moveDVH():
    """
    This function moves DVH images to a single file for forward
    """
    destFolder = os.path.join(visFolder, 'DVH')
    if not os.path.isdir(destFolder):
        os.mkdir(destFolder)
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        source = os.path.join(visFolder, patientName, 'DVH.png')
        dest = os.path.join(destFolder, patientName + '.png')
        command = 'cp {} {}'.format(source, dest)
        os.system(command)
        


if __name__ == '__main__':
    # allAnonymize()
    # anonymize()
    # visCTMR()
    # showAnatomy()
    # colorGen()
    # visDose()
    # visDVH()
    # checkCTrt()
    moveDVH()