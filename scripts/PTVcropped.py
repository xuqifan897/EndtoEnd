import os
import glob
import sys
from PyQt5 import QtWidgets
QtWidgets.QApplication(sys.argv)
import numpy as np
import json
import pydicom
import matplotlib.pyplot as plt
import pydicom
from rt_utils import RTStructBuilder
from scipy.io import loadmat
from drawDVH import ReSize, DVHgroup

globalFolder = '/data/datasets/UCLAPatients'
expFolder = os.path.join(globalFolder, 'experiment')
visFolder = os.path.join(globalFolder, 'visNew')
anonymousDataFolder = os.path.join(globalFolder, 'anonymousDataNew')

def reviseStructures():
    """
    Previously, we tried to use the PTV uncropped. However, 
    for some of the patients, the duodenum overlaps with the PTV,
    which makes it impossible to both impose high dose to PTV and
    sparse the OAR. So here we choose to use PTV cropped instead
    """

    PTVs = ['PTV_CROP', 'PTV_CROPPED', 'PTV_CROPPED', 'PTV_CROPPED', 
        'PTV_CROPPED', 'PTV_CROPPED', 'PTV_crop', 'PTV_crop']
    numPatients = 8

    # firstly, we print the structures
    if False:
        for i in range(numPatients):
            patientName = 'patient{}'.format(i+1)
            expPatFolder = os.path.join(expFolder, patientName)
            dataPatFolder = os.path.join(anonymousDataFolder, patientName)
            MRFolder = os.path.join(dataPatFolder, 'MR')
            MRrtFile = os.path.join(dataPatFolder, 'MRrt.dcm')

            rtstruct = RTStructBuilder.create_from(MRFolder, MRrtFile)
            names = rtstruct.get_roi_names()
            names_selected = [a for a in names if 'PTV' in a]
            assert PTVs[i] in names_selected
    
    # secondly, we revise the structures
    if True:
        for i in range(numPatients):
            patientName = 'patient{}'.format(i+1)
            expPatFolder = os.path.join(expFolder, patientName)
            optFolderNew = os.path.join(expPatFolder, 'optimizePTVcropped')
            if not os.path.isdir(optFolderNew):
                os.mkdir(optFolderNew)
            structuresFileOrg = os.path.join(expPatFolder, 'structures.json')
            with open(structuresFileOrg) as f:
                structuresOrg = json.load(f)
            structuresOrg['ptv'] = PTVs[i]
            # print(structuresOrg, '\n')
            structuresFileNew = os.path.join(optFolderNew, 'structures.json')
            with open(structuresFileNew, 'w') as f:
                json.dump(structuresOrg, f)
            print(patientName)
    
    # thirdly, we create a shorter beamlist for dose calculation. 
    # The main purpose of dose calculation here is to generate the 
    # new mask file.
    if True:
        for i in range(numPatients):
            patientName = 'patient{}'.format(i+1)
            expPatFolder = os.path.join(expFolder, patientName)
            beamlistOrg = os.path.join(expPatFolder, 'beamlist.txt')
            with open(beamlistOrg, 'r') as f:
                line = f.readline()  # we only read one line for simplicity
            optFolderNew = os.path.join(expPatFolder, 'optimizePTVcropped')
            beamlistNew = os.path.join(optFolderNew, 'beamlist.txt')
            with open(beamlistNew, 'w') as f:
                f.writelines(line)
            print(patientName)



def createCTanatomy():
    """
    Previously, there were something wrong with the PTV_CROP for some patients
    Here we propose to use a different way to redo this step using something new
    """
    numPatients = 8
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRFolder = os.path.join(dataPatFolder, 'MR')
        CTAlignFolder = os.path.join(dataPatFolder, 'CTAlign')
        MRrtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        
        # firstly, compare the differences in metadata 
        # of dicom files of CT and MR
        if False:
            exampleCTFile = os.listdir(CTAlignFolder)
            exampleCTFile.sort()
            exampleCTFile = exampleCTFile[0]
            exampleCTFile = os.path.join(CTAlignFolder, exampleCTFile)
            CTdata = pydicom.dcmread(exampleCTFile)
            
            exampleMRFile = os.listdir(MRFolder)
            exampleMRFile.sort()
            exampleMRFile = exampleMRFile[0]
            exampleMRFile = os.path.join(MRFolder, exampleMRFile)
            MRdata = pydicom.dcmread(exampleMRFile)

            CTSliceSize = CTdata.Rows * CTdata.PixelSpacing[0]
            MRSliceSize = MRdata.Rows * MRdata.PixelSpacing[0]
            print(CTSliceSize, MRSliceSize, '\n')
        
        # then, we copy the CT images
        if False:
            newCTFolder = os.path.join(dataPatFolder, 'CT_PTVcrop')
            if not os.path.isdir(newCTFolder):
                os.mkdir(newCTFolder)
            oldCTFolder = os.path.join(dataPatFolder, 'CTAlignResize')
            files = os.listdir(oldCTFolder)
            for file in files:
                if file == 'CTrtFilter.dcm':
                    continue
                fileORG = os.path.join(oldCTFolder, file)
                fileNew = os.path.join(newCTFolder, file)
                command = 'cp {} {}'.format(fileORG, fileNew)
                os.system(command)
            print(patientName)
        
        # then, try to read MR structures
        if False:
            RTstruct = RTStructBuilder.create_from(MRFolder, MRrtFile)
            names = RTstruct.get_roi_names()
            names_filtered = []
            for name in names:
                if name not in names_filtered:
                    names_filtered.append(name)
            # print(names_filtered, '\n')
            for name in names:
                try:
                    mask = RTstruct.get_roi_mask_by_name(name)
                except:
                    print('{} {} error'.format(patientName, name))
        
        # then, try to extract structures using CT folder and MRrt
        # got SOPInstanceUID problem
        if False:
            RTstructData = pydicom.dcmread(MRrtFile)
            break
            
        # then, check SOPInstanceUID


def doseResizeCropped():
    """
    Firstly, resize the dose after the BOO 
    optimization to the original shape
    """
    numPatients = 8
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimizePTVcropped')

        # get trailNO
        template = os.path.join(optFolder, 'polishResult*.mat')
        files = glob.glob(template)
        files.sort()
        polishResultFile = files[-1]
        polishResultFile_ = polishResultFile.split('/')[-1]
        digits = [a for a in polishResultFile_ if a.isdigit()]
        trailNO = ''.join(digits)
        trailNO = int(trailNO)

        # get original dose
        polishData = loadmat(polishResultFile)['polishResult'][0, 0]
        polishDose = polishData['dose']
        print(polishDose.shape)

        # get original dose shape
        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        MRdose = pydicom.dcmread(MRdoseFile).pixel_array
        MRdoseShape = MRdose.shape
        targetShape = (MRdoseShape[1], MRdoseShape[2], MRdoseShape[0])
        print(targetShape)

        resizedDose = ReSize(polishDose, targetShape)
        
        targetFile = os.path.join(optFolder, 'dose{}.npy'.format(trailNO))
        np.save(targetFile, resizedDose)
        print(patientName, '\n')


def visDVH_PTVcropped():
    """
    This function draws the DVH graph of the 
    optimization results using PTV cropped.
    """
    numPatients = 8
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimizePTVcropped')
        # get trailNO
        template = os.path.join(optFolder, 'dose*.npy')
        files = glob.glob(template)
        files.sort()
        BOOdoseFile = files[-1]
        BOOdoseFile_ = BOOdoseFile.split('/')[-1]
        digits = [a for a in BOOdoseFile_ if a.isdigit()]
        trailNO = ''.join(digits)
        trailNO = int(trailNO)
        BOOdose = np.load(BOOdoseFile)

        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        clinicalDose = pydicom.dcmread(MRdoseFile).pixel_array
        clinicalDose = np.transpose(clinicalDose, (1, 2, 0))

        # load structures
        structuresFile = os.path.join(optFolder, 'structures.json')
        with open(structuresFile, 'r') as f:
            structures = json.load(f)
        PTVname = structures['ptv']
        OARnames = structures['oar']
        ROIs = [PTVname, ] + OARnames[1:]  # to exclude skin/body

        # load masks
        MRFolder = os.path.join(dataPatFolder, 'MR')
        rtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=MRFolder, rt_struct_path=rtFile)
        masks = {name: rtstruct.get_roi_mask_by_name(name) for name in ROIs}

        # dose normalization. As the clinical Dose and 
        # BOOdose are not of the same scale
        # here we normalize using D98.
        PTVmask = masks[PTVname]
        BOOPTV = BOOdose[PTVmask]
        clinicalPTV = clinicalDose[PTVmask]
        thresh = 5
        BOOPercentile = np.percentile(BOOPTV, thresh)
        clinicalPercentile = np.percentile(clinicalPTV, thresh)
        clinicalDose = clinicalDose / clinicalPercentile * BOOPercentile
        
        # draw DVH
        colorList = DVHgroup(BOOdose, masks)
        DVHgroup(clinicalDose, masks, linestyle='--', colorList=colorList)
        visPatFolder = os.path.join(visFolder, patientName)
        outFile = os.path.join(visPatFolder, 'DVH_PTVcropped{}.png'.format(trailNO))
        plt.legend(ROIs)
        plt.savefig(outFile)
        plt.clf()
        print(patientName)


def DVHgroupAAPM(dose, masks, ax, linestyle='-', colorList=None):
    if colorList is None:
        colorList = []
        for mask, array in masks.items():
            doseMasked = dose[array]
            doseMasked = np.sort(doseMasked)
            doseMasked = np.insert(doseMasked, 0, 0)
            yAxis = 1 - np.arange(len(doseMasked)) / len(doseMasked)
            yAxis = yAxis * 100
            p = ax.plot(doseMasked, yAxis, linestyle=linestyle)
            colorList.append(p[0].get_color())
        return colorList
    else:
        count = 0
        for mask, array in masks.items():
            doseMasked = dose[array]
            doseMasked = np.sort(doseMasked)
            doseMasked = np.insert(doseMasked, 0, 0)
            yAxis = 1 - np.arange(len(doseMasked)) / len(doseMasked)
            yAxis = yAxis * 100
            p = ax.plot(doseMasked, yAxis, linestyle=linestyle, color=colorList[count])
            count += 1


def prepareAAPMSuppDoc():
    """
    My phd supervisor suggested me to submit the current results to AAPM
    Here we prepare the figures for the supporting documents.
    """
    numPatients = 8
    rowSize = 4
    colSize = int(numPatients / rowSize)
    fig, axs = plt.subplots(colSize, rowSize)
    fontSize = 16
    labelFontSize = 10
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimizePTVcropped')
        # get trailNO
        template = os.path.join(optFolder, 'dose*.npy')
        files = glob.glob(template)
        files.sort()
        BOOdoseFile = files[-1]
        BOOdoseFile_ = BOOdoseFile.split('/')[-1]
        digits = [a for a in BOOdoseFile_ if a.isdigit()]
        trailNO = ''.join(digits)
        trailNO = int(trailNO)
        BOOdose = np.load(BOOdoseFile)

        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        clinicalDose = pydicom.dcmread(MRdoseFile).pixel_array
        clinicalDose = np.transpose(clinicalDose, (1, 2, 0))

        # load structures
        structuresFile = os.path.join(optFolder, 'structures.json')
        with open(structuresFile, 'r') as f:
            structures = json.load(f)
        PTVname = structures['ptv']
        OARnames = structures['oar']
        ROIs = [PTVname, ] + OARnames[1:]  # to exclude skin/body

        # load masks
        MRFolder = os.path.join(dataPatFolder, 'MR')
        rtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=MRFolder, rt_struct_path=rtFile)
        masks = {name: rtstruct.get_roi_mask_by_name(name) for name in ROIs}

        # dose normalization. As the clinical Dose and 
        # BOOdose are not of the same scale
        # here we normalize using D98.
        PTVmask = masks[PTVname]
        BOOPTV = BOOdose[PTVmask]
        clinicalPTV = clinicalDose[PTVmask]
        thresh = 5
        BOOPercentile = np.percentile(BOOPTV, thresh)
        clinicalPercentile = np.percentile(clinicalPTV, thresh)
        clinicalDose = clinicalDose / clinicalPercentile * BOOPercentile
        
        # # draw DVH
        # colorList = DVHgroup(BOOdose, masks)
        # DVHgroup(clinicalDose, masks, linestyle='--', colorList=colorList)
        # visPatFolder = os.path.join(visFolder, patientName)
        # outFile = os.path.join(visPatFolder, 'DVH_PTVcropped{}.png'.format(trailNO))
        # plt.legend(ROIs)
        # plt.savefig(outFile)
        # plt.clf()
        # print(patientName)

        rowIdx = int(i / rowSize)
        colIdx = i % rowSize
        ax = axs[rowIdx, colIdx]
        colorList = DVHgroupAAPM(BOOdose, masks, ax)
        DVHgroupAAPM(clinicalDose, masks, ax, linestyle='--', colorList=colorList)
        # ax.legend(ROIs, bbox_to_anchor=(1.05, 1))
        # ax.legend(ROIs)
        ax.legend(ROIs, loc='upper right', fontsize=labelFontSize)
        ax.set_title(patientName, fontsize=fontSize)
        ax.set_xlabel('Dose [Gy]', fontsize=fontSize)
        ax.set_ylabel('Relative volume [%]', fontsize=fontSize)
        ax.set_ylim([0, 100])
        ax.tick_params(axis='both', labelsize=labelFontSize)
        print(patientName)
    fig.set_figheight(4.0*colSize)
    fig.set_figwidth(6.4*rowSize)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    plt.margins(y=0.3)
    # plt.rcParams.update({'font.size': 40})
    # plt.show()
    outFile = os.path.join(visFolder, 'DVH_AAPM_poster.png')
    plt.savefig(outFile)


def DVH_AUC_comp():
    """
    In preparing the AAPM abstract, I plan to compare the area under 
    curve (AUC) of the DVH figures between the 4-pi plan and clinical plan
    """
    numPatients = 8
    plus = 0
    minus = 0
    negative = 0
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimizePTVcropped')
        # get trailNO
        template = os.path.join(optFolder, 'dose*.npy')
        files = glob.glob(template)
        files.sort()
        BOOdoseFile = files[-1]
        BOOdoseFile_ = BOOdoseFile.split('/')[-1]
        digits = [a for a in BOOdoseFile_ if a.isdigit()]
        trailNO = ''.join(digits)
        trailNO = int(trailNO)
        BOOdose = np.load(BOOdoseFile)

        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        clinicalDose = pydicom.dcmread(MRdoseFile).pixel_array
        clinicalDose = np.transpose(clinicalDose, (1, 2, 0))

        # load structures
        structuresFile = os.path.join(optFolder, 'structures.json')
        with open(structuresFile, 'r') as f:
            structures = json.load(f)
        PTVname = structures['ptv']
        OARnames = structures['oar']
        ROIs = [PTVname, ] + OARnames[1:]  # to exclude skin/body

        # load masks
        MRFolder = os.path.join(dataPatFolder, 'MR')
        rtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=MRFolder, rt_struct_path=rtFile)
        masks = {name: rtstruct.get_roi_mask_by_name(name) for name in ROIs}

        # dose normalization. As the clinical Dose and 
        # BOOdose are not of the same scale
        # here we normalize using D98.
        PTVmask = masks[PTVname]
        BOOPTV = BOOdose[PTVmask]
        clinicalPTV = clinicalDose[PTVmask]
        thresh = 5
        BOOPercentile = np.percentile(BOOPTV, thresh)
        clinicalPercentile = np.percentile(clinicalPTV, thresh)
        clinicalDose = clinicalDose / clinicalPercentile * BOOPercentile

        AUC_BOO = AUCcalc(BOOdose, masks)
        AUC_clinical = AUCcalc(clinicalDose, masks)
        diff = {key: value - AUC_BOO[key] for key, value in AUC_clinical.items()}
        # OAR_diff = {key: value for key, value in diff.items() if 'PTV' not in key}
        OAR_diff = [value for key, value in diff.items() if 'PTV' not in key]
        plus_flag = np.array([a > 0 for a in OAR_diff])
        plus_local = np.sum(plus_flag)
        minus_local = len(plus_flag) - plus_local
        plus = plus + plus_local
        minus = minus + minus_local
        print(patientName, plus_local, minus_local)

        inferior_OARs = [key for key, value in diff.items() if ('PTV' not in key) and (value < 0)]
        print(inferior_OARs)
    print('the number of OARs that BOO outperforms clinical: {}'.format(plus))
    print('the number of OARs that clinical outperforms BOO: {}'.format(minus))
    print('number of total OARs: {}'.format(plus + minus))


def AUCcalc(dose, masks):
    """
    This function calculates the AUC of the DVH
    """
    results = {}
    for name, mask in masks.items():
        doseMasked = dose[mask]
        doseMasked = np.sort(doseMasked)
        doseMasked = np.insert(doseMasked, 0, 0)
        yAxis = 1 - np.arange(len(doseMasked)) / len(doseMasked)
        doseDiff = np.ediff1d(doseMasked, to_end=0)
        result = np.sum(doseDiff*yAxis)
        results[name] = result
    return results


def signTest():
    """
    the function DVH_AUC_comp() gives that BOO results outperform 
    clinical result in 76 out of 83 OARs. This function performs the
    one-tail sign test
    """
    wins = 76
    all = 83
    extremes = 0
    for i in range(wins, all+1):
        product = 1
        for j in range(i+1, all+1):
            product *= j
        for j in range(1, all - i + 1):
            product /= j
        print(product)
        extremes += product
    p_value = extremes / 2 ** all
    print('the p_value equals to {}'.format(p_value))


def PTVHomogeneity():
    """
    This function measures the dose homogeneity in the PTV area
    """
    numPatients = 8
    BOOHomogeneity = []
    clinicalHomogeneity = []
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimizePTVcropped')
        # get trailNO
        template = os.path.join(optFolder, 'dose*.npy')
        files = glob.glob(template)
        files.sort()
        BOOdoseFile = files[-1]
        BOOdoseFile_ = BOOdoseFile.split('/')[-1]
        digits = [a for a in BOOdoseFile_ if a.isdigit()]
        trailNO = ''.join(digits)
        trailNO = int(trailNO)
        BOOdose = np.load(BOOdoseFile)

        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        clinicalDose = pydicom.dcmread(MRdoseFile).pixel_array
        clinicalDose = np.transpose(clinicalDose, (1, 2, 0))

        # load structures
        structuresFile = os.path.join(optFolder, 'structures.json')
        with open(structuresFile, 'r') as f:
            structures = json.load(f)
        PTVname = structures['ptv']
        OARnames = structures['oar']
        ROIs = [PTVname, ] + OARnames[1:]  # to exclude skin/body

        # load masks
        MRFolder = os.path.join(dataPatFolder, 'MR')
        rtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=MRFolder, rt_struct_path=rtFile)
        masks = {name: rtstruct.get_roi_mask_by_name(name) for name in ROIs}

        # dose normalization. As the clinical Dose and 
        # BOOdose are not of the same scale
        # here we normalize using D98.
        PTVmask = masks[PTVname]
        BOOPTV = BOOdose[PTVmask]
        clinicalPTV = clinicalDose[PTVmask]
        thresh = 5
        BOOD95 = np.percentile(BOOPTV, thresh)
        clinicalD95 = np.percentile(clinicalPTV, thresh)
        clinicalDose = clinicalDose / clinicalD95 * BOOD95
        clinicalPTV = clinicalDose[PTVmask]  # update clinicalPTV after the normalization

        # now D95 of the two plans are equal, we measure the D5
        thresh = 95
        BOOD5 = np.percentile(BOOPTV, thresh)
        clinicalD5 = np.percentile(clinicalPTV, thresh)
        BOODHI = BOOD95 / BOOD5
        clinicalDHI = BOOD95 / clinicalD5
        BOOHomogeneity.append(BOODHI)
        clinicalHomogeneity.append(clinicalDHI)
        print(patientName)
    
    # result = np.array([a > b for a, b in zip(BOOHomogeneity, cliincalHomogeneity)])
    # print(np.sum(result))
    print('BOOHomogeneity average: {}, detail: {}'.format(
        np.mean(np.array(BOOHomogeneity)), BOOHomogeneity))
    print('clinicalHomogeneity average: {}, detail: {}'.format(
        np.mean(np.array(clinicalHomogeneity)), clinicalHomogeneity))


def OARDoseAnalysis():
    """
    This function analyze the OAR dose Dmax, D95 etc
    """
    numPatients = 8
    duodDmax = {}
    duodD95 = {}
    stomachDmax = {}
    stomachD95 = {}
    Bwel_SmDmax = {}
    Bwel_SmD95 = {}
    for i in range(numPatients):
        patientName = 'patient{}'.format(i+1)
        expPatFolder = os.path.join(expFolder, patientName)
        optFolder = os.path.join(expPatFolder, 'optimizePTVcropped')
        # get trailNO
        template = os.path.join(optFolder, 'dose*.npy')
        files = glob.glob(template)
        files.sort()
        BOOdoseFile = files[-1]
        BOOdoseFile_ = BOOdoseFile.split('/')[-1]
        digits = [a for a in BOOdoseFile_ if a.isdigit()]
        trailNO = ''.join(digits)
        trailNO = int(trailNO)
        BOOdose = np.load(BOOdoseFile)

        dataPatFolder = os.path.join(anonymousDataFolder, patientName)
        MRdoseFile = os.path.join(dataPatFolder, 'MRdose.dcm')
        clinicalDose = pydicom.dcmread(MRdoseFile).pixel_array
        clinicalDose = np.transpose(clinicalDose, (1, 2, 0))

        # load structures
        structuresFile = os.path.join(optFolder, 'structures.json')
        with open(structuresFile, 'r') as f:
            structures = json.load(f)
        PTVname = structures['ptv']
        OARnames = structures['oar']
        ROIs = [PTVname, ] + OARnames[1:]  # to exclude skin/body

        # load masks
        MRFolder = os.path.join(dataPatFolder, 'MR')
        rtFile = os.path.join(dataPatFolder, 'MRrt.dcm')
        rtstruct = RTStructBuilder.create_from(
            dicom_series_path=MRFolder, rt_struct_path=rtFile)
        masks = {name: rtstruct.get_roi_mask_by_name(name) for name in ROIs}

        # dose normalization. As the clinical Dose and 
        # BOOdose are not of the same scale
        # here we normalize using D98.
        PTVmask = masks[PTVname]
        BOOPTV = BOOdose[PTVmask]
        clinicalPTV = clinicalDose[PTVmask]
        thresh = 5
        BOOD95 = np.percentile(BOOPTV, thresh)
        clinicalD95 = np.percentile(clinicalPTV, thresh)
        clinicalDose = clinicalDose / clinicalD95 * BOOD95
        clinicalPTV = clinicalDose[PTVmask]  # update clinicalPTV after the normalization

        # calculate the statistics
        thresh = 5
        duodMask = None
        for name in OARnames:
            if 'duod' in name.lower():
                duodMask = masks[name]
                print(patientName, name)
                break
        if duodMask is not None:
            BOOduodDose = BOOdose[duodMask]
            clinicalduodDose = clinicalDose[duodMask]
            duodDmax[patientName] = \
                (np.max(BOOduodDose), np.max(clinicalduodDose))
            duodD95[patientName] = \
                (np.percentile(BOOduodDose, thresh), np.percentile(clinicalduodDose, thresh))

        stmcMask = None
        for name in OARnames:
            if 'stmc' in name.lower():
                stmcMask = masks[name]
                print(patientName, name)
                break
        if stmcMask is not None:
            BOOstmcDose = BOOdose[stmcMask]
            clinicalstmcDose = clinicalDose[stmcMask]
            stomachDmax[patientName] = (np.max(BOOstmcDose), np.max(clinicalstmcDose))
            stomachD95[patientName] = \
                (np.percentile(BOOstmcDose, thresh), np.percentile(clinicalstmcDose, thresh))
        
        Bwel_SmMask = None
        for name in OARnames:
            if 'bwel' in name.lower() and 'sm' in name.lower():
                Bwel_SmMask = masks[name]
                print(patientName, name)
                break
        if Bwel_SmMask is not None:
            BOO_Bwel_Sm_dose = BOOdose[Bwel_SmMask]
            clinical_Bwel_Sm_dose = clinicalDose[Bwel_SmMask]
            Bwel_SmDmax[patientName] = \
                (np.max(BOO_Bwel_Sm_dose), np.max(clinical_Bwel_Sm_dose))
            Bwel_SmD95[patientName] = \
                (np.percentile(BOO_Bwel_Sm_dose, thresh), np.percentile(clinical_Bwel_Sm_dose, thresh))

        print('\n')
    
    print(duodDmax)
    print(duodD95, '\n')
    print(stomachDmax)
    print(stomachD95, '\n')
    print(Bwel_SmDmax)
    print(Bwel_SmD95, '\n')


def drawTable():
    """
    This funciton draws the table for latex supporting document
    """
    if False:
        BOO_PTV_homogeneity = [0.8069704537975423, 0.831610374097704, 0.8141663893520457, 0.8477902369437897, \
            0.8152625310994703, 0.7965236405731738, 0.829527639796673, 0.8188531728013859]
        clinical_PTV_homogeneity = [0.8243131944603638, 0.6392917050717626, 0.6016216690477412, 0.7431962864940729, \
            0.5875315261144793, 0.652439080350446, 0.842159916926272, 0.7484852896893811]
        BOOPTVHomoLine = ''
        clinicalPTVHomoLine = ''
        for a, b in zip(BOO_PTV_homogeneity, clinical_PTV_homogeneity):
            if a <= b:
                BOOPTVHomoLine = BOOPTVHomoLine + '{:.2f}'.format(a) + ' & '
                clinicalPTVHomoLine = clinicalPTVHomoLine + '\\textbf{' + '{:.2f}'.format(b) + '} & '
            else:
                BOOPTVHomoLine = BOOPTVHomoLine + '\\textbf{' + '{:.2f}'.format(a) + '} & '
                clinicalPTVHomoLine = clinicalPTVHomoLine + '{:.2f}'.format(b) + ' & '
        BOOPTVHomoLine = BOOPTVHomoLine[:-2] + '\\\\'
        clinicalPTVHomoLine = clinicalPTVHomoLine[:-2] + '\\\\'
        print(BOOPTVHomoLine)
        print(clinicalPTVHomoLine)

    if False:
        # work on duodenum data
        DmaxData = {'patient1': (4.434961, 9.542294145400604), 'patient3': (20.11872, 16.32372906512959), \
            'patient4': (5.9481754, 7.462174300960876), 'patient5': (20.52946, 19.007588592070096), \
            'patient6': (20.721945, 18.105720650842812), 'patient7': (8.374871, 9.737669512265825), \
            'patient8': (20.050129, 15.89796739737968)}
        D95Data = {'patient1': (0.007169494638219476, 0.36532869675322555), 
            'patient3': (0.6078828811645508, 1.459492540149958), 
            'patient4': (0.07187967896461486, 0.5483072557735594), 
            'patient5': (0.5204100906848907, 1.092840708675186), 
            'patient6': (0.08204947113990785, 0.5828413865801028), 
            'patient7': (0.09270827397704125, 0.32978152364548996), 
            'patient8': (0.14333276748657225, 1.4260327845618623)}
    
    if False:
        # work on stomach data
        DmaxData = {'patient1': (20.067997, 12.925583943839994), 
            'patient2': (18.924007, 17.02968003032645), 
            'patient3': (19.78814, 16.71740797398583), 
            'patient4': (19.174734, 17.6213070551169), 
            'patient5': (5.796191, 3.612891711912347), 
            'patient6': (19.440487, 18.787120983234924), 
            'patient7': (20.281845, 14.27970569320833), 
            'patient8': (19.046717, 17.42650469676143)}
        D95Data = {'patient1': (0.11163704693317414, 0.8952022220709066), 
            'patient2': (0.3837432414293289, 0.7850860534687265), 
            'patient3': (0.07628622315824034, 0.46957983669611747), 
            'patient4': (1.2377982139587402, 0.531475283666129), 
            'patient5': (0.04227860681712627, 0.5072216599469099), 
            'patient6': (0.03479145020246506, 0.7260481190216057), 
            'patient7': (0.11516406089067459, 0.4199328547324983), 
            'patient8': (0.3579861521720886, 0.9727055541520212)}
    
    if True:
        # work on small intestine
        DmaxData = {'patient1': (16.483133, 13.289606729255407), 
            'patient2': (19.426186, 16.56037855567829), 
            'patient7': (11.974609, 14.225349743582342), 
            'patient8': (0.43506527, 0.6505590386806811)}
        D95Data = {'patient1': (0.007545600971207023, 0.23898177481980437), 
            'patient2': (0.3293981209397316, 0.3452779726002546), 
            'patient7': (0.059619066119194035, 0.3142039039356025), 
            'patient8': (0.0756531871855259, 0.4031254507677399)}
    
    # work on Dmax first
    fourPiDmax = ''
    clinicalDmax = ''
    for patient, data in DmaxData.items():
        BOOdata = data[0]
        clinicalData = data[1]
        if BOOdata < clinicalData:
            fourPiDmax = fourPiDmax + '\\textbf{' + '{:.2f}'.format(BOOdata) + '} & '
            clinicalDmax = clinicalDmax + '{:.2f}'.format(clinicalData) + ' & '
        else:
            fourPiDmax = fourPiDmax + '{:.2f}'.format(BOOdata) + ' & '
            clinicalDmax = clinicalDmax + '\\textbf{' + '{:.2f}'.format(clinicalData) + '} & '
    
    # then D95
    fourPiD95 = ''
    clinicalD95 = ''
    for patient, data in D95Data.items():
        BOOdata = data[0]
        clinicalData = data[1]
        if BOOdata < clinicalData:
            fourPiD95 = fourPiD95 + '\\textbf{' + '{:.2f}'.format(BOOdata) + '} & '
            clinicalD95 = clinicalD95 + '{:.2f}'.format(clinicalData) + ' & '
        else:
            fourPiD95 = fourPiD95 + '{:.2f}'.format(BOOdata) + ' & '
            clinicalD95 = clinicalD95 + '\\textbf{' + '{:.2f}'.format(clinicalData) + '} & '
    print(fourPiDmax)
    print(clinicalDmax)
    print(fourPiD95)
    print(clinicalD95)


def BOObeamFuse():
    """
    Last night, my Ph.D. supervisor asked me to add some plots of beams
    of representative cases. I generated these plots using matlab. This
    function is to fuse them into one image
    """
    numPatients = 8
    patientList = [1, 2, 3, 4]
    targetShape = (1380, 960)
    croppedImages = []
    for i in patientList:
        patientName = 'patient{}'.format(i)
        visPatFolder = os.path.join(visFolder, patientName)
        imageFile = os.path.join(visPatFolder, 'BOObeams.png')
        image = plt.imread(imageFile)
        # print(patientName, image.shape)
        
        shape = image.shape
        upperBound = 0
        lowerBound = targetShape[0]
        leftBound = int((shape[1] - targetShape[1])/2)
        rightBound = leftBound + targetShape[1]
        cropped = image[upperBound: lowerBound, leftBound: rightBound]
        croppedImages.append(cropped)
        print(patientName)
    result = np.concatenate(croppedImages, axis=1)
    outputPath = os.path.join(visFolder, 'AAPMSupBOOBeams.png')
    plt.imsave(outputPath, result)


if __name__ == '__main__':
    # reviseStructures()
    # createCTanatomy()
    # doseResizeCropped()
    # visDVH_PTVcropped()
    prepareAAPMSuppDoc()
    # DVH_AUC_comp()
    # signTest()
    # PTVHomogeneity()
    # OARDoseAnalysis()
    # drawTable()
    # BOObeamFuse()