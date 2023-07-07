import os
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import math

# below are the material properties.
# density are in g/cm^3
materials = {'adipose': [0.95, {'H': 0.114, 'C': 0.598, 
        'N': 0.007, 'O': 0.278, 'Na': 0.001, 'S': 0.001, 'Cl': 0.001, }],
    'bone': [1.85, {'H': 0.064, 'C': 0.278, 'N': 0.027, 'O': 0.41, 
        'Mg': 0.002, 'P': 0.07, 'S': 0.002, 'Ca': 0.147}],
    'lung': [1.040, {'H': 0.105, 'C': 0.083, 'N': 0.023, 
        'O': 0.779, 'Na': 0.002, 'P': 0.001, 'S': 0.002, 'Cl': 0.003, 'K': 0.002}],
    'muscle': [1.05, {'H': 0.102, 'C': 0.143, 'N': 0.034, 
        'O': 0.71, 'Na': 0.001, 'P': 0.002, 'S': 0.003, 'Cl': 0.001, 'K': 0.004}]}

# data are under 6MeV, units are in cm^2/g
TermaTable = {'H': 4.498e-2, 'C': 2.469e-2, 'N': 2.511e-02, 'O': 2.552e-02,
    'Na': 2.559e-02, 'Mg': 2.681e-02, 'P': 2.747e-02, 'S': 2.872e-02,
    'Cl': 2.798e-02, 'K': 2.915e-02, 'Ca':3.035e-02}

muTable = None

def muTableInit():
    """
    This function initializes the muTable
    """
    # firstly, check the unity of the material decomposition
    for key, item in materials.items():
        elementTable = item[1]
        fractions = list(elementTable.values())
        unity = np.sum(fractions)
        assert np.abs(unity - 1) < 1e-3, "the fractions do not add to unity"

    global muTable
    muTable = {}
    for key, value in materials.items():
        density = value[0]
        elements = value[1]
        mu = 0
        for element, fraction in elements.items():
            elementDensity = density * fraction
            elementMuRho = TermaTable[element]
            elementMu = elementDensity * elementMuRho
            mu += elementMu
        muTable[key] = mu

muTableInit()


def doseRead():
    """
    This function reads the dose
    """
    doseFolder = '/data/qifan/projects/EndtoEnd4/results' \
        '/InhomoJuly6/phantomDose'
    doseFile = os.path.join(doseFolder, 'array.bin')
    shape = (199, 199, 256)
    res = (0.1, 0.1, 0.1 )  # cm
    edep = np.fromfile(doseFile, dtype=np.double)
    edep = np.reshape(edep, shape)

    # construct density matrix
    density = np.zeros(shape)
    thicknesses = [16, 16, 16, 16, 96, 16, 16, 16, 16, 16, 16]
    MM = ['adipose', 'muscle', 'bone', 'muscle', 'lung', 'muscle', 
        'bone', 'adipose', 'bone', 'muscle', 'adipose']
    prev = 0
    for thick, mat in zip(thicknesses, MM):
        density[:, :, prev:prev+thick] = materials[mat][0]
        prev += thick
    
    depth = np.arange(shape[2]) * res[2]
    dose = edep / density
    file = os.path.join(doseFolder, 'IPB.png')
    plt.plot(depth, dose[99, 99, :])
    plt.savefig(file)
    plt.clf()

    # then, we simulate the 0.5cm, 1cm, and 2cm data
    kernelSizes = [5, 10, 20]
    for ks in kernelSizes:
        kernelSizePrecessing(dose, shape, ks, doseFolder)
    


def kernelSizePrecessing(dose, shape, kernelSize, doseFolder):
    newShape = (shape[0] + kernelSize - 1, shape[1] + kernelSize - 1, shape[2])
    newDose = np.zeros(newShape, dtype=np.double)
    kernel = np.ones((kernelSize, kernelSize))
    for i in range(shape[2]):
        newDose[:, :, i] = signal.convolve2d(dose[:, :, i], kernel)

    centerIdx = int((shape[0] + kernelSize) / 2 - 1)
    centralLine = newDose[centerIdx, centerIdx, :]
    centralLine = centralLine / np.max(centralLine)
    res = 0.1
    depth = np.arange(shape[2]) * res
    plt.plot(depth, centralLine)
    plt.xlabel('depth (cm)')
    plt.ylabel('dose (a.u.)')
    plt.title('field size {} cm'.format(kernelSize*res))
    file = os.path.join(doseFolder, 'window{}.png'.format(kernelSize*res))
    plt.savefig(file)
    plt.clf()


def waterBoneComp():
    """
    This function compares the dose distribution within water and bone
    """
    waterFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/waterDose'
    boneFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/boneDose'
    outputFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/phantomDose'
    
    shape = (199, 199, 256)
    waterEdepFile = os.path.join(waterFolder, 'array.bin')
    waterEdep = np.fromfile(waterEdepFile, np.double)
    waterEdep = np.reshape(waterEdep, shape)
    waterDensity = 1
    waterDose = waterEdep / waterDensity

    boneEdepFile = os.path.join(boneFolder, 'array.bin')
    boneEdep = np.fromfile(boneEdepFile, dtype=np.double)
    boneEdep = np.reshape(boneEdep, shape)
    boneDensity = 1.85
    boneDose = boneEdep / boneDensity

    res = 0.1
    centerIdx = 99
    waterDepth = np.arange(shape[2]) * res

    if True:
        # plot the naive dose comparison
        plt.plot(waterDepth, waterDose[centerIdx, centerIdx, :])
        plt.plot(waterDepth, boneDose[centerIdx, centerIdx, :])
        plt.legend(['water', 'bone'])
        plt.xlabel('depth')
        plt.ylabel('dose (a.u.)')
        figureFile = os.path.join(outputFolder, 'waterBoneVanilla.png')
        plt.savefig(figureFile)
        plt.clf()
    
    if True:
        # plot the scaled dose comparison
        boneDepth = np.arange(shape[2]) * boneDensity * res
        plt.plot(waterDepth, waterDose[centerIdx, centerIdx, :])
        plt.plot(boneDepth, boneDose[centerIdx, centerIdx, :])
        plt.legend(['water', 'bone'])
        plt.xlabel('depth')
        plt.ylabel('dose (a.u.)')
        figureFile = os.path.join(outputFolder, 'waterBoneCentralRescaled.png')
        plt.savefig(figureFile)
        plt.clf()
    
    if True:
        # plot the partial dose
        waterPartialDose = np.sum(waterDose, axis=(0, 1))
        bonePartialDose = np.sum(boneDose, axis=(0, 1))
        boneDepth = np.arange(shape[2]) * boneDensity * res
        plt.plot(waterDepth, waterPartialDose)
        plt.plot(boneDepth, bonePartialDose)
        plt.legend(['water', 'bone'])
        plt.xlabel('depth')
        plt.ylabel('dose (a.u.)')
        figureFile = os.path.join(outputFolder, 'waterBonePartialRescaled.png')
        plt.savefig(figureFile)
        plt.clf()


def transverseProfile():
    """
    This function shows the transverse profile of the dose
    """
    waterFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/waterDose'
    boneFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/boneDose'
    outputFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/phantomDose'
    
    shape = (199, 199, 256)
    waterEdepFile = os.path.join(waterFolder, 'array.bin')
    waterEdep = np.fromfile(waterEdepFile, np.double)
    waterEdep = np.reshape(waterEdep, shape)
    waterDensity = 1.
    waterDose = waterEdep / waterDensity

    boneEdepFile = os.path.join(boneFolder, 'array.bin')
    boneEdep = np.fromfile(boneEdepFile, dtype=np.double)
    boneEdep = np.reshape(boneEdep, shape)
    boneDensity = 1.85
    boneDose = boneEdep / boneDensity

    res = 0.1
    binScale = 10000
    depth = np.arange(shape[2]) * res
    # we only show dose profile within this range
    radialCutoff = 2.0
    # radialCutoff = 0.6

    # here we calculate the depth at the following depths:
    # 2cm, 5cm, 10cm, 15cm, 20cm
    depths = [20, 50, 100, 150, 200]
    if True:
        water2Show = []
        for depth in depths:
            waterSlice = waterDose[:, :, depth]

            waterDistance, waterProfile = slice2Tranverse(waterSlice, binScale)
            physicalDistance = np.array(waterDistance) / binScale * res

            # get the dose profile cutoff idx
            cutoffIdx = 0
            for i, r in enumerate(physicalDistance):
                if r > radialCutoff:
                    cutoffIdx = i
                    break
            
            distance2Show = physicalDistance[:cutoffIdx]
            waterProfile2Show = np.array(waterProfile[:cutoffIdx])
            hProfile = distance2Show * waterProfile2Show
            water2Show.append([distance2Show, hProfile])
            
        if False:
            # firstly, plot the water h profile
            for x_axis, y_axis in water2Show:
                plt.plot(x_axis[1:], y_axis[1:])
            plt.xlabel('distance (cm)')
            plt.ylabel('h profile (a.u.)')
            plt.legend(['2cm', '5cm', '10cm', '15cm', '20cm'])
            figureFile = os.path.join(outputFolder, 'hprofile_water.png')
            plt.savefig(figureFile)
            plt.clf()

        if False:
            # secondly, plot the normalized h profile
            for x_axis, y_axis in water2Show:
                y_axis_ = y_axis / np.max(y_axis)
                plt.plot(x_axis[1:], y_axis_[1:])
            plt.xlabel('distance (cm)')
            plt.ylabel('h profile (a.u.)')
            plt.legend(['2cm', '5cm', '10cm', '15cm', '20cm'])
            figureFile = os.path.join(outputFolder, 'hprofile_water_normalized.png')
            plt.savefig(figureFile)
            plt.clf()
    
    if True:
        bone2Show = []
        for depth in depths:
            # convert depth in water to depth in bone
            boneDepth = int(depth * waterDensity / boneDensity)
            boneSlice = boneDose[:, :, boneDepth]

            boneDistance, boneProfile = slice2Tranverse(boneSlice, binScale)
            physicalDistance = np.array(boneDistance) / binScale * res

            # get the dose profile cutoff idx
            cutoffIdx = 0
            for i, r in enumerate(physicalDistance):
                if r > radialCutoff:
                    cutoffIdx = i
                    break
            
            distance2Show = physicalDistance[:cutoffIdx]
            boneProfile2Show = np.array(boneProfile[:cutoffIdx])
            hProfile = distance2Show * boneProfile2Show
            bone2Show.append([distance2Show, hProfile])
        
        if False:
            # firstly, plot the water h profile
            for x_axis, y_axis in bone2Show:
                plt.plot(x_axis[1:], y_axis[1:])
            plt.xlabel('distance (cm)')
            plt.ylabel('h profile (a.u.)')
            plt.legend(['2cm', '5cm', '10cm', '15cm', '20cm'])
            figureFile = os.path.join(outputFolder, 'hprofile_bone.png')
            plt.savefig(figureFile)
            plt.clf()
        
        if False:
            # secondly, plot the normalized h profile
            for x_axis, y_axis in bone2Show:
                y_axis_ = y_axis / np.max(y_axis)
                plt.plot(x_axis[1:], y_axis_[1:])
            plt.xlabel('distance (cm)')
            plt.ylabel('h profile (a.u.)')
            plt.legend(['2cm', '5cm', '10cm', '15cm', '20cm'])
            figureFile = os.path.join(outputFolder, 'hprofile_bone_normalized.png')
            plt.savefig(figureFile)
            plt.clf()
    
    # compare water dose and bone dose
    # here we only take 10cm depth for example
    waterDistance10, waterProfile10 = water2Show[2]
    boneDistance10, boneProfile10 = bone2Show[2]
    # normalize dose
    roof = np.max(waterProfile10)
    waterProfile10 = waterProfile10 / roof
    boneProfile10 = boneProfile10 / roof

    # normalize distance
    boneDistance10 = boneDistance10 * boneDensity

    if False:
        plt.plot(waterDistance10[1:], waterProfile10[1:])
        plt.plot(boneDistance10[1:], boneProfile10[1:])
        plt.legend(['water', 'bone'])
        plt.xlabel('distance (cm)')
        plt.ylabel('h profile (a.u.)')
        plt.title('water and bone $h$ profile comparison')
        figureFile = os.path.join(outputFolder, 'waterBoneCompNorm.png')
        plt.savefig(figureFile)
    
    # then we calculate the ratio between water dose and bone dose
    boneDistance10_ = []
    for a in boneDistance10:
        if a < waterDistance10[-1]:
            boneDistance10_.append(a)
    boneDistance10_ = np.array(boneDistance10_)
    waterInterpolated = np.interp(boneDistance10_, waterDistance10, waterProfile10)
    
    ratio = boneProfile10[:len(boneDistance10_)] / waterInterpolated
    plt.plot(boneDistance10_, ratio)
    figureFile = os.path.join(outputFolder, 'waterBoneRatio.png')
    plt.savefig(figureFile)
    plt.clf()


def slice2Tranverse(slice, binScale=10000):
    """
    This function converts the 2d slice to a linear dose spread function, 
    as the dose distribution is symmetric
    """
    # we assume that the shape of the slice is square whose edge size is 
    # odd. So that the central pixel is the origin
    shape = slice.shape
    centerIdx = int((1 + shape[0]) / 2) - 1
    statistics = {}
    datapoints = {}
    for i in range(shape[0]):
        for j in range(shape[1]):
            distance = math.sqrt((i-centerIdx)**2 + (j-centerIdx)**2)
            distanceRescale = int(binScale*distance)
            if distanceRescale in statistics:
                statistics[distanceRescale] += slice[i, j]
                datapoints[distanceRescale] += 1
            else:
                statistics[distanceRescale] = slice[i, j]
                datapoints[distanceRescale] = 1
    # averaging
    for key in statistics:
        statistics[key] /= datapoints[key]
    # convert to list
    statistics_ = [(a, b) for a, b in statistics.items()]
    statistics_.sort(key=lambda a : a[0])
    distance = [a[0] for a in statistics_]
    dose = [a[1] for a in statistics_]
    return distance, dose


if __name__ == '__main__':
    doseRead()
    # waterBoneComp()
    # transverseProfile()