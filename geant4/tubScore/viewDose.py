import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

phantom1 = [("adipose", 0.8),
            ("muscle", 0.8),
            ("bone", 0.8),
            ("muscle", 0.8),
            ("lung", 4.8),
            ("muscle", 0.8),
            ("bone", 0.8),
            ("adipose", 0.8),
            ("bone", 0.8),
            ("muscle", 0.8),
            ("adipose", 0.8)]

density = {"adipose": 0.92, "muscle": 1.04, "bone": 1.85, "lung": 0.25}

def getArray(folder, shape=None):
    """
    this function gets the array in the folder
    """
    file = os.path.join(folder, 'SD.bin')
    array = np.fromfile(file, dtype=np.float64)
    if shape is None:
        shape = (256, 200)
    array = np.reshape(array, shape)
    return array


def getDensity(shape, resZ=0.05):
    """
    This function gets the density of the phantom
    """
    den = np.zeros(shape)
    count = 0
    resZ = 0.05
    for material, thickness in phantom1:
        d = density[material]
        dim = round(thickness / resZ)
        for i in range(count, count+dim):
            den[i, :] = d
        count += dim
    return den


def readDose():
    """
    This function reads the dose from the Monte Carlo simulation
    """
    resultFolder = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14'
    array = getArray(resultFolder)
    shape = array.shape

    # normalize the dose array against the density
    den = getDensity(shape)
    
    dose = array / den
    resultFolder = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14'
        
    # view central dose
    slice = dose[:, 0]
    depth = np.arange(shape[0]) * 0.1  # cm
    plt.plot(depth, slice)
    plt.xlabel('depth / cm')
    plt.ylabel('dose (a.u.)')
    plt.title('centerline dose')
    file = os.path.join(resultFolder, 'centerlineDose.png')
    plt.savefig(file)
    plt.clf()

    # view partial dose
    partialDose = np.sum(dose, axis=1)
    plt.plot(depth, partialDose)
    plt.xlabel('depth / cm')
    plt.ylabel('(a.u.)')
    plt.title('partial dose')
    file = os.path.join(resultFolder, 'partialDose.png')
    plt.savefig(file)
    plt.clf()


def cylindricalToCartesian(array_cylindrical):
    """
    This function converts the dose distribution in a cylindrical grid to a cartisian grid
    We assume the result, array_cartesian is of the same dimension as the the input 
    array_cylindrical
    """
    shape_cylindrical = array_cylindrical.shape
    nLayers, nRings = shape_cylindrical

    # we intend to use up sampling strategy.
    up_sampling = 10
    res = 1. / up_sampling
    halfGridSize =  nRings * up_sampling
    matrixSize = (nLayers, halfGridSize, halfGridSize)
    matrix = np.zeros(matrixSize, dtype=array_cylindrical.dtype)

    for i in range(halfGridSize):
        x = (i + 0.5) * res
        xSquare = x ** 2
        for j in range(halfGridSize):
            y = (j + 0.5) * res
            ySquare = y ** 2
            rSquare = xSquare + ySquare
            r = np.sqrt(rSquare)
            r_idx = int(np.floor(r))
            if r_idx >= nRings:
                continue
            matrix[:, i, j] = array_cylindrical[:, r_idx]
        print('progress: {}/{}'.format(i+1, halfGridSize))
    
    # then, we downsample the matrix to get a quarter of the result
    quarter = np.zeros((nLayers, nRings, nRings), dtype=array_cylindrical.dtype)
    for i in range(nRings):
        for j in range(nRings):
            sample = matrix[:, i*up_sampling:(i+1)*up_sampling, j*up_sampling:(j+1)*up_sampling]
            quarter[:, i, j] = np.sum(sample, axis=(1, 2))
    
    # then, we concatenate the matrix to get the full result
    secondQuarter = np.flip(quarter, axis=1)
    quarter12 = np.concatenate((secondQuarter, quarter), axis=1)
    quarter34 = np.flip(quarter12, axis=2)
    result = np.concatenate((quarter34, quarter12), axis=2)
    return result.copy()


def testCylindricalToCartesian():
    """
    This function tests the functionality of the function above
    """

    resultFolder = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14'
    array = getArray(resultFolder)
    shape = array.shape

    # the array here is the energy deposition. Firstly, we convert it to dose
    # Firstly, we calculate the volume of the individual rings
    resR = 0.05  # cm
    resZ = 0.1  # cm
    nRings = shape[1]
    volumes = np.zeros(nRings)
    RoldSquare = 0.
    for i in range(nRings):
        Rnew = (i + 1) * resR
        RnewSquare = Rnew ** 2
        volumes[i] = np.pi * (RnewSquare - RoldSquare) * resZ
        RoldSquare = RnewSquare
    
    # Then, we multiplies the volume matrix with the density matrix to get the mass matrix
    den = getDensity(shape)
    mass = np.zeros_like(den)
    for i in range(shape[0]):
        mass[i, :] = den[i, :] * volumes
    
    dose = array / mass

    if False:
        # take a look at the dose, we should show the centerline dose
        centerline = dose[:, 0]
        depth = np.arange(shape[0]) * resZ
        plt.plot(depth, centerline)
        plt.show()

        # we then take a look at the dose along the radius dimension
        line = dose[100, :]
        radius = np.arange(shape[1]) * resR
        plt.plot(radius, line)
        plt.show()

        # we then take a look at the energy deposition
        line = array[100, :]
        plt.plot(radius, line)
        plt.show()
    
    dose_Cartesian = cylindricalToCartesian(dose)
    resultFolder = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14'
    file = os.path.join(resultFolder, 'doseCartisian.npy')
    np.save(file, dose_Cartesian)


def TakeALook6MeV():
    """
    In the test function above, we tested the functionality to convert 
    a cylindrical matrix to a cartesian matrix. Now we take a look at the
    cartesian matrix
    """
    file = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14/doseCartisian.npy'
    array = np.load(file)

    # slice profile
    slice = array[100, :, :]
    plt.imshow(slice)
    plt.show()

    # line profile
    line = slice[200, :]
    xAxis = np.arange(400) - 199.5
    plt.plot(xAxis, line)
    plt.show()

    # partial dose
    partialDose = np.sum(array, axis=(1, 2))
    depth = np.arange(256) * 0.1
    plt.plot(depth,partialDose)
    plt.show()


def polyChromatic():
    """
    Here we combine the different monochromatic kernels to 
    form a polychromatic kernel up to some coefficients.
    """
    dataFolder = '/data/qifan/projects/EndtoEnd/results/spec6MV'
    file = os.path.join(dataFolder, 'spec_6mv.spec')
    validLines = 14
    validCols = 4
    with open(file, 'r') as f:
        lines = f.readlines()
    lines = lines[:validLines]
    parseMat = np.zeros((validLines, validCols))
    for i, line in enumerate(lines):
        line = line.split(' ')
        for j in range(validCols):
            parseMat[i, j] = float(line[j])

    if False:
        # just to take take a look at this data.
        # the first colume is photon energy,
        # the second colume is energy fluence,
        # the third colume is mass attenuation coefficient
        # the forth colume is mass energy transfer coefficient
        print(parseMat)
        colume0 = parseMat[:, 0]
        colume1 = parseMat[:, 1]
        sum = np.sum(colume1)
        print('sum the second colume: {}'.format(sum))

        colume2 = parseMat[:, 2]
        sum = np.sum(colume2)
        print('sum the third colume: {}'.format(sum))

        colume3 = parseMat[:, 3]
        sum = np.sum(colume3)
        print('sum the forth colume: {}'.format(sum))

        testColume1 = count2Energy(colume1, colume0)
        print(testColume1)
        testColume2 = count2Energy(colume2, colume0)
        print(testColume2)
    
    energy = parseMat[:, 0]
    fluence = parseMat[:, 1]
    photonCount = fluence / energy
    photonCountSum = np.sum(photonCount)
    photonCount = photonCount / photonCountSum
    
    dataFolder = '/data/qifan/projects/EndtoEnd/results/spec6MV'
    energies = "0.2 0.3 0.4 0.5 0.6 0.8 1.00 1.25 1.50 2.00 3.00 4.00 5.00 6.00"
    energies = energies.split(' ')
    nSlices = 256
    nRings = 200
    shape = (nSlices, nRings)
    Array = np.zeros(shape, dtype=np.float64)
    for i, en in enumerate(energies):
        folder = os.path.join(dataFolder, 'E{}'.format(en))
        array = getArray(folder)
        Array += array * photonCount[i]
    
    resR = 0.05
    resZ = 0.1
    volumes = np.zeros(nRings)
    RoldSquare = 0.
    for i in range(nRings):
        Rnew = (i + 1) * resR
        RnewSquare = Rnew ** 2
        volumes[i] = np.pi * (RnewSquare - RoldSquare) * resZ
        RoldSquare = RnewSquare
    
    den = getDensity(shape)
    mass = np.zeros_like(den)
    for i in range(shape[0]):
        mass[i, :] = den[i, :] * volumes
    
    doseCylindrical = Array / mass
    doseCartesian = cylindricalToCartesian(doseCylindrical)
    outputFile = '/data/qifan/projects/EndtoEnd/results/spec6MV/polyIPB.npy'
    np.save(outputFile, doseCartesian)


def count2Energy(count, energy):
    """
    This function converts the photoncount to energy deposition
    count is photon count, energy is single photon energy
    """
    semi = count * energy
    total = np.sum(semi)
    result = semi / total
    return result


def viewPoly():
    resultFolder = '/data/qifan/projects/EndtoEnd/results/spec6MV'
    outputFile = os.path.join(resultFolder, 'polyIPB.npy')
    polyDose = np.load(outputFile)

    # to take a look
    if False:
        centerLine = polyDose[:, 199:201, 199:201]
        centerLine = np.sum(centerLine, axis=(1, 2))
        depth = np.arange(256) * 0.1
        plt.plot(depth, centerLine)
        plt.show()
    
    if False:
        # try to colvolve the kernel with native python code
        # however, too slow
        nSlices = 256
        dimOrg = 400
        # window size: 0.5cm, 1cm, 2cm, corresponding to 10, 20, 40 pixels
        pixels = [10, 20, 40]
        for p in pixels:
            dimPad = dimOrg + p - 1
            outputShape = (nSlices, dimPad, dimPad)
            kernel = np.ones((p, p), dtype=polyDose.dtype)
            output = np.zeros(outputShape, dtype=polyDose.dtype)
            for i in range(nSlices):
                output[i, :, :] = signal.convolve2d(kernel, polyDose[i, :, :])
                print('progress: {}/{}'.format(i+1, nSlices))
            file = os.path.join(resultFolder, 'kernel{}.npy'.format(p))
            np.save(file, output)
    
    if False:
        # export it to binary file to be used by C++ code
        outputFile = os.path.join(resultFolder, 'polyIPB.bin')
        polyDose.tofile(outputFile)
        print('shape: {}'.format(polyDose.shape))
    
    nSlices = 256
    dimOrg = 400
    # window size: 0.5cm, 1cm, 2cm, corresponding to 10, 20, 40 pixels
    pixels = [10, 20, 40]
    for p in pixels:
        # firstly, convolve along axis 1
        result1 = np.zeros((nSlices, dimOrg+p-1, dimOrg), dtype=np.float64)
        convKernel1d = np.ones((p, 1), dtype=np.float64)
        for i in range(nSlices):
            result1[i, :, :] = signal.convolve2d(polyDose[i, :, :], convKernel1d)
            print('progress 1: {}/{}'.format(i+1, nSlices))
        
        # secondly, convolve along axis 2
        result2 = np.zeros((nSlices, dimOrg+p-1, dimOrg+p-1), dtype=np.float64)
        convKernel1d = np.ones((1, p), dtype=np.float64)
        for i in range(nSlices):
            result2[i, :, :] = signal.convolve2d(result1[i, :, :], convKernel1d)
            print('progress 2: {}/{}'.format(i+1, nSlices))
        
        file = os.path.join(resultFolder, 'window{}.npy'.format(p))
        np.save(file, result2)


def viewCenterline():
    """
    This function views the centerline dose
    """
    resultFolder = '/data/qifan/projects/EndtoEnd/results/spec6MV'
    windowSizes = [10, 20, 40]
    windowRes = 0.05
    nSlices = 256
    resZ = 0.1
    depth = np.arange(nSlices) * resZ
    for w in windowSizes:
        file = os.path.join(resultFolder, 'window{}.npy'.format(w))
        data = np.load(file)
        dmax = np.max(data)
        data /= dmax
        dimNow = data.shape[1]
        centerLine = round((dimNow-1)/2)
        line = data[:, centerLine, centerLine]
        plt.plot(depth, line)
        plt.xlabel('depth / cm')
        plt.ylabel('centerline dose (a.u.)')
        plt.title('6MV, window size {}'.format(w*windowRes))
        file = os.path.join(resultFolder, 'window{}.png'.format(w*windowRes))
        plt.savefig(file)
        plt.clf()

        if False:
            # in the figure below, we show the whole range
            idxZ = 100
            lateral = data[idxZ, centerLine, :]
            xAxis = (np.arange(dimNow, dtype=np.float64)) - centerLine
            xAxis *= windowRes
            plt.plot(xAxis, lateral)
            plt.xlabel('off-axis distance / cm')
            plt.ylabel('lateral dose (a.u.)')
            plt.title('6MV, window size {}'.format(w*windowRes))
            file = os.path.join(resultFolder, 'window{}lateral.png'.format(w*windowRes))
            plt.savefig(file)
            plt.clf()
        
        # now, we only show a window size of 2.5 cm
        validWindow = round(5.0 / 0.05) + 1
        margin = round((dimNow - validWindow) / 2)
        idxZ = 100
        lateral = data[idxZ, centerLine, margin:margin+validWindow]
        lateral = lateral.copy()
        # normalize
        dmax = np.max(data)
        xAxis = (np.arange(validWindow) - int((validWindow-1)/2)) * windowRes
        plt.plot(xAxis, lateral)
        plt.xlabel('off-axis distance / cm')
        plt.ylabel('dose (a.u.)')
        # plt.ylim(0.55)
        file = os.path.join(resultFolder, 'window{}lateral.png'.format(w*windowRes))
        plt.savefig(file)
        plt.clf()


def homoDose():
    """
    We got the dose distribution in homogeneous phantoms.
    Here, we are going to test our hypothesis that the dose 
    at the same radiological coordinate is proportional to 
    density square
    """
    figureFolder = './figures'
    if not os.path.isdir(figureFolder):
        os.mkdir(figureFolder)
    waterDensity = 1.0
    boneDensity = 1.85
    resR = 0.05
    resZ = 0.1
    waterResult='/data/qifan/projects/EndtoEnd/results/waterIPB/mono6MeV'
    boneResult='/data/qifan/projects/EndtoEnd/results/waterIPB/bone6MeV'
    shape = (256, 200)
    
    waterEdep = np.fromfile(os.path.join(waterResult, 'SD.bin'))
    boneEdep = np.fromfile(os.path.join(boneResult, 'SD.bin'))
    waterEdep = np.reshape(waterEdep, shape)
    boneEdep = np.reshape(boneEdep, shape)
    # we then normalize the dose with the maximum dose in water
    EdepMax = np.max(waterEdep)
    waterEdep /= EdepMax
    boneEdep /= EdepMax

    if False:
        # examine the centerline dose
        depth = np.arange(shape[0]) * resZ
        plt.plot(depth, waterEdep[:, 0])
        plt.plot(depth, boneEdep[:, 0])
        plt.xlabel('depth / cm')
        plt.ylabel('dose (a.u.)')
        plt.title('centerline dose (normalized by D_{max} in water)')
        figureFile = os.path.join(figureFolder, 'centerLine.png')
        plt.savefig(figureFile)
        plt.clf()
    
    if False:
        # Then we compare the lateral dose profile at different depths
        indices = [10, 20, 50, 100, 150, 250]  # the z index in water phantom
        rangeR = 25
        for idx in indices:
            # compute the z index in bone phantom
            idxBone = (idx + 0.5) * waterDensity / boneDensity
            idxBone = int(idxBone)

            sliceWater = waterEdep[idx, :rangeR]
            sliceBone = boneEdep[idxBone, :rangeR]
            sliceWater = sliceWater.copy()
            sliceBone = sliceBone.copy()
            sliceWater /= waterDensity ** 2
            sliceBone /= boneDensity ** 2

            radiusWater = (np.arange(rangeR) + 0.5) * resR * waterDensity
            radiusBone = (np.arange(rangeR) + 0.5) * resR * boneDensity
            plt.plot(radiusWater, sliceWater)
            plt.plot(radiusBone, sliceBone)
            plt.legend(['water', 'bone'])
            plt.xlabel('radius (g/cm^2)')
            plt.ylabel('energy deposition / density^2')
            plt.legend("water depth {} cm".format(idx*resZ))
            figureFile = os.path.join(figureFolder, 'depth{}.png'.format(idx))
            plt.savefig(figureFile)
            plt.clf()
    
    if True:
        # Then we compare the lateral dose profile at different depths
        indices = [10, 20, 50, 100, 150, 250]  # the z index in water phantom
        figIndices = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
        rangeR = 25
        fig, axs = plt.subplots(2, 3, figsize=(12, 8))
        for idx, figIdx in zip(indices, figIndices):
            # compute the z index in bone phantom
            idxBone = (idx + 0.5) * waterDensity / boneDensity
            idxBone = int(idxBone)

            sliceWater = waterEdep[idx, :rangeR]
            sliceBone = boneEdep[idxBone, :rangeR]
            sliceWater = sliceWater.copy()
            sliceBone = sliceBone.copy()
            sliceWater /= waterDensity ** 2
            sliceBone /= boneDensity ** 2

            radiusWater = (np.arange(rangeR) + 0.5) * resR * waterDensity
            radiusBone = (np.arange(rangeR) + 0.5) * resR * boneDensity
            axs[*figIdx].plot(radiusWater, sliceWater, color='g')
            axs[*figIdx].plot(radiusBone, sliceBone, color='b', linestyle='--')
            axs[*figIdx].legend(['water', 'bone'])
            axs[*figIdx].set_xlabel('radius (g/cm^3)')
            axs[*figIdx].set_ylabel('Edep / density^2')
            axs[*figIdx].set_title('depth {} g/cm^2'.format(idx * resZ))
        
        plt.tight_layout()
        figureFile = os.path.join(figureFolder, 'lateral.png')
        plt.savefig(figureFile)
        plt.clf()


def waterLateralDepth():
    """
    This function studies the dependency of lateral dose profile on depth
    """
    arrayPath = '/data/qifan/projects/EndtoEnd/results/waterIPB/mono6MeV/SD.bin'
    shape = (256, 200)
    array = np.fromfile(arrayPath)
    array = np.reshape(array, shape)
    figureFolder = './figures'

    resR = 0.05
    resZ = 0.1
    if False:
        # Firstly, we study the relationship between the centerline 
        # dose and the partial dose
        centerLine = array[:, 0].copy()
        centerLine /= np.max(centerLine)

        partial = np.sum(array, axis=1)
        partial /= np.max(partial)

        ratio = centerLine / partial

        depth = np.arange(shape[0]) * resZ

        fig, ax1 = plt.subplots()
        ax1.plot(depth, centerLine)
        ax1.plot(depth, partial)
        ax1.set_ylabel('Edep (a.u.)')
        ax1.legend(['centerline Edep', 'partialEdep'])

        ax2 = ax1.twinx()

        ax2.plot(depth, ratio, color='g')
        ax2.set_ylim([0, 1])
        ax2.set_ylabel('ratio')
        ax2.legend(['ratio'])
        plt.title('Centerline Edep vs partial Edep')
        figureFile = os.path.join(figureFolder, 'cenPar.png')
        plt.savefig(figureFile)
        plt.clf()

    if True:
        depths = [1, 5, 10, 20, 50, 100, 150, 200]
        rangeR = 25
        radius = (np.arange(rangeR) + 0.5) * resR
        for d in depths:
            slice = array[d, :rangeR]
            # normalize the slice against its maximum value
            slice /= np.max(slice)
            plt.plot(radius, slice)
        legend = ['depth {} cm'.format(d*resZ) for d in depths]
        plt.legend(legend)
        plt.xlabel('radius (cm)')
        plt.ylabel('lateral Edep (a.u.)')
        plt.title('lateral water dose profile at different depths')
        figureFile = os.path.join(figureFolder, 'latWater.png')
        plt.savefig(figureFile)
        plt.clf()


def inhomoExamine():
    """
    This function demonstrates the effect of inhomogeneity to lateral dose profile
    """
    figureFolder = './figures'
    arrayFile = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14/SD.bin'
    shape = (256, 200)
    array = np.fromfile(arrayFile, dtype=np.float64)
    array = np.reshape(array, shape)

    # Here we study the interface between the muscle and lung
    critical = 64
    Range = 6
    rangeR = 50
    resR = 0.05
    resZ = 0.1
    radius = (np.arange(rangeR) + 0.5) * resR

    if False:
        # the left figure, show the lateral dose profile at two sides of the interface
        for i in range(-Range+1, Range+1):
            layer = array[critical+i, :rangeR]
            plt.plot(radius, layer)
        legend = ['depth {:.1f} cm'.format((critical+i)*resZ) 
            for i in range(-Range+1, Range+1)]
        plt.legend(legend)
        plt.xlabel('radius (cm)')
        plt.ylabel('Edep')
        plt.title('dose around the muscle-lung interface')
        figureFile = os.path.join(figureFolder, 'muscleLung.png')
        plt.savefig(figureFile)

    if True:
        # muscle to lung
        fig, axs = plt.subplots(2, 2, figsize=(11, 8))
        Range = 6
        for i in range(-Range, Range):
            layer = array[critical + i, :rangeR]
            axs[0, 0].plot(radius, layer)
        legend = ['depth {:.1f} cm'.format((critical+i)*resZ) 
                  for i in range(-Range, Range)]
        axs[0, 0].legend(legend)
        axs[0, 0].set_xlabel('radius (cm)')
        axs[0, 0].set_ylabel('Edep (a.u.)')
        axs[0, 0].set_title('Edep near the muscle-lung interface')
        Range = 16
        for i in range(Range):
            layer = array[critical + i, :rangeR]
            axs[0, 1].plot(radius, layer)
        legend = ['depth {:.1f} cm'.format((critical+i)*resZ) 
                  for i in range(Range)]
        axs[0, 1].legend(legend)
        axs[0, 1].set_xlabel('radius (cm)')
        axs[0, 1].set_ylabel('Edep (a.u.)')
        axs[0, 1].set_title('Edep in the lung side of the muscle-lung interface')
        # plt.tight_layout()
        # figureFile = os.path.join(figureFolder, 'muscleLung.png')
        # plt.savefig(figureFile)

        # lung to muscle
        critical = 160
        Range = 6
        for i in range(-Range, Range):
            layer = array[critical + i, :rangeR]
            axs[1, 0].plot(radius, layer)
        legend = ['depth {:.1f} cm'.format((critical+i)*resZ) 
                  for i in range(-Range, Range)]
        axs[1, 0].legend(legend)
        axs[1, 0].set_xlabel('radius (cm)')
        axs[1, 0].set_ylabel('Edep (a.u.)')
        axs[1, 0].set_title('Edep near the lung-muscle interface')
        Range = 16
        for i in range(Range):
            layer = array[critical + i, :rangeR]
            axs[1, 1].plot(radius, layer)
        legend = ['depth {:.1f} cm'.format((critical+i)*resZ) 
                  for i in range(Range)]
        axs[1, 1].legend(legend)
        axs[1, 1].set_xlabel('radius (cm)')
        axs[1, 1].set_ylabel('Edep (a.u.)')
        axs[1, 1].set_title('Edep in the muscle side of the lung-muscle interface')
        plt.tight_layout()
        figureFile = os.path.join(figureFolder, 'lungMuscle.png')
        plt.savefig(figureFile)
        


if __name__ == '__main__':
    # readDose()
    # testCylindricalToCartesian()
    # TakeALook6MeV()
    # polyChromatic()
    # viewPoly()
    # viewCenterline()
    # homoDose()
    # waterLateralDepth()
    inhomoExamine()