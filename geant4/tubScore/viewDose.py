import os
import numpy as np
import matplotlib.pyplot as plt

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


def testCylindricalToCartisian():
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
    outputFile = '/data/qifan/projects/EndtoEnd/results/spec6MV/polyIPB.npy'
    polyDose = np.load(outputFile)

    # to take a look
    if False:
        centerLine = polyDose[:, 199:201, 199:201]
        centerLine = np.sum(centerLine, axis=(1, 2))
        depth = np.arange(256) * 0.1
        plt.plot(depth, centerLine)
        plt.show()
    
    nSlices = 256
    dimOrg = 400
    # window size: 0.5cm, 1cm, 2cm, corresponding to 10, 20, 40 pixels
    pixels = [10, 20, 40]
    for p in pixels:
        dimPad = dimOrg + p - 1
        outputShape = (nSlices, dimPad, dimPad)
        output = np.zeros(outputShape, dtype=np.float64)
        for i in range(nSlices):
            # output[i, :, :] = 


if __name__ == '__main__':
    # readDose()
    # testCylindricalToCartisian()
    # TakeALook6MeV()
    # polyChromatic()
    viewPoly()