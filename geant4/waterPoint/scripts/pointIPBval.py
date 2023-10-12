import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def point2Cartesian():
    """
    This function synthesizes the water dose using the point kernel we obtained,
    and compare it against the IPB dose kernel obtained through the Monte Carlo 
    calculation
    """
    kernelFolder = '/data/qifan/projects/EndtoEnd/results/Sept1Point/pointKernel'
    kernelEvents = 1171944
    kernelDim = (55, 50, 50)
    PhantomSZ = 5
    PhantomBottom = 50
    datatype = np.double
    kernelFile = os.path.join(kernelFolder, 'SD.bin')
    kernel = np.fromfile(kernelFile, dtype=datatype)
    kernel = np.reshape(kernel, kernelDim)

    # For 6 MeV photon
    MuOverRho = 2.77e-2  # cm^2/g
    resolution = 0.1
    totalDepth = 256
    depth = (np.arange(totalDepth) + 0.5) * resolution
    Terma = np.exp(- depth * MuOverRho)

    if False:
        plt.plot(depth, Terma)
        file = 'Terma.png'
        plt.savefig(file)
        plt.clf()
    
    IPBDepth = totalDepth + kernelDim[0] - 1
    IPB = np.zeros((IPBDepth, kernelDim[1], kernelDim[2]), dtype=datatype)
    for i in range(kernelDim[1]):
        for j in range(kernelDim[2]):
            IPB[:, i, j] = np.convolve(Terma, kernel[:, i, j])
    # trim
    IPB = IPB[PhantomSZ-1: -PhantomBottom, :, :]
    
    expFolder = '/data/qifan/projects/EndtoEnd/results/waterIPB/PointEval'
    arrayFile = os.path.join(expFolder, 'Point2Cartesian.npy')
    np.save(arrayFile, IPB)


def readIPB():
    """
    This function reads the IPB dose. In the calculation, a Cartesian phantom was used.
    """
    IPBFolder = '/data/qifan/projects/EndtoEnd/results/waterIPB/mono6MeV'
    IPBFile = os.path.join(IPBFolder, 'SD.bin')
    shape = (256, 200)
    datatype = np.double
    IPB = np.fromfile(IPBFile, dtype=datatype)
    IPB = np.reshape(IPB, shape)

    IPBCart = Cyd2Cart(IPB)

    expFolder = '/data/qifan/projects/EndtoEnd/results/waterIPB/PointEval'
    if not os.path.isdir(expFolder):
        os.mkdir(expFolder)
    arrayFile = os.path.join(expFolder, 'Cylinder2Cartesian.npy')
    np.save(arrayFile, IPBCart)


def Cyd2Cart(input):
    """
    This function converts a cylindrical phantom into a Cartesian phantom
    """
    inputShape = input.shape
    outputShape = (inputShape[0], inputShape[1], inputShape[1])
    output = np.zeros(outputShape, dtype=input.dtype)

    # Then, we prepare the indece for interpolation
    resolution = 0.1  # unit: cm
    middle = (inputShape[1] - 1) / 2
    idxX = np.arange(inputShape[1]) - middle
    idxY = idxX.copy()
    idxX = np.expand_dims(idxX, axis=1)
    idxY = np.expand_dims(idxY, axis=0)
    dist = np.sqrt(np.square(idxX) + np.square(idxY)) * resolution
    dist = dist.flatten()

    if False:
        figureFile = './dist.png'
        plt.imsave(figureFile, dist)
        plt.clf()
    
    # convert cylindrical value to Cartesian value. Because 
    # cylindrical energy is cummulative, we have to devide 
    # it by the volume it occupies
    area = np.zeros(inputShape[1])
    for i in range(inputShape[1]):
        area[i] = (i+1)**2 - i**2
    area = np.expand_dims(area, axis=0)
    input_convert = input / area

    paddingSize = 400
    inputSlicePadding = np.zeros(paddingSize, dtype=input_convert.dtype)
    xTicks = (np.arange(paddingSize) + 0.6) * resolution / 2
    for i in range(outputShape[0]):
        inputSlice = input_convert[i, :]
        inputSlice = inputSlice.copy()
        inputSlicePadding[:input_convert.shape[1]] = inputSlice
        f = interpolate.interp1d(xTicks, inputSlicePadding)
        outputSlice = f(dist)
        outputSlice = np.reshape(outputSlice, (outputShape[1], outputShape[2]))
        output[i, :, :] = outputSlice
    return output


def compare():
    """
    This function compares the two IPB dose from the two methods,
    tub phantom and point convolution
    """
    expFolder = '/data/qifan/projects/EndtoEnd/results/waterIPB/PointEval'
    resolution = 0.1  # cm
    pointFile = os.path.join(expFolder, 'Point2Cartesian.npy')
    IPBPoint = np.load(pointFile)
    IPBPoint /= np.max(IPBPoint)

    cylinderFile = os.path.join(expFolder, 'Cylinder2Cartesian.npy')
    IPBCylinder = np.load(cylinderFile)
    IPBCylinder /= np.max(IPBCylinder)

    if True:
        # Plot centerline dose
        centerIdx = int(IPBPoint.shape[1] / 2)
        centerlinePoint = IPBPoint[:, centerIdx, centerIdx]
        depth = (np.arange(IPBPoint.shape[0]) + 0.5) * resolution
        plt.plot(depth, centerlinePoint)

        centerIdx = int(IPBCylinder.shape[1] / 2)
        centerlineCylinder = IPBCylinder[:, centerIdx, centerIdx]
        plt.plot(depth, centerlineCylinder)
        plt.legend(['Terma kernel convolution', 'Direct MC'])
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.title('Point kernel validation in water')

        file = os.path.join(expFolder, 'IPBvsConvLongitudinal.png')
        plt.savefig(file)
        plt.clf()
    
    if True:
        # Plot transverse profile
        depths = [20, 50, 100, 150, 200]  # unitless
        width = 16
        pointWidth = IPBPoint.shape[1]
        pointStart = int((pointWidth - width) / 2)
        pointEnd = int((pointWidth + width) / 2)
        pointMiddle = int(pointWidth / 2)
        cylinderWidth = IPBCylinder.shape[1]
        cylinderStart = int((cylinderWidth - width) / 2)
        cylinderEnd = int((cylinderWidth + width) / 2)
        cylinderMiddle = int(cylinderWidth / 2)
        axis = np.arange(width) - (width-1) / 2
        axis *= resolution

        for dep in depths:
            pointLine = IPBPoint[dep, pointMiddle, pointStart: pointEnd]
            cylinderLine = IPBCylinder[dep, cylinderMiddle, cylinderStart: cylinderEnd]
            line, = plt.plot(axis, pointLine)
            color = 'red'
            plt.plot(axis, cylinderLine, color, linestyle='--')
            plt.legend(['Direct MC', 'Terma kernel convolution'])
            plt.xlabel('lateral distance (cm)')
            plt.ylabel('dose (a.u.)')
            plt.ylim(0., 1.2)
            plt.title('lateral dose profile at depth {}cm'.format(dep * resolution))
            file = os.path.join(expFolder, 'IPBvsConvDepth{}.png'.format(dep * resolution))
            plt.savefig(file)
            plt.clf()


if __name__ == '__main__':
    # point2Cartesian()
    # readIPB()
    compare()