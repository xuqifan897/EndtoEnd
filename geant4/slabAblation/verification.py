import os
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import math

KernelLong = None
KernelTrans = None
KernelPartial = None
phantomTheory = None
centerLine_pre = None
layers = [(16, 0.92), (16, 1.04), (16, 1.85), (16, 1.04), (96, 0.25),
    (16, 1.04), (16, 1.85), (16, 0.92), (16, 1.85), (16, 1.04), (16, 0.92)]


def kernelInit():
    """
    This function initializes both the longitudinal and the 
    transverse kernels
    """
    waterFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/waterDose'
    waterFile = os.path.join(waterFolder, 'array.bin')
    res = 0.1  # cm
    density = 1  # g/cm^3
    shape = (199, 199, 256)
    water = np.fromfile(waterFile, dtype=np.float64)
    water = np.reshape(water, shape)
    
    global KernelLong
    global KernelTrans
    global KernelPartial

    # initialize longitudinal kernel
    KernelLong = water[99, 99, :].copy()

    # initialize partial kernel
    KernelPartial = np.sum(water, axis=(0, 1))
    KernelPartial /= KernelPartial[100]
    
    # initialize lateral kernel, here we use the depth 10 cm kernel
    slice = water[:, :, 100]
    binScale = 10000
    distance, dose = slice2Tranverse(slice, binScale)
    # normalize dose with centerline dose
    dose = np.array(dose)
    distance = np.array(distance) / binScale * res
    KernelTrans = [distance, dose]


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


def phantomTheoryInit():
    """
    This funciton initializes the theoretical dose in phantom
    """
    # firstly, construct the phantom. Layers in (thickness, density)
    sizeZ = 0
    radioZ = 0.
    res = 0.1  # cm
    waterDensity = 1.0
    for a, b in layers:
        sizeZ += a
        radioZ += a * b

    # get the radiological depth
    radiologicalDepth = np.zeros(sizeZ)
    depth = -1
    rD = -1.
    for layer in layers:
        for depthLocal in range(layer[0]):
            depth += 1
            rD += layer[1]
            radiologicalDepth[depth] = rD
    
    # this part is to use centerline dose
    # interpolate the centerline from the radological depth
    global centerLine_pre
    centerLine_pre = np.zeros(sizeZ)
    for i in range(sizeZ):
        rD = radiologicalDepth[i]
        lower = int(np.floor(rD))
        higher = lower + 1
        lowerWeight = higher - rD
        higherWeight = rD - lower
        
        if False:
            # deal with the first element
            if lower < 0:
                centerLine_pre[0] = higherWeight * KernelLong[0]
            
            centerLine_pre[i] = lowerWeight * KernelLong[lower] + \
                higherWeight * KernelLong[higher]
        if True:
            if lower < 0:
                centerLine_pre[0] = higherWeight * KernelPartial[0]
            
            centerLine_pre[i] = lowerWeight * KernelPartial[lower] + \
                higherWeight * KernelPartial[higher]

    # then we obtain the full phantom dose by applying the transverse kernel
    shape = (199, 199, 256)
    radiologicalRadius = np.zeros((shape[0], shape[1]))
    for i in range(shape[0]):
        for j in range(shape[1]):
            radiologicalRadius[i, j] = math.sqrt((i-99)**2+(j-99)**2)
    
    global phantomTheory
    phantomTheory = np.zeros(shape)
    depthZ = -1
    for thickness, density in layers:
        # prepare the transverse kernel
        radiRadi = radiologicalRadius * res * density
        kernelLocal = np.interp(radiRadi, KernelTrans[0], KernelTrans[1])
        kernelLocal *= (density / waterDensity) ** 2

        for i in range(thickness):
            depthZ += 1
            phantomTheory[:, :, depthZ] = centerLine_pre[depthZ] * kernelLocal
    

def theoryExpComp():
    """
    This function compares the theoretical dose with the Monte Carlo phantom dose
    """
    phantomFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/phantomCorrect'
    resultFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly6/phantomDose'
    phantomFile = os.path.join(phantomFolder, 'array.bin')
    shape = (199, 199, 256)
    res = 0.1
    phantom = np.fromfile(phantomFile, dtype=np.float64)
    phantom = np.reshape(phantom, shape)

    # get the density map
    densityMap = np.zeros(shape[2])
    current = 0
    for thickness, density in layers:
        densityMap[current: current+thickness] = density
        current += thickness
    
    # calculate the dose to phantoms
    for i in range(shape[2]):
        phantom[:, :, i] /= densityMap[i]

    depth = np.arange(shape[2]) * res

    if True:
        # get the partial dose
        phantomPartial = np.sum(phantom, axis=(0, 1))
        theoryPartial = np.sum(phantomTheory, axis=(0, 1))
        plt.plot(depth, phantomPartial)
        plt.plot(depth, theoryPartial)
        plt.legend(['Monte Carlo', 'theory'])
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.title('partial dose comparison')
        file = os.path.join(resultFolder, 'slabPartial.png')
        plt.savefig(file)
        plt.clf()
    
    if True:
        # normalize
        phantomCenter = phantom[99, 99, :]
        theoryCenter = phantomTheory[99, 99, :]
        plt.plot(depth, phantomCenter)
        plt.plot(depth, theoryCenter)
        plt.legend(['Monte Carlo', 'theory'])
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.title('centerline dose comparison')
        file = os.path.join(resultFolder, 'slabCentral.png')
        plt.savefig(file)
        plt.clf()


if __name__ == '__main__':
    kernelInit()
    phantomTheoryInit()
    theoryExpComp()