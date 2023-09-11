import os
import numpy as np
import matplotlib.pyplot as plt

def readSmallMat():
    """
    For the first experiment, we found that Geant4 can now only work 
    with a small number of detectors. When the number of detectors 
    become too large, it just get stuck. Here we process and analyze 
    the data generated with phantom dimension (50, 50, 100)
    """
    resultFolder = "/data/qifan/projects/EndtoEnd/results/Sept1Point/initialDet"
    resultFile = os.path.join(resultFolder, 'SD.bin')
    dimension = (100, 50, 50)
    sourceDisplacement = 5
    resolution = 0.1  # cm
    array = np.fromfile(resultFile, dtype=np.float64)
    array = np.reshape(array, dimension)

    # view centerline dose
    centerLine = array[:, 25, 25]
    depth = (np.arange(dimension[0]) - sourceDisplacement) * resolution
    plt.plot(depth, centerLine)
    figureFile = os.path.join(resultFolder, 'centerline.png')
    plt.savefig(figureFile)
    plt.clf()

def readPointKernel():
    """
    We just generated an example point kernel file. Take a look
    """
    resultFolder = '/data/qifan/projects/EndtoEnd/results/Sept1Point/pointKernel'
    resultFile = os.path.join(resultFolder, 'SD.bin')
    dimension = (55, 50, 50)
    sourceD = 5
    resolution = 0.1  # cm
    array = np.fromfile(resultFile, dtype=np.double)
    # z, y, x
    array = np.reshape(array, dimension)

    depth = (np.arange(dimension[0]) - sourceD) * resolution

    # plot centerline profile
    centerLineIdxY = int(dimension[1] / 2)
    centerLineIdxX = int(dimension[2] / 2)
    centerLineKernel = array[:, centerLineIdxY, centerLineIdxX]
    # normalize
    centerLineKernel = centerLineKernel / np.max(centerLineKernel)
    plt.plot(depth, centerLineKernel)

    # plot partial profile
    partialKernel = np.sum(array, axis=(1, 2))
    # normalize
    partialKernel = partialKernel / np.max(partialKernel)
    plt.plot(depth, partialKernel)

    plt.legend(["Center line profile", "Partial profile"])
    plt.xlabel("depth (cm)")
    plt.ylabel("Energy deposition (a.u.)")
    figureFile = os.path.join(resultFolder, 'pointKernel.png')
    plt.savefig(figureFile)
    plt.clf()

if __name__ == '__main__':
    # readSmallMat()
    readPointKernel()