import os
import numpy as np
import matplotlib.pyplot as plt


def readDose():
    """
    This function reads the dose from the Monte Carlo simulation
    """
    resultFolder = '/data/qifan/projects/EndtoEnd4/results/InhomoJuly20/slab2'
    file1 = os.path.join(resultFolder, 'SD1.bin')
    file2 = os.path.join(resultFolder, 'SD2.bin')
    dim1 = 128
    dim2 = 128
    array1 = np.fromfile(file1, dtype=np.float64)
    array2 = np.fromfile(file2, dtype=np.float64)
    
    waterDensity = 1
    boneDensity = 1.85
    array1 /= waterDensity
    array2 /= boneDensity
    array = np.concatenate((array1, array2), axis=0)

    res = 0.1  # cm
    depth = np.arange(dim1 + dim2) * res
    plt.plot(depth, array)
    plt.legend(['water', 'bone'])
    plt.xlabel('depth (cm)')
    plt.ylabel('Edep / radiological depth')
    plt.title('a 2-layer slab phantom demo (water: 0-12.8cm, bone: 12.8-25.6cm)')
    
    outputFile = os.path.join(resultFolder, 'slab2.png')
    plt.savefig(outputFile)
    plt.clf()


if __name__ == '__main__':
    readDose()