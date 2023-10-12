import os
import numpy as np
import matplotlib.pyplot as plt

def viewCenterline():
    """
    This function takes a closer look at the centerline dose profile
    """
    expFolder = '/data/qifan/projects/EndtoEnd/results/Sept1Point/slabPoint'
    if not os.path.isdir(expFolder):
        os.mkdir(expFolder)
    
    kernelFile = '/data/qifan/projects/EndtoEnd/results/Sept1Point/pointKernel/SD.bin'
    datatype = np.double
    shape = (55, 50, 50)
    kernel = np.fromfile(kernelFile, dtype=datatype)
    kernel = np.reshape(kernel, shape)
    centerIdx = int(shape[1] / 2)
    centerLine  = kernel[:, centerIdx, centerIdx]
    indices = np.arange(shape[0])

    rangeStart = 0
    rangeEnd = 10
    indices = indices[rangeStart: rangeEnd]
    centerLine = centerLine[rangeStart: rangeEnd]
    plt.plot(indices, centerLine)
    file = os.path.join(expFolder, 'kernelDDP.png')
    plt.savefig(file)


if __name__ == '__main__':
    viewCenterline()