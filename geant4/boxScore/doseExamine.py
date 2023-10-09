import os
import numpy as np
import matplotlib.pyplot as plt

def simpleView():
    shape = (16, 256, 256)
    folder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly'
    data = os.path.join(folder, 'SD001.bin')
    array = np.fromfile(data, dtype=np.double)
    array = np.reshape(array, shape)

    centerLine = array[:, 128, 128]
    plt.plot(np.arange(shape[0]), centerLine)
    figureFile = os.path.join(folder, 'ddp001.png')
    plt.savefig(figureFile)
    plt.clf()


def doseConcat():
    """
    Previsouly, we successfully (probably) implemented the dose 
    calculation of the slab phantom with a diverging pencil beam.
    Now we concatenate the energy depostitions and calculate the dose.
    """
    folder = '/data/qifan/projects/EndtoEnd/results/slabBench/patient1_g1'
    shape = (256, 256, 256)
    datatype = np.double
    energyDeposition = np.zeros(shape, dtype=datatype)
    slabShape = (16, 256, 256)
    resolution = 0.1  # cm
    nSlabs = int(shape[0] / slabShape[0])

    for i in range(nSlabs):
        fileName = 'SD{:03d}.bin'.format(i+1)
        fileName = os.path.join(folder, fileName)
        slabArray = np.fromfile(fileName, dtype=datatype)
        slabArray = np.reshape(slabArray, slabShape)
        offset = i * slabShape[0]
        energyDeposition[offset: offset+slabShape[0], :, :] = slabArray
    
    depth = np.arange(shape[0]) * resolution  # cm

    if False:
        plt.plot(depth, energyDeposition[:, 128, 128])
        plt.xlabel('depth [cm]')
        plt.ylabel('Energy deposition')
        figureFile = os.path.join(folder, 'Edep.png')
        plt.savefig(figureFile)
        plt.clf()

    slabs = [("adipose", 0.8), ("muscle", 0.8), ("bone", 0.8), ("muscle", 0.8), 
             ("lung", 4.8), ("muscle", 0.8), ("bone", 0.8), ("adipose", 0.8),
             ("bone", 0.8), ("muscle", 0.8), ("adipose", 0.8)]
    LUTmat = {"water": 1., "adipose": 0.92, "muscle": 1.04, 
              "bone": 1.85, "lung": 0.25}
    layers = [(a, int(np.round(2*b/resolution)), LUTmat[a]) for a, b in slabs]
    
    localOffset = 0
    dose = energyDeposition.copy()
    for mat, lay, den in layers:
        dose[localOffset: localOffset+lay] /= den
        localOffset += lay
    
    plt.plot(depth, dose[:, 128, 128])
    plt.xlabel('depth (cm)')
    plt.ylabel('dose (a.u.)')
    figureFile = os.path.join(folder, 'dose.png')
    plt.savefig(figureFile)
    plt.clf()


def specVerify():
    """
    In this function, we simply do a sanity check of the spectrum
    """
    flunce = [0.001, 0.01, 0.02, 0.03, 0.068, 0.09, 0.101, 0.100, 0.131, 0.188, 0.140, 0.090, 0.030, 0.005]
    total = 0.
    for a in flunce:
        total += a
    print(total)
    

if __name__ == '__main__':
    simpleView()
    # doseConcat()
    # specVerify()