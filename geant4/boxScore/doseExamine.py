import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import signal

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
    folder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly'
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


def MC_CCCS_comp():
    """
    This function compares the CCCS and Geant4 based Monte Carlo dose calculation
    """
    CCCS_folder = '/data/qifan/projects/EndtoEnd/results/' \
        'slabBench/slab_dosecalc_9_0.5'
    MC_folder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_0.8_200'
    MCShape = (256, 64, 64)
    # MCShape = (256, 256, 256)

    CCCSShape = (256, 256, 256)
    size = CCCSShape[0] * CCCSShape[1] * CCCSShape[2]
    datatype = np.double
    FluenceDim = 9
    resolution = 0.1  # cm

    CCCSDoseFile = os.path.join(CCCS_folder, 'Dose_Coefficients.h5')
    beamletIdx = int((FluenceDim - 1) / 2 * FluenceDim + (FluenceDim - 1) / 2)
    dataset = h5py.File(CCCSDoseFile, 'r')
    dataset = dataset['beams']['data']
    beam = dataset['beam_00000']
    BeamletKey = 'beamlet_{:05d}'.format(beamletIdx)
    beamlet = beam[BeamletKey]
    coeffs = beamlet['coeffs'][()]
    lindex = beamlet['lindex'][()]

    CCCSDose = np.zeros(size, dtype=datatype)
    for key, value in zip(lindex, coeffs):
        CCCSDose[key] = value
    CCCSDose = np.reshape(CCCSDose, CCCSShape)

    # find the location for the beamlet by calculating the mass center
    centralIdx = int(CCCSShape[1] / 2)
    slice = CCCSDose[:, centralIdx, :]
    result1 = 0.
    for i in range(slice.shape[0]):
        result1 += i * np.sum(slice[i, :])
    result1 /= np.sum(slice)
    result2 = 0.
    for i in range(slice.shape[1]):
        result2 += i * np.sum(slice[:, i])
    result2 /= np.sum(slice)
    result1, result2 = int(result1), int(result2)
    CCCSDoseLine = CCCSDose[result1, :, result2]

    depth = np.arange(CCCSShape[1]) * resolution
    if False:
        plt.plot(depth, CCCSDoseLine)
        figureFile = os.path.join(MC_folder, 'CCCSDose.png')
        plt.savefig(figureFile)
        plt.clf()

    # The other part
    energyDeposition = np.zeros(MCShape, dtype=datatype)
    slabShape = (16, MCShape[1], MCShape[1])
    resolution = 0.1  # cm
    nSlabs = int(MCShape[0] / slabShape[0])
    for i in range(nSlabs):
        fileName = 'SD{:03d}.bin'.format(i+1)
        fileName = os.path.join(MC_folder, fileName)
        slabArray = np.fromfile(fileName, dtype=datatype)
        slabArray = np.reshape(slabArray, slabShape)
        offset = i * slabShape[0]
        energyDeposition[offset: offset+slabShape[0], :, :] = slabArray
    slabs = [("adipose", 0.8), ("muscle", 0.8), ("bone", 0.8), ("muscle", 0.8), 
             ("lung", 4.8), ("muscle", 0.8), ("bone", 0.8), ("adipose", 0.8),
             ("bone", 0.8), ("muscle", 0.8), ("adipose", 0.8)]
    LUTmat = {"water": 1., "adipose": 0.92, "muscle": 1.04, 
              "bone": 1.85, "lung": 0.25}
    layers = [(a, int(np.round(2*b/resolution)), LUTmat[a]) for a, b in slabs]
    
    localOffset = 0
    MCDose = energyDeposition.copy()
    for mat, lay, den in layers:
        MCDose[localOffset: localOffset+lay] /= den
        localOffset += lay
    
    MCCenterLineIdx = int(MCShape[1] / 2)
    MCDoseLine = MCDose[:, MCCenterLineIdx, MCCenterLineIdx]

    # normalize to minimize the L2 distance
    if True:
        CCCSDoseLine = CCCSDoseLine.copy()
        CCCSDoseLine /= np.max(CCCSDoseLine)
        MCDoseLine = MCDoseLine.copy()
        sum_MC_square = np.sum(np.square(MCDoseLine))
        sum_CCCS_MC = np.sum(CCCSDoseLine * MCDoseLine)
        scale = sum_CCCS_MC / sum_MC_square
        MCDoseLine *= scale

        plt.plot(depth, CCCSDoseLine)
        plt.plot(depth, MCDoseLine)
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        figureFile = os.path.join(MC_folder, 'DoseComp.png')
        plt.savefig(figureFile)
        plt.clf()
    

def IPBConvExam():
    """
    According to the results associated with the code above, we observed that
    the higher SAD is, the better the Monte Carlo results matches the CCCS result.
    So we hypothesize that the CCCS dose could be calculated using a non-diverging 
    Terma distribution
    """
    IPBFolder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_0.0_100'
    MCShape = (256, 63, 63)
    datatype = np.double

    energyDeposition = np.zeros(MCShape, dtype=datatype)
    slabShape = (16, MCShape[1], MCShape[1])
    resolution = 0.1  # cm
    nSlabs = int(MCShape[0] / slabShape[0])
    for i in range(nSlabs):
        fileName = 'SD{:03d}.bin'.format(i+1)
        fileName = os.path.join(IPBFolder, fileName)
        slabArray = np.fromfile(fileName, dtype=datatype)
        slabArray = np.reshape(slabArray, slabShape)
        offset = i * slabShape[0]
        energyDeposition[offset: offset+slabShape[0], :, :] = slabArray
    slabs = [("adipose", 0.8), ("muscle", 0.8), ("bone", 0.8), ("muscle", 0.8), 
             ("lung", 4.8), ("muscle", 0.8), ("bone", 0.8), ("adipose", 0.8),
             ("bone", 0.8), ("muscle", 0.8), ("adipose", 0.8)]
    LUTmat = {"water": 1., "adipose": 0.92, "muscle": 1.04, 
              "bone": 1.85, "lung": 0.25}
    layers = [(a, int(np.round(2*b/resolution)), LUTmat[a]) for a, b in slabs]
    
    localOffset = 0
    IPBDose = energyDeposition.copy()
    for mat, lay, den in layers:
        IPBDose[localOffset: localOffset+lay] /= den
        localOffset += lay
    
    convWindowSize = 10
    convKernel = np.ones((convWindowSize, convWindowSize), dtype=datatype)
    beamDose = np.zeros_like(IPBDose)
    for i in range(MCShape[0]):
        beamDose[i, :, :] = signal.convolve2d(IPBDose[i, :, :], convKernel, 'same')
    

    # prepare for CCCS dose
    CCCS_folder = '/data/qifan/projects/EndtoEnd/results/' \
        'slabBench/slab_dosecalc_9_0.5'

    CCCSShape = (256, 256, 256)
    size = CCCSShape[0] * CCCSShape[1] * CCCSShape[2]
    datatype = np.double
    FluenceDim = 9
    resolution = 0.1  # cm

    CCCSDoseFile = os.path.join(CCCS_folder, 'Dose_Coefficients.h5')
    beamletIdx = int((FluenceDim - 1) / 2 * FluenceDim + (FluenceDim - 1) / 2)
    dataset = h5py.File(CCCSDoseFile, 'r')
    dataset = dataset['beams']['data']
    beam = dataset['beam_00000']
    BeamletKey = 'beamlet_{:05d}'.format(beamletIdx)
    beamlet = beam[BeamletKey]
    coeffs = beamlet['coeffs'][()]
    lindex = beamlet['lindex'][()]

    CCCSDose = np.zeros(size, dtype=datatype)
    for key, value in zip(lindex, coeffs):
        CCCSDose[key] = value
    CCCSDose = np.reshape(CCCSDose, CCCSShape)


    # plot
    # find the location for the beamlet by calculating the mass center
    centralIdx = int(CCCSShape[1] / 2)
    slice = CCCSDose[:, centralIdx, :]
    result1 = 0.
    for i in range(slice.shape[0]):
        result1 += i * np.sum(slice[i, :])
    result1 /= np.sum(slice)
    result2 = 0.
    for i in range(slice.shape[1]):
        result2 += i * np.sum(slice[:, i])
    result2 /= np.sum(slice)
    result1, result2 = int(result1), int(result2)
    CCCSDoseLine = CCCSDose[result1, :, result2]
    CCCSDoseLine = CCCSDoseLine.copy()
    CCCSDoseLine /= np.max(CCCSDoseLine)

    beamDoseWidth = beamDose.shape[1]
    idx = int((beamDoseWidth-1)/2)
    beamDoseLine = beamDose[:, idx, idx]
    beamDoseLine = beamDoseLine.copy()

    # # normalize by its maximum value
    # beamDoseLine /= np.max(beamDoseLine)

    # take into account the inverse square law
    SAD = 200  # cm
    r0 = SAD - 12.8  # cm
    scale = np.arange(MCShape[0]) * resolution + r0
    scale = SAD**2 / scale**2
    beamDoseLine *= scale

    # # normalize it w.r.t. the maximum value
    # beamDoseLine /= np.max(beamDoseLine)

    # normalize it to minimize the L2 distance
    sum_beamDoseLine_square = np.sum(np.square(beamDoseLine))
    sum_beamDoseLine_CCCSDoseLine = np.sum(beamDoseLine * CCCSDoseLine)
    scale = sum_beamDoseLine_CCCSDoseLine / sum_beamDoseLine_square
    beamDoseLine *= scale

    depth = (np.arange(MCShape[0]) + 0.5) * resolution
    plt.plot(depth, CCCSDoseLine)
    plt.plot(depth, beamDoseLine)
    # plt.plot(depth, scale)
    figureFile = os.path.join(IPBFolder, 'IPBconvInverseSquare.png')
    plt.savefig(figureFile)
    plt.clf()


def IPBG4Exam():
    """
    According to the results associated with the code above, we observed that
    the higher SAD is, the better the Monte Carlo results matches the CCCS result.
    So we hypothesize that the CCCS dose could be calculated using a non-diverging 
    Terma distribution
    """
    IPBFolder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_1.0_250_1e8'
    MCShape = (256, 63, 63)
    datatype = np.double

    energyDeposition = np.zeros(MCShape, dtype=datatype)
    slabShape = (16, MCShape[1], MCShape[1])
    resolution = 0.1  # cm
    nSlabs = int(MCShape[0] / slabShape[0])
    for i in range(nSlabs):
        fileName = 'SD{:03d}.bin'.format(i+1)
        fileName = os.path.join(IPBFolder, fileName)
        slabArray = np.fromfile(fileName, dtype=datatype)
        slabArray = np.reshape(slabArray, slabShape)
        offset = i * slabShape[0]
        energyDeposition[offset: offset+slabShape[0], :, :] = slabArray
    slabs = [("adipose", 0.8), ("muscle", 0.8), ("bone", 0.8), ("muscle", 0.8), 
             ("lung", 4.8), ("muscle", 0.8), ("bone", 0.8), ("adipose", 0.8),
             ("bone", 0.8), ("muscle", 0.8), ("adipose", 0.8)]
    LUTmat = {"water": 1., "adipose": 0.92, "muscle": 1.04, 
              "bone": 1.85, "lung": 0.25}
    layers = [(a, int(np.round(2*b/resolution)), LUTmat[a]) for a, b in slabs]
    
    localOffset = 0
    IPBDose = energyDeposition.copy()
    for mat, lay, den in layers:
        IPBDose[localOffset: localOffset+lay] /= den
        localOffset += lay

    # prepare for CCCS dose
    CCCS_folder = '/data/qifan/projects/EndtoEnd/results/' \
        'slabBench/slab_dosecalc_9_0.5'

    CCCSShape = (256, 256, 256)
    size = CCCSShape[0] * CCCSShape[1] * CCCSShape[2]
    datatype = np.double
    FluenceDim = 9
    resolution = 0.1  # cm

    CCCSDoseFile = os.path.join(CCCS_folder, 'Dose_Coefficients.h5')
    beamletIdx = int((FluenceDim - 1) / 2 * FluenceDim + (FluenceDim - 1) / 2)
    dataset = h5py.File(CCCSDoseFile, 'r')
    dataset = dataset['beams']['data']
    beam = dataset['beam_00000']
    BeamletKey = 'beamlet_{:05d}'.format(beamletIdx)
    beamlet = beam[BeamletKey]
    coeffs = beamlet['coeffs'][()]
    lindex = beamlet['lindex'][()]

    CCCSDose = np.zeros(size, dtype=datatype)
    for key, value in zip(lindex, coeffs):
        CCCSDose[key] = value
    CCCSDose = np.reshape(CCCSDose, CCCSShape)


    # plot
    # find the location for the beamlet by calculating the mass center
    centralIdx = int(CCCSShape[1] / 2)
    slice = CCCSDose[:, centralIdx, :]
    result1 = 0.
    for i in range(slice.shape[0]):
        result1 += i * np.sum(slice[i, :])
    result1 /= np.sum(slice)
    result2 = 0.
    for i in range(slice.shape[1]):
        result2 += i * np.sum(slice[:, i])
    result2 /= np.sum(slice)
    result1, result2 = int(result1), int(result2)
    CCCSDoseLine = CCCSDose[result1, :, result2]
    CCCSDoseLine = CCCSDoseLine.copy()
    CCCSDoseLine /= np.max(CCCSDoseLine)

    IPBDoseWidth = IPBDose.shape[1]
    idx = int((IPBDoseWidth-1)/2)
    IPBDoseLine = IPBDose[:, idx, idx]
    IPBDoseLine = IPBDoseLine.copy()
    # # normalize w.r.t the maximum value
    # IPBDoseLine /= np.max(IPBDoseLine)

    # normalize to minimize the mse
    sum_IPBDoseLine_square = np.sum(np.square(IPBDoseLine))
    sum_IPBDoseLine_CCCSDoseLine = np.sum(IPBDoseLine * CCCSDoseLine)
    scale = sum_IPBDoseLine_CCCSDoseLine / sum_IPBDoseLine_square
    IPBDoseLine *= scale

    depth = (np.arange(MCShape[0]) + 0.5) * resolution
    plt.plot(depth, CCCSDoseLine)
    plt.plot(depth, IPBDoseLine)
    figureFile = os.path.join(IPBFolder, 'CCCSMCLong.png')
    plt.savefig(figureFile)
    plt.clf()


if __name__ == '__main__':
    # simpleView()
    # doseConcat()
    # specVerify()
    # MC_CCCS_comp()
    IPBConvExam()
    # IPBG4Exam()