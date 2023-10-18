import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy import signal

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
    
    convWindowSize = 15
    convKernel = np.ones((convWindowSize, convWindowSize), dtype=datatype)
    beamDose = np.zeros_like(IPBDose)
    for i in range(MCShape[0]):
        beamDose[i, :, :] = signal.convolve2d(IPBDose[i, :, :], convKernel, 'same')
    

    # prepare for CCCS dose
    CCCS_folder = '/data/qifan/projects/EndtoEnd/results/' \
        'slabBench/slab_dosecalc_9_2.0_0.25'

    CCCSShape = (103, 103, 103)
    size = CCCSShape[0] * CCCSShape[1] * CCCSShape[2]
    datatype = np.double
    FluenceDim = 9
    resolutionCCCS = 0.25  # cm

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

    # take into account the inverse square law
    SAD = 200  # cm
    r0 = SAD - 12.8  # cm
    scale = np.arange(MCShape[0]) * resolution + r0
    scale = SAD**2 / scale**2
    beamDoseLine *= scale

    if False:
        # normalize it against the maximum value
        beamDoseLine /= np.max(beamDoseLine)
        depthCCCS = (np.arange(CCCSShape[1]) + 0.5 ) * resolutionCCCS
        depth = (np.arange(MCShape[0]) + 0.5) * resolution
        plt.plot(depthCCCS, CCCSDoseLine)
        plt.plot(depth, beamDoseLine)
        figureFile = os.path.join(IPBFolder, 'IPBconvInverseSquareBeamlet2SAD100.png')
        plt.savefig(figureFile)

    if True:
        # interpolate the CCCSDoseLine into the shape of beamDoseLine
        depthCCCS = (np.arange(CCCSShape[1]) + 0.5 ) * resolutionCCCS
        depth = (np.arange(MCShape[0]) + 0.5) * resolution
        CCCSDoseLineInterp = np.interp(depth, depthCCCS, CCCSDoseLine)
        # normalize to minimize the L2 distance
        sum_beamDoseLine_square = np.sum(np.square(beamDoseLine))
        sum_beamDoseLine_CCCSDoseLine = np.sum(beamDoseLine * CCCSDoseLineInterp)
        scale = sum_beamDoseLine_CCCSDoseLine / sum_beamDoseLine_square
        beamDoseLine *= scale
        plt.plot(depth, CCCSDoseLineInterp)
        plt.plot(depth, beamDoseLine)
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.legend(['CCCS', 'Monte Carlo'])
        plt.title('depth dose beamlet 2.0 cm')
        figureFile = os.path.join(IPBFolder, 'IPBconvInverseSquareBeamlet2SAD100.png')
        plt.savefig(figureFile)
        plt.clf()


def getConvDose(windowSize, IPBFolder, MCShape):
    """
    This function gets the convolution of the IPB dose
    """
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
    
    convKernel = np.ones((windowSize, windowSize), dtype=datatype)
    beamDose = np.zeros_like(IPBDose)
    for i in range(MCShape[0]):
        beamDose[i, :, :] = signal.convolve2d(IPBDose[i, :, :], convKernel, 'same')
    return beamDose


def getCCCSDose(file, shape):
    """
    This function gets the CCCS dose
    """
    size = shape[0] * shape[1] * shape[2]
    datatype = np.double
    FluenceDim = 9

    beamletIdx = int((FluenceDim - 1) / 2 * FluenceDim + (FluenceDim - 1) / 2)
    dataset = h5py.File(file, 'r')
    dataset = dataset['beams']['data']
    beam = dataset['beam_00000']
    BeamletKey = 'beamlet_{:05d}'.format(beamletIdx)
    beamlet = beam[BeamletKey]
    coeffs = beamlet['coeffs'][()]
    lindex = beamlet['lindex'][()]

    CCCSDose = np.zeros(size, dtype=datatype)
    for key, value in zip(lindex, coeffs):
        CCCSDose[key] = value
    CCCSDose = np.reshape(CCCSDose, shape)

    centralIdx = int(shape[1] / 2)
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
    return CCCSDose, result1, result2


def lateralProfileExamineConv():
    """
    This function examines the lateral dose profile of the both 
    the CCCS dose and the diverging dose obtained by covolving 
    the IPB dose with the fluence map, then undergoes an inverse 
    square transformation
    """
    CCCS_folder = '/data/qifan/projects/EndtoEnd/results/' \
        'slabBench/slab_dosecalc_9_0.5'
    CCCS_file = os.path.join(CCCS_folder, 'Dose_Coefficients.h5')
    shape = (256, 256, 256)
    CCCSDose, coord1, coord2 = getCCCSDose(CCCS_file, shape)
    # normalize
    CCCSDose /= np.max(CCCSDose)
    CCCSCenterLine = CCCSDose[coord1, :, coord2]

    IPBFolder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_0.0_100'
    MCShape = (256, 63, 63)
    IPBResolution = 0.1
    convWindowSize = 10
    parallelBeam = getConvDose(convWindowSize, IPBFolder, MCShape)
    centerIdx = int(MCShape[1] / 2)
    parallelCenterline = parallelBeam[:, centerIdx, centerIdx]

    SAD = 200
    r0 = SAD - 12.8
    depth = r0 + (np.arange(MCShape[0]) + 0.5) * IPBResolution
    scale = SAD**2 / depth**2
    parallelCenterline = parallelCenterline * scale
    # then we normalize to minimize the L2 distance
    sum_parallelCenterline_square = np.sum(np.square(parallelCenterline))
    sum_parallelCenterline_CCCSCenterLine = np.sum(parallelCenterline*CCCSCenterLine)
    factor = sum_parallelCenterline_CCCSCenterLine / sum_parallelCenterline_square

    if False:
        parallelCenterline *= factor
        plt.plot(depth, CCCSCenterLine)
        plt.plot(depth, parallelCenterline)
        file = os.path.join(IPBFolder, 'lateralSAD200Window10.png')
        plt.savefig(file)
        plt.clf()
    
    parallelBeam *= factor

    # then we examine the lateral profile
    depthCCCS = 100
    depthParallel = 100
    CCCSLateral = CCCSDose[coord1, depthCCCS, :]
    parallelLateral = parallelBeam[depthParallel, centerIdx, :]
    r = r0 + depthParallel*IPBResolution
    scale_ = (SAD / r)**2
    parallelLateral = parallelLateral * scale_

    CCCSAxis = (np.arange(shape[2]) - coord2 - 0.5) * IPBResolution * 2
    parallelAxis = (np.arange(MCShape[2]) - (MCShape[2]-1) / 2) * IPBResolution
    # scale the parallelAxis
    parallelAxis *= r / SAD

    intervalCCCS = 20
    intervalParallel = 40
    # CCCSPre = int((shape[2] - intervalCCCS) / 2)
    # CCCSPost = int((shape[2] + intervalCCCS) / 2)
    CCCSPre = int(coord2 - intervalCCCS / 2 + 1)
    CCCSPost = int(coord2 + intervalCCCS / 2 + 1)
    parallelPre = int((MCShape[2] - intervalParallel) / 2)
    parallelPost = int((MCShape[2] + intervalParallel) / 2)
    plt.plot(CCCSAxis[CCCSPre: CCCSPost], CCCSLateral[CCCSPre: CCCSPost])
    plt.plot(parallelAxis[parallelPre: parallelPost], parallelLateral[parallelPre: parallelPost])
    plt.xlabel('lateral distance (cm)')
    plt.ylabel('dose (a.u.)')
    plt.legend(['CCCS', 'Monte Carlo'])
    plt.title('lateral beamlet 1.0 cm')
    file = os.path.join(IPBFolder, 'lateralSAD200Window10.png')
    plt.savefig(file)
    plt.clf()


def getMCDose(MC_folder):
    """
    This function gets the Monte Carlo dose
    """
    MCShape = (256, 63, 63)
    datatype = np.double
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
    return MCDose


def beamlet2_0Longitudinal():
    """
    This function processes the dose comparison of beamlet size 2 cm
    """
    CCCSFolder = '/data/qifan/projects/EndtoEnd/results/slabBench/slab_dosecalc_9_2.0_0.25'
    CCCSFile = os.path.join(CCCSFolder, 'Dose_Coefficients.h5')
    CCCSShape = (103, 103, 103)
    CCCSResolution = 0.25
    CCCSDose, coord1, coord2 = getCCCSDose(CCCSFile, CCCSShape)
    CCCSDose /= np.max(CCCSDose)
    CCCSCenterline = CCCSDose[coord1, :, coord2].copy()
    CCCSDepth = (np.arange(CCCSShape[0]) + 0.5) * CCCSResolution

    if False:
        plt.plot(CCCSDepth, CCCSCenterline)
        figureFile = os.path.join(CCCSFolder, 'longitudinal.png')
        plt.savefig(figureFile)
    
    # prepare Monte Carlo dose
    MCFolder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_1.5_200_1e8'
    MCShape = (256, 63, 63)
    MCResolution = 0.1
    MCDose = getMCDose(MCFolder)
    MCCenterIdx = int(MCShape[1] / 2)
    MCCenterline = MCDose[:, MCCenterIdx, MCCenterIdx].copy()
    MCDepth = (np.arange(MCShape[0]) + 0.5) * MCResolution
    
    # firstly, compare the centerline dose
    CCCSCenterlineInterpolated = np.interp(MCDepth, CCCSDepth, CCCSCenterline)
    sum_MCCenterline_Square = np.sum(np.square(MCCenterline))
    sum_MCCenterline_CCCSCenterline = np.sum(MCCenterline * CCCSCenterlineInterpolated)
    factor = sum_MCCenterline_CCCSCenterline / sum_MCCenterline_Square

    if True:
        MCCenterline *= factor
        plt.plot(MCDepth, CCCSCenterlineInterpolated)
        plt.plot(MCDepth, MCCenterline)
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.title('depth dose beamlet 2.0 cm')
        plt.legend(['CCCS', 'Monte Carlo'])
        figureFile = os.path.join(MCFolder, 'CCCSMCLong.png')
        plt.savefig(figureFile)
        plt.clf()
    
    MCDose *= factor
    depth = 10.0  # cm
    CCCSDepthIdx = int(depth / CCCSResolution)
    CCCSLateral = CCCSDose[coord1, CCCSDepthIdx, :]
    CCCSAxis = (np.arange(CCCSShape[2]) - coord2 - 1) * CCCSResolution
    CCCSinterval = 16
    CCCSbegin = int(coord2 + 1 - CCCSinterval / 2)
    CCCSend = int(coord2 + 1 + CCCSinterval / 2)

    MCDepthIdx = int(depth / MCResolution)
    MCLateral = MCDose[MCDepthIdx, MCCenterIdx, :]
    MCAxis = np.arange(MCShape[2]) - (MCShape[2]-1) / 2
    MCAxis *= MCResolution * 20 / 15
    MCinterval = 30
    MCbegin = int((MCShape[2] - 1 - MCinterval) / 2)
    MCend = int((MCShape[2] - 1 + MCinterval) / 2)
    plt.plot(CCCSAxis[CCCSbegin:CCCSend], CCCSLateral[CCCSbegin:CCCSend])
    plt.plot(MCAxis[MCbegin:MCend], MCLateral[MCbegin:MCend])
    plt.legend(['CCCS', 'Monte Carlo'])
    plt.xlabel('lateral distance (cm)')
    plt.ylabel('dose (a.u.)')
    plt.title('lateral beamlet 2.0 cm')
    figureFile = os.path.join(MCFolder, 'lateralSAD200window15.png')
    plt.savefig(figureFile)


def lateralProfileExamineConv_beamlet2():
    """
    This function examines the lateral dose profile of the both 
    the CCCS dose and the diverging dose obtained by covolving 
    the IPB dose with the fluence map, then undergoes an inverse 
    square transformation
    """
    CCCS_folder = '/data/qifan/projects/EndtoEnd/results/' \
        'slabBench/slab_dosecalc_9_2.0_0.25'
    CCCS_file = os.path.join(CCCS_folder, 'Dose_Coefficients.h5')
    shape = (103, 103, 103)
    CCCSResolution = 0.25
    CCCSDose, coord1, coord2 = getCCCSDose(CCCS_file, shape)
    # normalize
    CCCSDose /= np.max(CCCSDose)
    CCCSCenterLine = CCCSDose[coord1, :, coord2]

    IPBFolder = '/data/qifan/projects/EndtoEnd/results/slabBench/slabPoly_0.0_100'
    MCShape = (256, 63, 63)
    IPBResolution = 0.1
    convWindowSize = 15
    parallelBeam = getConvDose(convWindowSize, IPBFolder, MCShape)
    centerIdx = int(MCShape[1] / 2)
    parallelCenterline = parallelBeam[:, centerIdx, centerIdx]

    SAD = 200
    r0 = SAD - 12.8
    depth = r0 + (np.arange(MCShape[0]) + 0.5) * IPBResolution
    scale = SAD**2 / depth**2
    parallelCenterline = parallelCenterline * scale
    IPBdepth = (np.arange(MCShape[0]) + 0.5) * IPBResolution
    CCCSdepth = (np.arange(shape[1]) + 0.5) * CCCSResolution
    CCCSCenterLineInterpolated = np.interp(IPBdepth, CCCSdepth, CCCSCenterLine)
    # then we normalize to minimize the L2 distance
    sum_parallelCenterline_square = np.sum(np.square(parallelCenterline))
    sum_parallelCenterline_CCCSCenterLine = np.sum(parallelCenterline*CCCSCenterLineInterpolated)
    factor = sum_parallelCenterline_CCCSCenterLine / sum_parallelCenterline_square

    if False:
        parallelCenterline *= factor
        plt.plot(depth, CCCSCenterLine)
        plt.plot(depth, parallelCenterline)
        file = os.path.join(IPBFolder, 'lateralSAD200Window10.png')
        plt.savefig(file)
        plt.clf()
    
    parallelBeam *= factor

    # then we examine the lateral profile
    depth = 10  # cm
    depthCCCS = int(depth / CCCSResolution)
    depthParallel = int(depth / IPBResolution)
    CCCSLateral = CCCSDose[coord1, depthCCCS, :]
    parallelLateral = parallelBeam[depthParallel, centerIdx, :]
    r = r0 + depthParallel*IPBResolution
    scale_ = (SAD / r)**2
    parallelLateral = parallelLateral * scale_

    CCCSAxis = (np.arange(shape[2]) - coord2 - 1) * CCCSResolution
    parallelAxis = (np.arange(MCShape[2]) - (MCShape[2]-1) / 2) * IPBResolution * 20 / 15
    # scale the parallelAxis
    parallelAxis *= r / SAD

    intervalCCCS = 16
    intervalParallel = 40
    # CCCSPre = int((shape[2] - intervalCCCS) / 2)
    # CCCSPost = int((shape[2] + intervalCCCS) / 2)
    CCCSPre = int(coord2 - intervalCCCS / 2 + 1)
    CCCSPost = int(coord2 + intervalCCCS / 2 + 1)
    parallelPre = int((MCShape[2] - intervalParallel) / 2)
    parallelPost = int((MCShape[2] + intervalParallel) / 2)
    plt.plot(CCCSAxis[CCCSPre: CCCSPost], CCCSLateral[CCCSPre: CCCSPost])
    plt.plot(parallelAxis[parallelPre: parallelPost], parallelLateral[parallelPre: parallelPost])
    plt.xlabel('lateral distance (cm)')
    plt.ylabel('dose (a.u.)')
    plt.legend(['CCCS', 'Monte Carlo'])
    plt.title('lateral beamlet 1.0 cm')
    file = os.path.join(IPBFolder, 'lateralSAD200Window20.png')
    plt.savefig(file)
    plt.clf()


if __name__ == '__main__':
    # IPBConvExam()
    # lateralProfileExamineConv()
    # beamlet2_0Longitudinal()
    lateralProfileExamineConv_beamlet2()