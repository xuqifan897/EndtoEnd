import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from pointIPBval import Cyd2Cart
from scipy.signal import convolve2d


def RadDepthConstruct():
    """
    This function returns a 1-d array representing the radiological depth
    """
    resolution = 0.05
    slabs = [("adipose", 0.8), ("muscle", 0.8), ("bone", 0.8), ("muscle", 0.8), 
             ("lung", 4.8), ("muscle", 0.8), ("bone", 0.8), ("adipose", 0.8),
             ("bone", 0.8), ("muscle", 0.8), ("adipose", 0.8)]
    # slabs = [("water", 12.8)]
    LUTmat = {"water": 1., "adipose": 0.92, "muscle": 1.04, 
              "bone": 1.85, "lung": 0.25}
    # material, thickness, density
    layers = [(a, int(np.round(b/resolution)), LUTmat[a]) for a, b in slabs]
    DimZ = 0
    for mat, thick, dens in layers:
        DimZ += thick
    denseMap = np.zeros(DimZ, dtype=np.double)
    offset = 0
    for mat, thick, density in layers:
        denseMap[offset: offset+thick] = density
        offset += thick
    result = np.zeros_like(denseMap)
    prev = 0.
    for i in range(DimZ):
        result[i] = prev
        prev += denseMap[i]
    return result, denseMap


def ConvertMonoCylinder2Cartesian():
    """
    This function converts the 6 MeV monoenergetic slab dose 
    profile in cylindrical coordinates to Cartesian coordinates
    """
    cylinderFile = '/data/qifan/projects/EndtoEnd/results/spec6MV/E6.00/SD.bin'
    shape = (256, 200)
    datatype = np.double
    cylinder = np.fromfile(cylinderFile, dtype=datatype)
    cylinder = np.reshape(cylinder, shape)
    cartesian = Cyd2Cart(cylinder)
    expFolder = '/data/qifan/projects/EndtoEnd/results/Sept1Point/slabPoint'
    cartesianFile = os.path.join(expFolder, 'IPBFromCylinder.npy')
    np.save(cartesianFile, cartesian)


def slabPointSynthesis():
    """
    Here we implement the correct version of the slab IPB synthesis from point kernel
    """
    expFolder = '/data/qifan/projects/EndtoEnd/results/Sept1Point/slabPoint'
    if not os.path.isdir(expFolder):
        os.mkdir(expFolder)
    # Firstly, read the point kernel
    pointKernelPath = '/data/qifan/projects/EndtoEnd/' \
        'results/Sept1Point/pointKernel/SD.bin'
    pointKernelShape = (55, 50, 50)
    interactionPoint = (5, 24.5, 24.5)
    phantomSZ = 5
    phantomBottom = 50
    datatype = np.double
    pointKernel = np.fromfile(pointKernelPath, dtype=datatype)
    pointKernel = np.reshape(pointKernel, pointKernelShape)

    DepthRad, DenseMap = RadDepthConstruct()
    DimZ = DepthRad.shape[0]
    
    # construct interpolation function
    AxisZ = np.linspace(-phantomSZ, phantomBottom-1, phantomSZ+phantomBottom)
    AxisY = np.linspace(-interactionPoint[1], interactionPoint[1], pointKernelShape[1])
    AxisX = np.linspace(-interactionPoint[2], interactionPoint[2], pointKernelShape[2])
    interpolator = RegularGridInterpolator((AxisZ, AxisY, AxisX), 
        pointKernel, bounds_error=False, fill_value=0)
    
    def DepthRad2Idx(dd):
        """
        This function converts the radiological depth into index
        """
        result = 0
        while(DepthRad[result] < dd and result < DimZ-1):
            result += 1
        return result

    gridY, gridX = np.meshgrid(AxisY, AxisX)
    globalRes = np.zeros((DimZ, pointKernelShape[1], pointKernelShape[2]), dtype=datatype)
    mu = 2.770E-02  # attenuation coefficient of 6 MeV Xray in water
    resolution = 0.1  # cm
    # start convolution
    for i in range(DimZ):
        CenterDepthRad = DepthRad[i]
        LeadingDepthRad = CenterDepthRad - phantomSZ
        LeadingDepthRad = max(0., LeadingDepthRad)
        LeadingIdx = DepthRad2Idx(LeadingDepthRad)
        EndingDepthRad = CenterDepthRad + phantomBottom
        EndingDepthRad = min(DepthRad[-1], EndingDepthRad)
        EndingIdx = DepthRad2Idx(EndingDepthRad)

        SizeZ = EndingIdx - LeadingIdx
        localRes = np.zeros((SizeZ, pointKernelShape[1], 
            pointKernelShape[2]), dtype=datatype)
        for j in range(SizeZ):
            idx_j = LeadingIdx + j
            if idx_j != i:
                DeltaZPhys = idx_j - i
                DeltaZRad = DepthRad[idx_j] - DepthRad[i]
                factor = DeltaZRad / DeltaZPhys
                gridY_local = gridY * factor
                gridY_local = np.expand_dims(gridY_local, 2)
                gridX_local = gridX * factor
                gridX_local = np.expand_dims(gridX_local, 2)
                gridZ_local = np.ones_like(gridY_local) * DeltaZRad
                coords = np.concatenate((gridZ_local, gridY_local, gridX_local), axis=2)
            else:
                localDensity = DenseMap[idx_j]
                gridY_local = gridY * localDensity
                gridY_local = np.expand_dims(gridY_local, 2)
                gridX_local = gridX * localDensity
                gridX_local = np.expand_dims(gridX_local, 2)
                gridZ_local = np.zeros_like(gridY_local)
                coords = np.concatenate((gridZ_local, gridY_local, gridX_local), axis=2)

            slice = interpolator(coords)
            # then apply inverse square law
            squareRad = np.sum(np.square(coords), axis=2)
            squarePhys = np.square(gridY) + np.square(gridX) + (i-idx_j)**2
            slice *= squareRad / squarePhys
            localRes[j, :, :] = slice
        if False:
            # for debug purposes
            debugFile = os.path.join(expFolder, 'sourceDepth{:03d}.npy'.format(i))
            np.save(debugFile, localRes)
        attenuation = np.exp(-mu * DepthRad[i] * resolution)
        globalRes[LeadingIdx: EndingIdx, :, :] += localRes * attenuation * DenseMap[i]
        print('source depth: {}'.format(i))
    arrayFile = os.path.join(expFolder, 'IPBSynthedWithPoint.npy')
    np.save(arrayFile, globalRes)


def viewCenterline():
    """
    This function views and compares the centerlines
    """
    expFolder = '/data/qifan/projects/EndtoEnd/results/Sept1Point/slabPoint'

    if False:
        # show the centerline dose
        IPBSynthFile = os.path.join(expFolder, 'IPBSynthedWithPoint.npy')
        IPBSynth = np.load(IPBSynthFile)
        IPBSynthShape = IPBSynth.shape
        IPBSynthCenterIdx = int(IPBSynthShape[1]/2)
        IPBSynthCenterline = IPBSynth[:, IPBSynthCenterIdx, IPBSynthCenterIdx]
        IPBSynthCenterline = IPBSynthCenterline / np.max(IPBSynthCenterline)

        IPBMCFile = os.path.join(expFolder, 'IPBFromCylinder.npy')
        IPBMC = np.load(IPBMCFile)
        DepthRad, DenseRad = RadDepthConstruct()
        DenseRad = np.expand_dims(DenseRad, axis=(1, 2))
        IPBMC /= DenseRad
        IPBMCShape = IPBMC.shape
        IPBMCCenterIdx = int(IPBMCShape[1]/2)
        IPBMCCenterline = IPBMC[:, IPBMCCenterIdx, IPBMCCenterIdx]
        IPBMCCenterline = IPBMCCenterline / np.max(IPBMCCenterline)

        DimZ = IPBSynthShape[0]
        depth = np.arange(DimZ) * 0.1  # cm
        plt.plot(depth, IPBMCCenterline)
        plt.plot(depth, IPBSynthCenterline)
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.legend(['Monte Carlo', 'Synthesized with point dose kernel'])
        figureFile = os.path.join(expFolder, 'IPBSynthCenterline.png')
        plt.savefig(figureFile)
        plt.clf()
    if True:
        # show beam with finite width, by convolution
        IPBSynthFile = os.path.join(expFolder, 'IPBSynthedWithPoint.npy')
        IPBSynth = np.load(IPBSynthFile)
        IPBSynthShape = IPBSynth.shape

        IPBMCFile = os.path.join(expFolder, 'IPBFromCylinder.npy')
        IPBMC = np.load(IPBMCFile)
        DepthRad, DenseRad = RadDepthConstruct()
        DenseRad = np.expand_dims(DenseRad, axis=(1, 2))
        IPBMC /= DenseRad
        IPBMCShape = IPBMC.shape

        convWindowShape = (10, 10)
        convWindow = np.ones(convWindowShape, dtype=np.double)
        IPBSynthConv = np.zeros_like(IPBSynth)
        for i in range(IPBSynthShape[0]):
            IPBSynthConv[i, :, :] = convolve2d(IPBSynth[i, :, :], convWindow, 'same')
        IPBMCConv = np.zeros_like(IPBMC)
        for i in range(IPBMCShape[0]):
            IPBMCConv[i, :, :] = convolve2d(IPBMC[i, :, :], convWindow, 'same')
        IPBSynthCenterIdx = int(IPBSynthShape[1]/2)
        IPBSynthCenterline = IPBSynthConv[:, IPBSynthCenterIdx, IPBSynthCenterIdx]
        IPBSynthNorm = np.max(IPBSynthCenterline)
        IPBSynthConv /= IPBSynthNorm
        IPBMCCenterIdx = int(IPBMCShape[1]/2)
        IPBMCCenterline = IPBMCConv[:, IPBMCCenterIdx, IPBMCCenterIdx]
        IPBMCNorm = np.max(IPBMCCenterline)
        IPBMCConv /= IPBMCNorm

        # depth dose curve
        depth = np.arange(IPBSynthShape[0]) * 0.1  # cm
        plt.plot(depth, IPBMCCenterline)
        plt.plot(depth, IPBSynthCenterline)
        plt.xlabel('depth (cm)')
        plt.ylabel('dose (a.u.)')
        plt.legend(['Monte Carlo', 'Synthesized with point dose kernel'])
        figureFile = os.path.join(expFolder, 'IPBSynthConv.png')
        plt.savefig(figureFile)
        plt.clf()

        # lateral dose profile at 10 cm depth
        DepthIdx = 100
        IPBSynthProfile = IPBSynthConv[DepthIdx, IPBSynthCenterIdx, :]
        IPBSynthAxis = (np.arange(IPBSynthShape[2]) - (IPBSynthShape[2]-1)/2) * 0.1  # cm
        IPBMCProfile = IPBMCConv[DepthIdx, IPBMCCenterIdx, :]
        IPBMCAxis = (np.arange(IPBMCShape[2]) - (IPBMCShape[2]-1)/2) * 0.1  # cm
        plt.plot(IPBSynthAxis, IPBSynthProfile)
        plt.plot(IPBMCAxis, IPBMCProfile)
        plt.xlabel('lateral distance (cm)')
        plt.ylabel('dose (a.u.)')
        plt.legend(['Monte Carlo', 'Synthesized with point dose kernel'])
        figureFile = os.path.join(expFolder, 'IPBSynthLateral.png')
        plt.savefig(figureFile)
        plt.clf()


if __name__ == '__main__':
    # slabPointSynthesis()
    # ConvertMonoCylinder2Cartesian()
    viewCenterline()