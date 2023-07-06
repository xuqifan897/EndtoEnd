import os
import numpy as np
import matplotlib.pyplot as plt

# below are the material properties.
# density are in g/cm^3
materials = {'adipose': [0.95, {'H': 0.114, 'C': 0.598, 
        'N': 0.007, 'O': 0.278, 'Na': 0.001, 'S': 0.001, 'Cl': 0.001, }],
    'bone': [1.85, {'H': 0.064, 'C': 0.278, 'N': 0.027, 'O': 0.41, 
        'Mg': 0.002, 'P': 0.07, 'S': 0.002, 'Ca': 0.147}],
    'lung': [1.040, {'H': 0.105, 'C': 0.083, 'N': 0.023, 
        'O': 0.779, 'Na': 0.002, 'P': 0.001, 'S': 0.002, 'Cl': 0.003, 'K': 0.002}],
    'muscle': [1.05, {'H': 0.102, 'C': 0.143, 'N': 0.034, 
        'O': 0.71, 'Na': 0.001, 'P': 0.002, 'S': 0.003, 'Cl': 0.001, 'K': 0.004}]}

# data are under 6MeV, units are in cm^2/g
TermaTable = {'H': 4.498e-2, 'C': 2.469e-2, 'N': 2.511e-02, 'O': 2.552e-02,
    'Na': 2.559e-02, 'Mg': 2.681e-02, 'P': 2.747e-02, 'S': 2.872e-02,
    'Cl': 2.798e-02, 'K': 2.915e-02, 'Ca':3.035e-02}

muTable = None

def muTableInit():
    """
    This function initializes the muTable
    """
    # firstly, check the unity of the material decomposition
    for key, item in materials.items():
        elementTable = item[1]
        fractions = list(elementTable.values())
        unity = np.sum(fractions)
        assert np.abs(unity - 1) < 1e-3, "the fractions do not add to unity"

    global muTable
    muTable = {}
    for key, value in materials.items():
        density = value[0]
        elements = value[1]
        mu = 0
        for element, fraction in elements.items():
            elementDensity = density * fraction
            elementMuRho = TermaTable[element]
            elementMu = elementDensity * elementMuRho
            mu += elementMu
        muTable[key] = mu
    # print(muTable)

muTableInit()

def showAnalytical():
    """
    This function shows the dose from the proposed analytical method
    """
    # first step, phantom construction. The thickness unit is in mm
    layers = [('adipose', 16), ('muscle', 16), ('bone', 16), ('muscle', 16),
        ('lung', 96), ('muscle', 16), ('bone', 16), ('adipose', 16),
        ('bone', 16), ('muscle', 16), ('adipose', 16)]
    millimeters = 0
    for a in layers:
        millimeters += a[1]
    depth = np.arange(millimeters) * 0.1

    densityMap = np.zeros(millimeters)
    muMap = np.zeros(millimeters)
    cumuThickness = 0
    for material, thickness in layers:
        mu = muTable[material]
        muMap[cumuThickness: cumuThickness + thickness] = mu

        density = materials[material][0]
        densityMap[cumuThickness: cumuThickness + thickness] = density
        cumuThickness += thickness
    
    # # for debug purposes, to show the ratio between mu and density
    # ratio = muMap / densityMap
    # plt.plot(depth, ratio)
    # figureFile = '/data/qifan/projects/EndtoEnd4/results/' \
    #     'slab6MeVAccu/muOverDens.png'
    # plt.savefig(figureFile)
    # plt.clf()
    # exit()
    
    muMapIntegral = np.cumsum(muMap) * 0.1  # thickness: 0.1cm
    Terma = np.exp(-muMapIntegral)
    TermaMu = Terma * muMap
    # which reflects the equivalent radiology path lengh in water, in mm
    densityMapIntegral = np.cumsum(densityMap)

    # # for debug purposes, show termaMu
    # plt.plot(depth, Terma)
    # plt.savefig('/data/qifan/projects/EndtoEnd4/results/slab6MeVAccu/Terma.png')
    # plt.clf()
    # exit()

    # load Terma to dose kernel
    kernelPath = '/data/qifan/projects/EndtoEnd4/results/point6MeV1e7Step/array.bin'
    kernelShape = (199, 199, 400)
    kernelData = np.fromfile(kernelPath, dtype=np.float64)
    kernelData = np.reshape(kernelData, kernelShape)
    reducedKernelData = np.sum(kernelData, axis=(0, 1))
    # kernelDepth = np.arange(kernelShape[2]) * 0.1
    # plt.plot(kernelDepth, reducedKernelData)
    # plt.show()

    kernelLeading = 50  # the leading 100 pixels are backscattering
    kernelSize = kernelShape[2]
    radDisMat = np.zeros((millimeters, millimeters))
    for i in range(millimeters):
        radDisMat[i, :] = densityMapIntegral - densityMapIntegral[i]
    radDisMat += kernelLeading
    weights = np.zeros_like(radDisMat)
    for i in range(millimeters):
        for j in range(millimeters):
            coord = radDisMat[i, j]
            coord_lower = int(np.floor(coord))
            coord_higher = coord_lower + 1
            coeff_lower = coord_higher - coord
            coeff_higher = coord - coord_lower
            if coord_lower >= 0 and coord_lower < kernelSize:
                lower_value = reducedKernelData[coord_lower]
            else:
                lower_value = 0
            if coord_higher >= 0 and coord_higher < kernelSize:
                higher_value = reducedKernelData[coord_higher]
            else:
                higher_value = 0
            value = lower_value * coeff_lower + higher_value * coeff_higher
            weights[i, j] = value
    
    DoseContri = np.zeros_like(weights)
    for i in range(millimeters):
        DoseContri[i, :] = weights[i, :] * TermaMu[i]
    AnalyticalDose = np.sum(DoseContri, axis=0)
    # AnalyticalDose *= muMap / densityMap
    
    doseFile = '/data/qifan/projects/EndtoEnd4/results/slab6MeVAccu/array.bin'
    shape = (199, 199, 256)
    array = np.fromfile(doseFile, dtype=np.float64)
    array = np.reshape(array, shape)
    MonteCarloEdep = np.sum(array, axis=(0, 1))
    MonteCarloDose = MonteCarloEdep / densityMap

    # normalize
    AnalyticalDose /= AnalyticalDose[100]
    MonteCarloDose /= MonteCarloDose[100]

    plt.plot(depth, MonteCarloDose)
    plt.plot(depth, AnalyticalDose)
    plt.legend(['Monte Carlo', 'Analytical'])
    plt.xlabel('depth (cm)')
    plt.ylabel('dose (a.u.)')
    # plt.show()
    graphPath = '/data/qifan/projects/EndtoEnd4/results/slab6MeVAccu/doseComp.png'
    plt.savefig(graphPath)
    plt.clf()


def examineDose():
    # first step, phantom construction. The thickness unit is in mm
    layers = [('adipose', 16), ('muscle', 16), ('bone', 16), ('muscle', 16),
        ('lung', 96), ('muscle', 16), ('bone', 16), ('adipose', 16),
        ('bone', 16), ('muscle', 16), ('adipose', 16)]
    millimeters = 0
    for a in layers:
        millimeters += a[1]
    depth = np.arange(millimeters) * 0.1

    densityMap = np.zeros(millimeters)
    muMap = np.zeros(millimeters)
    cumuThickness = 0
    for material, thickness in layers:
        mu = muTable[material]
        muMap[cumuThickness: cumuThickness + thickness] = mu

        density = materials[material][0]
        densityMap[cumuThickness: cumuThickness + thickness] = density
        cumuThickness += thickness
    
    Folder = '/data/qifan/projects/EndtoEnd4/results/slab6MeVAccu'
    doseFile = os.path.join(Folder, 'array.bin')
    shape = (199, 199, 256)
    array = np.fromfile(doseFile, dtype=np.float64)
    array = np.reshape(array, shape)

    # process = 'marginal'
    process = 'central'

    if process == 'marginal':
        array = np.sum(array, axis=(0, 1))
        figurePath = os.path.join(Folder, 'doseMarginal.png')
    elif process == 'central':
        array = array[99, 99, :]
        figurePath = os.path.join(Folder, 'doseCentral.png')
    array /= densityMap
    depth = np.arange(256) * 0.1
    plt.plot(depth, array)
    plt.savefig(figurePath)
    plt.clf()


def doseExamine():
    """
    After the first trail, we found that the analytical dose does not 
    match the Monte Carlo dose in the bone regions, where there are bumps 
    in the analytical dose, but not in the Monte Carlo dose. Here we are 
    going to examine the reason or problem behind this.
    """
    # first step, phantom construction. The thickness unit is in mm
    layers = [('adipose', 16), ('muscle', 16), ('bone', 16), ('muscle', 16),
        ('lung', 96), ('muscle', 16), ('bone', 16), ('adipose', 16),
        ('bone', 16), ('muscle', 16), ('adipose', 16)]
    millimeters = 0
    for a in layers:
        millimeters += a[1]
    depth = np.arange(millimeters) * 0.1

    densityMap = np.zeros(millimeters)
    muMap = np.zeros(millimeters)
    cumuThickness = 0
    for material, thickness in layers:
        mu = muTable[material]
        muMap[cumuThickness: cumuThickness + thickness] = mu

        density = materials[material][0]
        densityMap[cumuThickness: cumuThickness + thickness] = density
        cumuThickness += thickness
    muMapIntegral = np.cumsum(muMap) * 0.1  # thickness: 0.1cm
    Terma = np.exp(-muMapIntegral)
    
    resultFolder = '/data/qifan/projects/EndtoEnd4/results/slab6MeVAccu'

    phase = 3
    if phase == 1:
        TermaFigPath = os.path.join(resultFolder, 'Terma.png')
        plt.plot(depth, Terma)
        plt.savefig(TermaFigPath)
        plt.clf()
    elif phase == 2:
        TermaDiff = np.zeros_like(Terma)
        TermaDiff[:-1] = - np.diff(Terma)
        TermaDiff[-1] = TermaDiff[-2]
        TermaDiff /= np.max(TermaDiff)

        TermaMu = Terma * muMap
        TermaMu /= np.max(TermaMu)
        plt.plot(depth, TermaDiff)
        plt.plot(depth, TermaMu)
        plt.legend(['TermaDiff', 'TermaMu'])
        figurePath = os.path.join(resultFolder, 'TermaDiff.png')
        plt.savefig(figurePath)
        plt.clf()
    elif phase == 3:
        TermaMuRho = Terma * muMap / densityMap
        TermaMuRhoPath = os.path.join(resultFolder, 'TermaMuRho.png')
        plt.plot(depth, TermaMuRho)
        plt.savefig(TermaMuRhoPath)
        plt.clf()



if __name__ == '__main__':
    # showAnalytical()
    examineDose()
    # doseExamine()