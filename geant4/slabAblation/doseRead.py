import os
import numpy as np
from scipy import signal
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

muTableInit()


def doseRead():
    """
    This function reads the dose
    """
    doseFolder = '/data/qifan/projects/EndtoEnd4/results' \
        '/InhomoJuly6/phantomDose'
    doseFile = os.path.join(doseFolder, 'array.bin')
    shape = (199, 199, 256)
    res = (0.1, 0.1, 0.1 )  # cm
    edep = np.fromfile(doseFile, dtype=np.double)
    edep = np.reshape(edep, shape)

    # construct density matrix
    density = np.zeros(shape)
    thicknesses = [16, 16, 16, 16, 96, 16, 16, 16, 16, 16, 16]
    MM = ['adipose', 'muscle', 'bone', 'muscle', 'lung', 'muscle', 
        'bone', 'adipose', 'bone', 'muscle', 'adipose']
    prev = 0
    for thick, mat in zip(thicknesses, MM):
        density[:, :, prev:prev+thick] = materials[mat][0]
        prev += thick
    
    depth = np.arange(shape[2]) * res[2]
    dose = edep / density
    file = os.path.join(doseFolder, 'IPB.png')
    plt.plot(depth, dose[99, 99, :])
    plt.savefig(file)
    plt.clf()

    # then, we simulate the 0.5cm, 1cm, and 2cm data
    kernelSizes = [5, 10, 20]
    for ks in kernelSizes:
        kernelSizePrecessing(dose, shape, ks, doseFolder)
    


def kernelSizePrecessing(dose, shape, kernelSize, doseFolder):
    newShape = (shape[0] + kernelSize - 1, shape[1] + kernelSize - 1, shape[2])
    newDose = np.zeros(newShape, dtype=np.double)
    kernel = np.ones((kernelSize, kernelSize))
    for i in range(shape[2]):
        newDose[:, :, i] = signal.convolve2d(dose[:, :, i], kernel)

    centerIdx = int((shape[0] + kernelSize) / 2 - 1)
    centralLine = newDose[centerIdx, centerIdx, :]
    res = 0.1
    file = os.path.join(doseFolder, 'window{}.png'.format(kernelSize*res))
    depth = np.arange(shape[2]) * res
    plt.plot(depth, centralLine)
    plt.savefig(file)
    plt.clf()


if __name__ == '__main__':
    doseRead()