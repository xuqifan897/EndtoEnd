import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

args = None
numAngles = 48

class Kernel():
    def __init__(self, length=numAngles, step=1.875*np.pi/180):
        self.angles = (np.arange(length, dtype=np.float32) * 2 + 1) * step
        self.Atheta = np.zeros(length, dtype=np.float32)
        self.Btheta = np.zeros(length, dtype=np.float32)
        self.atheta = np.zeros(length, dtype=np.float32)
        self.btheta = np.zeros(length, dtype=np.float32)

kernel4MeV = Kernel()
kernel6MeV = Kernel()
kernel10MeV = Kernel()
kernel15MeV = Kernel()

def getargs():
    parser = argparse.ArgumentParser(description='the options for calculating pencil kernel from point kernel')
    # parser.add_argument('--folder', type=str, default='./kernels', help='the folder that contains the kernelfiles')
    parser.add_argument('--folder', type=str, default='.', help='the folder that contains the kernelfiles')
    parser.add_argument('--Atheta-file', type=str, default='upperATheta.csv', help='Atheta')
    parser.add_argument('--Btheta-file', type=str, default='upperBTheta.csv', help='Btheta')
    parser.add_argument('--atheta-file', type=str, default='lowerATheta.csv', help='atheta')
    parser.add_argument('--btheta-file', type=str, default='lowerBTheta.csv', help='btheta')
    parser.add_argument('--FCBB-file', type=str, default='FCBBkernel.csv')
    parser.add_argument('--muRho', type=float, nargs='+', default=[3.403E-02, 2.770E-02, 2.219E-02, 1.941E-02, 0],
                        help='mass linear attenuation coefficient, in cm^2/g, for energies 4, 6, 10, 15, 24')
    parser.add_argument('--rho', type=float, default=1, help='the density of water, in g/cm^3')
    parser.add_argument('--Ainit', type=float, default=0.99, help='the initialization value of A')
    parser.add_argument('--Binit', type=float, default=0.01, help='the initialization value of B')
    parser.add_argument('--ainit', type=float, default=6.64, help='the initialization value of a')
    parser.add_argument('--binit', type=float, default=0.158, help='the initialization value of b')
    global args
    _args = parser.parse_args()
    _args.mu = [a * _args.rho for a in _args.muRho]
    args = _args


def init_kernels():
    Atheta_file = os.path.join(args.folder, args.Atheta_file)
    Btheta_file = os.path.join(args.folder, args.Btheta_file)
    atheta_file = os.path.join(args.folder, args.atheta_file)
    btheta_file = os.path.join(args.folder, args.btheta_file)

    global kernel4MeV, kernel6MeV, kernel10MeV, kernel15MeV
    kernels = [kernel4MeV, kernel6MeV, kernel10MeV, kernel15MeV]

    with open(Atheta_file, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        angle, *values, value24 = line.split(',')
        for kernel, value in zip(kernels, values):
            kernel.Atheta[i] = value

    with open(Btheta_file, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        angle, *values, value24 = line.split(',')
        for kernel, value in zip(kernels, values):
            kernel.Btheta[i] = float(value)

    with open(atheta_file, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        angle, *values, value24 = line.split(',')
        for kernel, value in zip(kernels, values):
            kernel.atheta[i] = value

    with open(btheta_file, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        angle, *values, value24 = line.split(',')
        for kernel, value in zip(kernels, values):
            kernel.btheta[i] = value


def normCheck(kernel):
    angleStep = 3.75 * np.pi / 180
    num_angles = 48
    angles0 = np.arange(num_angles, dtype=np.float32) * angleStep
    angles1 = angles0 + angleStep
    solidAngles = 2 * np.pi * (np.cos(angles0) - np.cos(angles1))

    eps = 1e-5
    first_term = kernel.Atheta / (kernel.atheta + eps)
    first_term[kernel.Atheta == 0] = 0
    second_term = kernel.Btheta / kernel.btheta
    value = solidAngles * (first_term + second_term)
    value = np.sum(value)
    print(np.sum(solidAngles) / (4 * np.pi), value)


def main():
    getargs()
    init_kernels()
    kernels = [kernel4MeV, kernel6MeV, kernel10MeV, kernel15MeV]
    for kernel in kernels:
        normCheck(kernel)


if __name__ == '__main__':
    main()