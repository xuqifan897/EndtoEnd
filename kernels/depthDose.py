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


class pencilKernel():
    def __init__(self, A=0, B=0, a=0, b=0):
        self.A = A
        self.B = B
        self.a = a
        self.b = b

    def get(self, lateral_distance):
        result = self.A * np.exp(-self.a * lateral_distance) + \
                 self.B * np.exp(-self.b * lateral_distance)
        return result


kernel4MeV = Kernel()
kernel6MeV = Kernel()
kernel10MeV = Kernel()
kernel15MeV = Kernel()

pencilKernel4MeV = pencilKernel()
pencilKernel6MeV = pencilKernel()
pencilKernel10MeV = pencilKernel()
pencilKernel15MeV = pencilKernel()

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


def init_pencil_kernels():
    FCBB_file = os.path.join(args.folder, args.FCBB_file)
    with open(FCBB_file, 'r') as f:
        lines = f.readlines()
    # global pencilKernel4MeV, pencilKernel6MeV, pencilKernel10MeV, pencilKernel15MeV
    pencilKernels = [pencilKernel4MeV, pencilKernel6MeV, pencilKernel10MeV, pencilKernel15MeV]
    for i in range(4):
        line = lines[i]
        line = line.split(',')
        line = [float(a) for a in line]
        kernel = pencilKernels[i]
        kernel.A, kernel.B, kernel.a, kernel.b = line


def calcPencilKernel(mu, kernel, depth):
    # mu = args.mu[1]
    # kernel = kernel6MeV
    # depth = np.arange(0.1, 20.1, 0.1, dtype=np.float32)
    depth = np.expand_dims(depth, axis=(1, 2))
    lateral_distance = np.arange(0.05, 1.65, 0.05, dtype=np.float32)
    lateral_distance = np.expand_dims(lateral_distance, axis=(0, 2))
    angle = kernel.angles
    angle = np.expand_dims(angle, axis=(0, 1))
    shape = (depth.size, lateral_distance.size, angle.size)

    source_depth = depth - lateral_distance * np.tan(np.pi/2 - angle)
    source_depth_prev = np.zeros_like(source_depth)
    source_depth_prev[:, :, 1:] = source_depth[:, :, :-1]
    source_depth_prev[source_depth_prev < 0] = 0
    terma = np.exp(-mu * source_depth_prev) - np.exp(-mu * source_depth)
    terma[source_depth < 0] = 0

    distance = lateral_distance / np.cos(np.pi/2 - angle)

    Atheta, Btheta, atheta, btheta = kernel.Atheta, kernel.Btheta, kernel.atheta, kernel.btheta
    Atheta, Btheta, atheta, btheta = np.expand_dims(Atheta, axis=(0, 1)), np.expand_dims(Btheta, axis=(0, 1)), \
        np.expand_dims(atheta, axis=(0, 1)), np.expand_dims(btheta, axis=(0, 1))

    dose = terma * (Atheta * np.exp(- atheta * distance) + Btheta * np.exp(- btheta * distance)) / (distance * distance)
    dose = np.sum(dose, axis=2)
    lateral_distance = lateral_distance.squeeze(axis=2)
    pencilKernel = dose * lateral_distance
    return pencilKernel


def main():
    getargs()
    init_kernels()
    init_pencil_kernels()
    kernels = [kernel4MeV, kernel6MeV, kernel10MeV, kernel15MeV]
    pencilKernels = [pencilKernel4MeV, pencilKernel6MeV, pencilKernel10MeV, pencilKernel15MeV]
    lateral_distance = np.arange(0.05, 1.65, 0.05, dtype=np.float32)
    depth = np.arange(0.1, 40.1, 0.1, dtype=np.float32)
    depthDoses = []
    pks = []
    for i in range(4):
        pk = calcPencilKernel(args.mu[i], kernels[i], depth)
        pks.append(pk)
        pencilKernel = pencilKernels[i]
        pencilKernel_ = pencilKernel.get(lateral_distance)
        pk_sum = np.sum(pk, axis=1)
        pencilKernel_sum = np.sum(pencilKernel_)
        depthDose = pk_sum / pencilKernel_sum
        depthDoses.append(depthDose)
        plt.plot(depth, depthDose)
    plt.show()

    # pk = pks[1]
    # depthDose6MeV = depthDoses[1]
    # pencilKernel_6MeV = pencilKernel6MeV.get(lateral_distance)
    # indices = [20, 50, 100, 150, 200]
    # for i in range(5):
    #     idx = indices[i]
    #     CAX = depthDose6MeV[idx]
    #     pk_slice = pk[idx, :]
    #     plt.plot(lateral_distance, pk_slice)
    #     plt.plot(lateral_distance, pencilKernel_6MeV * CAX)
    #     plt.show()

    depthDosefile = os.path.join(args.folder, 'depthDose.csv')
    lines = ''
    for i in range(depth.size):
        line = '{},{},{},{},{}\n'.format(depth[i], depthDoses[0][i], depthDoses[1][i],
                                       depthDoses[2][i], depthDoses[3][i])
        lines = lines + line
    lines = lines[:-1]
    with open(depthDosefile, 'w') as f:
        f.writelines(lines)


if __name__ == '__main__':
    main()