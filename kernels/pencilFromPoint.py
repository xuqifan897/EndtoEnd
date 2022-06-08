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


def calcPencilKernel(mu, kernel):
    # mu = args.mu[1]
    # kernel = kernel6MeV
    depth = np.arange(1, 21, dtype=np.float32)
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

    # slice = dose[10, :]
    # lateral_distance = lateral_distance.squeeze()
    # plt.plot(lateral_distance, slice)
    # plt.show()

    return dose * lateral_distance.squeeze(axis=2)


def viewKernel(pencilKernel):
    slice2 = pencilKernel[1, :]
    slice5 = pencilKernel[4, :]
    slice10 = pencilKernel[9, :]
    slice15 = pencilKernel[14, :]
    slice20 = pencilKernel[19, :]

    lateral_distance = np.arange(0.05, 1.65, 0.05, dtype=np.float32)
    plt.plot(lateral_distance, slice2)
    plt.plot(lateral_distance, slice5)
    plt.plot(lateral_distance, slice10)
    plt.plot(lateral_distance, slice15)
    plt.plot(lateral_distance, slice20)
    plt.xlabel('lateral distance / cm')
    plt.ylabel('dose (a.u.)')
    plt.legend(['depth = 2cm', 'depth = 5cm', 'depth = 10cm', 'depth = 15cm', 'depth = 20cm'])
    plt.title('lateral k profile')
    # plt.show()
    plt.show()

    slice2 /= np.max(slice2)
    slice5 /= np.max(slice5)
    slice10 /= np.max(slice10)
    slice15 /= np.max(slice15)
    slice20 /= np.max(slice20)

    plt.plot(lateral_distance, slice2)
    plt.plot(lateral_distance, slice5)
    plt.plot(lateral_distance, slice10)
    plt.plot(lateral_distance, slice15)
    plt.plot(lateral_distance, slice20)
    plt.xlabel('lateral distance / cm')
    plt.ylabel('dose (a.u.)')
    plt.legend(['depth = 2cm', 'depth = 5cm', 'depth = 10cm', 'depth = 15cm', 'depth = 20cm'])
    plt.title('lateral k profile (normalized)')
    plt.show()


def averageKernel(pencilKernel):
    slice2 = pencilKernel[1, :]
    slice5 = pencilKernel[4, :]
    slice10 = pencilKernel[9, :]
    slice15 = pencilKernel[14, :]
    slice20 = pencilKernel[19, :]
    slice2 /= np.max(slice2)
    slice5 /= np.max(slice5)
    slice10 /= np.max(slice10)
    slice15 /= np.max(slice15)
    slice20 /= np.max(slice20)
    slice = (slice2 + slice5 + slice10 + slice15 + slice20) / 5
    return slice


def fitKernel(pencilKernel):
    averagedKernel = averageKernel(pencilKernel)
    lateral_distance = np.arange(0.05, 1.65, 0.05, dtype=np.float32)
    Ainit, Binit, ainit, binit = \
        solve(lateral_distance, averagedKernel, args.Ainit, args.Binit, args.ainit, args.binit)
    fitted = Ainit * np.exp(-ainit * lateral_distance) + Binit * np.exp(-binit * lateral_distance)
    AInit, BInit = Ainit / (Ainit + Binit), Binit / (Ainit + Binit)
    print('{},{},{},{}'.format(AInit, BInit, ainit, binit))
    plt.plot(lateral_distance, averagedKernel)
    plt.plot(lateral_distance, fitted)
    plt.xlabel('lateral distaince / cm')
    plt.ylabel('dose (a.u.)')
    plt.legend(['ground truth', 'fit'])
    plt.title('the fitting of the lateral kernel')
    plt.show()


def solve(xvalues, yvalues, Ainit, Binit, ainit, binit):
    num_iters = 5000
    errors = []
    # AB_step = 0.001
    # ab_step = 0.001
    for i in range(num_iters):
        first_part = np.exp(-ainit * xvalues)
        second_part = np.exp(-binit * xvalues)
        yestimates = Ainit * first_part + Binit * second_part
        error = np.sum((yvalues - yestimates)**2)
        # print(i, error)
        errors.append(error)
        A_grad = np.sum((yestimates - yvalues) * first_part)
        B_grad = np.sum((yestimates - yvalues) * second_part)
        a_grad = np.sum((yestimates - yvalues) * (-Ainit * xvalues * first_part))
        b_grad = np.sum((yestimates - yvalues) * (-Binit * xvalues * second_part))
        ABnorm = np.sqrt(A_grad ** 2 + B_grad ** 2)
        abnorm = np.sqrt(a_grad ** 2 + b_grad ** 2)

        if i < 200:
            AB_step = 0.005
            ab_step = 0.005
        elif i < 2000:
            AB_step = 0.002
            ab_step = 0.002
        elif i < 3000:
            AB_step = 0.001
            ab_step = 0.001
        elif i < 5000:
            AB_step = 0.0005
            ab_step = 0.0005

        Ainit -= A_grad / ABnorm * AB_step
        Binit -= B_grad / ABnorm * AB_step
        ainit -= a_grad / abnorm * ab_step
        binit -= b_grad / abnorm * ab_step

    # print(Ainit, Binit, ainit, binit)
    # plt.plot(list(range(num_iters)), errors)
    # plt.xlabel('iterations')
    # plt.ylabel('loss')
    # plt.yscale('log')
    # plt.show()

    return Ainit, Binit, ainit, binit


def main():
    getargs()
    init_kernels()

    kernels = [kernel4MeV, kernel6MeV, kernel10MeV, kernel15MeV]
    for i in range(4):
        pencilKernel = calcPencilKernel(args.mu[i], kernels[i])
        fitKernel(pencilKernel)


if __name__ == '__main__':
    main()