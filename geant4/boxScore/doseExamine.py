import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    shape = (20, 256, 256)
    folder = '/data/qifan/projects/EndtoEnd/results/slabBench/patient1_g4'
    data = os.path.join(folder, 'SD.bin')
    array = np.fromfile(data, dtype=np.double)
    array = np.reshape(array, shape)

    print(np.sum(array))

    # centerLine = array[:, 128, 128]
    # plt.plot(np.arange(shape[0]), centerLine)
    # figureFile = os.path.join(folder, 'ddp.png')
    # plt.savefig(figureFile)
    # plt.clf()

if __name__ == '__main__':
    main()