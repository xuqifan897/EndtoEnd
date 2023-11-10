import os
import numpy as np

def viewDens():
    resultFolder = '/data/qifan/projects/EndtoEnd/results/CCCSclone/results'
    texFile = os.path.join(resultFolder, 'texDens_log.bin')
    resultFile = os.path.join(resultFolder, 'density.bin')
    dataType = np.float32
    texArray = np.fromfile(texFile, dtype=dataType)
    resultArray = np.fromfile(resultFile, dtype=dataType)
    diff = np.sum(np.abs(texArray - resultArray))
    print(diff)


if __name__ == '__main__':
    viewDens()