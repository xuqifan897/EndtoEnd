import os
import numpy as np
import matplotlib.pyplot as plt

phantom1 = [("adipose", 0.8),
            ("muscle", 0.8),
            ("bone", 0.8),
            ("muscle", 0.8),
            ("lung", 4.8),
            ("muscle", 0.8),
            ("bone", 0.8),
            ("adipose", 0.8),
            ("bone", 0.8),
            ("muscle", 0.8),
            ("adipose", 0.8)]

density = {"adipose": 0.92, "muscle": 1.04, "bone": 1.85, "lung": 0.25}

def readDose():
    """
    This function reads the dose from the Monte Carlo simulation
    """
    resultFolder = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14'
    file = os.path.join(resultFolder, 'SD.bin')
    array = np.fromfile(file, dtype=np.float64)
    shape = (256, 200)
    array = np.reshape(array, shape)

    # normalize the dose array against the density
    den = np.zeros_like(array)
    count = 0
    resZ = 0.05
    for material, thickness in phantom1:
        d = density[material]
        dim = round(thickness / resZ)
        for i in range(count, count+dim):
            den[i, :] = d
        count += dim
    
    dose = array / den
    resultFolder = '/data/qifan/projects/EndtoEnd/results/InhomoJuly20/slabAug14'
        
    # view central dose
    slice = dose[:, 0]
    depth = np.arange(shape[0]) * 0.1  # cm
    plt.plot(depth, slice)
    plt.xlabel('depth / cm')
    plt.ylabel('dose (a.u.)')
    plt.title('centerline dose')
    file = os.path.join(resultFolder, 'centerlineDose.png')
    plt.savefig(file)
    plt.clf()

    # view partial dose
    partialDose = np.sum(dose, axis=1)
    plt.plot(depth, partialDose)
    plt.xlabel('depth / cm')
    plt.ylabel('(a.u.)')
    plt.title('partial dose')
    file = os.path.join(resultFolder, 'partialDose.png')
    plt.savefig(file)
    plt.clf()
    

if __name__ == '__main__':
    readDose()