import os
from os.path import join as pjoin
import numpy as np
from math import ceil

from volume import Volume
from fmaps import Fmaps, Beam
from materials import mat_map

if __name__ == "__main__":
    # CUSTOMIZABLE PARAMETERS
    # center of phantom surface at origin assumed
    # beam fired in +Y direction assumed
    OUTPUT_DIR = './phantom_output'
    blt_size  = (10, 10)        # beamlet size (X, Z)          [unit: cm]
    fmap_size = (1, 1)          # fluence map dimensions (X, Z)
    penumbra  = (2, 2)          # size outside of field (X, Z) [unit: cm]
    voxelsize = (0.1, 0.1, 0.1) # (X,Y,Z)                      [unit: cm]
    length    = 25              # length of phantom (Y-axis)   [unit: cm]
    iso       = (0, 0, 0)       # (X,Y,Z)                      [unit: cm]
    sad       = 100             # source-iso distance          [unit: cm]
    density   = 1.0             # uniform density              [unit: g/cm^3]
    fluence   = 1.0             # uniform fluence              [unit: arb]

    # construct density array
    fieldsize = np.multiply(blt_size, fmap_size)
    size2d = np.add(fieldsize, np.multiply(2, penumbra))
    size = (ceil(size2d[0]/voxelsize[0]),  ceil(length/voxelsize[1]),  ceil(size2d[1]/voxelsize[2]))
    dens = density*np.ones((size[::-1]))
    vol = Volume.CenterAt(dens, center=(0,0,0), voxelsize=voxelsize)
    vol.start[1] = 0 # place surface at origin
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    vol.generate(pjoin(OUTPUT_DIR, 'density.h5'))
    print(vol)

    # create beam and fluence map
    fmaps = Fmaps()
    fmap_wts = fluence*np.ones(fmap_size)
    beam = Beam(fmap_wts, iso=iso, sad=sad, beamlet_size=blt_size)
    fmaps.addBeam(beam)
    fmaps.generate(pjoin(OUTPUT_DIR, 'fmaps.h5'))
    print(beam)

