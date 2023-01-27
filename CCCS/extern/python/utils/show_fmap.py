#!/usr/bin/env python

import os
import sys
import math
import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

path = './fmo_results.h5'
if len(sys.argv) > 1:
    path = sys.argv[1]

fmaps = []
beam_indices = []
couch_gantry = []
with h5py.File(path, 'r') as h5file:
    for beam_idx, beam_meta in sorted(h5file['beams']['metadata'].items(), key=lambda t: int(t[0])):
        fmap = beam_meta.attrs['fmap_weights']
        fmaps.append(fmap)
        beam_indices.append(int(beam_idx))
        couch_gantry.append( (beam_meta.attrs['beam_specs']['couch_rot_rad'], beam_meta.attrs['beam_specs']['gantry_rot_rad']) )

nbeams = len(fmaps)
dimlong= math.ceil(math.sqrt(nbeams))
dim1d = math.ceil(nbeams/dimlong)
fig = plt.figure()
for ii in range(len(fmaps)):
    ax = fig.add_subplot(dim1d, dimlong, ii+1)
    ax.imshow(fmaps[ii], cmap='viridis')
    ax.set_title('beam {:d}\nc: {:0.0f} g: {:0.0f}'.format(beam_indices[ii], couch_gantry[ii][0]*180/math.pi, couch_gantry[ii][1]*180/math.pi))
    ax.set_xticks([])
    ax.set_yticks([])
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.tight_layout()
plt.show()
