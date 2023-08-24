"""fmo_io.py

Defines data I/O for use by treatment planning software suite
"""
import sys
import os.path
import argparse
import re
import h5py
import numpy as np
import scipy.io
import math
from read_dose_data import data_mat_parser

def create_fluence_map(wts, beamlet_idx, fmap_shape):
    """take co-registerd lists of beamlet indices and accompanying fluence intensities and construct a
    dense fluence map

    Inputs:
        wts (list[N]):         fluence intensities
        beamlet_idx (list[N]): beamlet row-major index into fmap
        fmap_shape (tuple[2]): dims of fluence map
    Output:
        (np.ndarray[r,c]):     2d array representing dense fluence map
    """
    if not len(wts) == len(beamlet_idx):
        raise RuntimeError("wts and beamlet_idx must have equal lengths")

    fmap = np.zeros((fmap_shape), dtype=np.float64)
    for ii in range(len(wts)):
        r = int(beamlet_idx[ii]) // fmap_shape[1]
        c = int(beamlet_idx[ii]) % fmap_shape[1]
        fmap[r, c] = wts[ii]
    import matplotlib.pyplot as plt
    plt.imshow(fmap)
    plt.show()
    return fmap


def write_fmo_results(outfile, fm_wts, col_idx, meta, selected_beams=None, overwrite=False):
    """writes result of FMO to h5 file with metadata necessary for running final dose calculation
    on selected beams and their optimized fluence maps

    Inputs:
        outfile (str):               output filename
        fm_wts (np.ndarray[N,1]):    concatenated vector of fluence intensities for all beamlets
        col_idx (list(tuple[2])[N]): tuple_i[0]: beam_idx; tuple[1]: beamlet_idx
        meta (dict{...}):            metadata including calc_specs and beamt_meta
    Optional Inputs:
        selected_beams (list[P]):     list of beam_indices for selected beams
        overwrite (bool):            overwrite existing outfile?
    """
    # clean filename
    outfile = os.path.splitext(outfile)[0] + '.h5'

    sparse_wts = bool(col_idx)

    # open output file
    with h5py.File(outfile, 'w' if overwrite else 'w-') as h5file:
        # write calculation metadata
        calc_meta_grp = h5file.create_group('calc_specs')
        for name, att in meta['calc_specs'].items():
            calc_meta_grp.attrs.create(name, att)

        if sparse_wts:
            # fmap wts is a concatenated sparse vector registered to col_idx (beam#, beamlet#)
            # unpack
            unpack = tuple(zip(*col_idx))
            beam_indices = list(unpack[0])
            beamlet_indices = list(unpack[1])
            # get ordered set of beam indices as list - match order of columns in A-matrix
            ordered_beam_indices = []
            last = None
            for idx in beam_indices:
                if idx != last: ordered_beam_indices.append(idx)
                last = idx
        else:
            ordered_beam_indices = selected_beams

        # store beam metadata
        all_beams_grp = h5file.create_group("beams")
        beam_meta_grp = all_beams_grp.create_group("metadata")
        ptr = 0
        nbeams = 0
        for b, beam_idx in enumerate(ordered_beam_indices):
            beam_meta = meta["beams"][beam_idx]

            if sparse_wts:
                nbeamlets = beam_meta['N_beamlets']

                # only include selected_beams if provided
                if (selected_beams and not int(beam_idx) in selected_beams):
                    ptr += nbeamlets
                    continue

                beam_grp = beam_meta_grp.create_group(str(beam_idx))

                # optimized fmap intensities
                wts = fm_wts[ptr:ptr+nbeamlets]
                beamlet_idx = beamlet_indices[ptr:ptr+nbeamlets]
                fmap = create_fluence_map(wts, beamlet_idx, tuple(beam_meta['beam_specs']['fmap_dims']))
                beam_grp.attrs.create('fmap_weights', fmap, dtype=np.dtype('<f4'))
                ptr += nbeamlets
                nbeams += 1

            else:
                beam_grp = beam_meta_grp.create_group(str(beam_idx))
                fmap = fm_wts[b]
                beam_grp.attrs.create('fmap_weights', fmap, dtype=np.dtype('<f4'))
                nbeams += 1

            # general beam meta
            for name, att in beam_meta.items():
                # skip some items not meant for h5
                if name == 'fmap_weights' or name == 'beam_specs_dict' or name == 'beamlets': continue
                beam_grp.attrs.create(name, att)

        # write true N_beams to "calc_specs"
        calc_meta_grp.attrs.modify('N_beams', nbeams)

def check_file(f):
    return os.path.exists(f) and os.path.isfile(f)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='converts fluence maps in .npz or .mat format into --fmaps= compatible inputs for dosecalc-preprocess')
    parser.add_argument('fluence_file', type=str, help='path to .npz or .mat file containing a fluence map or collection of fluence maps in their full size')
    parser.add_argument('dose_file', type=str, help='path to dosecalc output (h5) file')
    parser.add_argument('--out', '-o', type=str, help='output prefix or path to fmaps compatible file', default='fmaps.h5')
    args = parser.parse_args()

    if not check_file(args.fluence_file):
        raise FileNotFoundError("{} was not a valid fluence file".format(args.fluence_file))
    if not check_file(args.dose_file):
        raise FileNotFoundError("{} was not a valid dose output file".format(args.dose_file))

    fm_wts = []
    selected_beams = []
    meta, _ = data_mat_parser(args.dose_file, only_metadata=True)

    # handle different types of fluence map input
    ext = os.path.splitext(args.fluence_file)[1]
    if ext == '.mat':
        try:
            d = scipy.io.loadmat(args.fluence_file, appendmat=False, chars_as_strings=True)
        except NotImplementedError as err:
            d = {}
            with h5py.File(args.fluence_file, 'r') as f:
                for k, v in f.items(): d[k] = np.array(v)
        for k, v in d.items():
            m = re.search(r'(?:beam_?(\d+))', k)
            if m:
                if v.ndim != 2: raise ValueError('fluence map input must have 2 dims, not {}'.format(v.ndim))
                beam_idx = int(m.group(1))
                selected_beams.append(beam_idx)
                fm_wts.append(v)

                import matplotlib.pyplot as plt
                plt.imshow(v)
                plt.show()
        if not selected_beams or not fm_wts:
            raise RuntimeError("No valid fluence weight arrays were found in fluence_file")

    if ext == '.npz':
        raise NotImplementedError(".npz inputs are not yet supported")

    write_fmo_results(args.out, fm_wts, None, meta, selected_beams=selected_beams, overwrite=True)
