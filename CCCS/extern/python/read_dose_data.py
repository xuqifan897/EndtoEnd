""" Python data parser to construct the a sparse A matrix from the hdf5 output

This function provides a parser for the HDF5 data structure created by the beamlet-based dose calcluation
    program. Notably, this function accepts as input the path to the hdf5 file and returns a sparsely
    constructed dose coefficient ("A") matrix which has its columns ordered in increasing beam number then
    increasing beamlet number, and its rows ordered in the linear index of the target voxel. This function
    also returns the ordered list of beam number and beamlet number associated with each column such that
    the correspondence is not lost when its time to map fluence weights back to each selected beam.

Authors: Ryan Neph & Cheng Ouyang
Revisions:
    - 09 Aug.  2017: Initial
    - 26 Feb.  2018: support for new output format that supports dataset linking (in dosecalc v0.6.3+)
"""

import h5py
import sys
import os
import argparse
import numpy as np
import scipy.sparse as sps
import scipy.io as sio
import math
import gc
import psutil

def memory_usage():
    process = psutil.Process(os.getpid())
    return str(process.memory_info().rss)

def sparse_matrix_size_bytes(sparsemat):
    return sum([sparsemat.data.nbytes, sparsemat.indices.nbytes, sparsemat.indptr.nbytes])

def data_mat_parser(input_fid, verbose=False, only_metadata=False):
    """ parse the beamlet information and generate a big scipy sparse A matrix"""
    try:
        full_file = h5py.File(input_fid, 'r')["/"]
    except FileNotFoundError as e:
        print(e, e.__traceback__)
        raise FileNotFoundError("ERROR: {!s} does not exist".format(input_fid))

    meta_dict = {
        "calc_specs" : {},
        "beams" : {}
    }
    # a matrix dictionary
    # beam_idx: beamlet_idx: matrix. sort and combine them at the end
    Amat_dict = {} # beam id -> beamlet id -> vector
    num_beamlet = 0
    num_nonzero = 0 # number of non-zero entries for the whole dataset

    # parse calculation metadata
    calc_meta_group = full_file['calc_specs']
    for attr, value in calc_meta_group.attrs.items():
        meta_dict["calc_specs"][attr] = value
    if verbose: print("calc_specs has been parsed!")

    # iteration through beams
    for beam_name, beam_group in full_file['beams']['metadata'].items():
        if verbose: print(beam_name)

        # first store the metadata for a beam
        if beam_name not in meta_dict["beams"].keys():
            beam_idx = int(beam_name.split("_")[1])
            meta_dict["beams"][beam_idx] = {}
            Amat_dict[beam_idx] = {}
            for attr_name, attr_value in beam_group.attrs.items():
                # beam specs is compound datatype
                if "beam_specs" in attr_name:
                    meta_dict["beams"][beam_idx]["beam_specs_dict"] = dict(zip(attr_value.dtype.names, attr_value)) #convert to dict
                    meta_dict["beams"][beam_idx][attr_name] = attr_value # retain compound type as well
                else:
                    meta_dict["beams"][beam_idx][attr_name] = attr_value

        # now deal with beamlets
    for beam_name, beam_group in full_file['beams']['data'].items():
        beam_idx = int(beam_name.split("_")[1])
        Amat_dict[beam_idx]["beamlets"] = {}
        meta_dict["beams"][beam_idx]["beamlets"] = {}
        for blt_name, blt_group in beam_group.items():
            num_beamlet += 1
            blt_idx = int(blt_name.split("_")[1])
            Amat_dict[beam_idx]["beamlets"][blt_idx] = {}
            # write dose vector
            Amat_dict[beam_idx]["beamlets"][blt_idx]["coeffs"] = blt_group["coeffs"]
            Amat_dict[beam_idx]["beamlets"][blt_idx]["lindex"] = blt_group["lindex"]
            Amat_dict[beam_idx]["beamlets"][blt_idx]["num_dose_vox"] = blt_group["lindex"].shape[0]
            meta_dict["beams"][beam_idx]["beamlets"][blt_idx] = {}
            num_nonzero += blt_group["coeffs"].shape[0]
            # write beamlet metadata
            for blt_attr_name, blt_attr_value in blt_group.attrs.items():
                meta_dict["beams"][beam_idx]["beamlets"][blt_idx][blt_attr_name] = blt_attr_value

    # now construct the sparse matrix
    # get the size of the whole a matrix
    reduced = ("roi_order" in meta_dict["calc_specs"])
    if (reduced):
        if (verbose): print("  *Reduced matrix data detected*")
        num_vox = int( np.sum( meta_dict["calc_specs"]["row_block_capacities"] ) )
    else:
        num_vox = int( np.prod( meta_dict["calc_specs"]["full_dicom_size"] ) )

    # build col_idx array
    col_idx = [] # column index of beams and beamlets to be stored later for column indexing
    for beam_idx, beam_info in sorted( Amat_dict.items(), key = lambda t : int(t[0]) ) :
        for beamlet_idx in sorted(beam_info["beamlets"].keys(), key = lambda t: int(t) ) :
            col_idx.append( (beam_idx, beamlet_idx) )

    if only_metadata:
        return meta_dict, col_idx

    rows = np.zeros(num_nonzero, dtype = np.uint64)
    cols = np.zeros(num_nonzero, dtype = np.uint64)
    elems = np.zeros(num_nonzero, dtype = np.float32)
    beamlet_counter = 0
    # sorting
    vec_idx = 0 # starting index for recording current vector builing point
    for beam_idx, beam_info in sorted( Amat_dict.items(), key = lambda t : int(t[0]) ) :
        for beamlet_idx, beamlet_info in sorted(beam_info["beamlets"].items(), key = lambda t: int(t[0]) ) :
            if verbose: print("now visiting beam %s beamlet %s ..." %(beam_idx,beamlet_idx))
#            rows = np.append(rows, beamlet_info["lindex"])
#            elems = np.append(elems, beamlet_info["coeffs"])
#            cols = np.append(cols, np.full(beamlet_info["coeffs"].shape, beamlet_counter))
            rows[vec_idx : vec_idx + beamlet_info["num_dose_vox"]] = np.array(beamlet_info["lindex"],dtype = np.uint64)
            elems[vec_idx : vec_idx + beamlet_info["num_dose_vox"]] = np.array(beamlet_info["coeffs"], dtype = np.float32)
            cols[vec_idx : vec_idx + beamlet_info["num_dose_vox"]] = np.full(beamlet_info["num_dose_vox"], beamlet_counter, dtype = np.uint64)
            beamlet_counter += 1
            if verbose: print(vec_idx)
            vec_idx += beamlet_info["num_dose_vox"]

    if (rows.shape != cols.shape) or (rows.shape != elems.shape):
        print("number of elements and coordinates are not consistant. Exited with error! ")
        exit()
    # now make the sparse matrix
    if verbose:
        print("before recycling, memory usage: " + memory_usage())
    # clear unnecessary memory
    Amat_dict.clear()
    del Amat_dict
    del full_file
    gc.enable()
    gc.collect()
    if verbose:
        print("now memory usage: " + memory_usage())
        print("now writing the sparse matrix ... ")
    Amat = sps.csc_matrix( (elems, (rows, cols)), shape = (num_vox, num_beamlet), dtype = np.float32)
    if verbose:
        print("sparse matrix has been constructed.")
        bytes = sparse_matrix_size_bytes(Amat)
        if bytes < 1024:
            storage_size = bytes;
            storage_units = 'B';
        elif bytes < 1048576:
            storage_size = bytes/1024;
            storage_units = 'KB';
        elif bytes < 1073741824:
            storage_size = bytes/1048576;
            storage_units = 'MB';
        else:
            storage_size = bytes/1073741824;
            storage_units = 'GB';
        print('Total memory usage for sparse matrix: {:f} {:s}'.format(storage_size, storage_units))
    return Amat, meta_dict, col_idx

def save_npz(fid_prefix, Amat, meta = None, col_idx = None):
    """ save the sparse array to npz file  """
    if meta is not None:
        np.savez(fid_prefix + "meta.npz", meta)
        print("metadata has been saved as fid_prefix" + "meta.npz")

    if col_idx is not None:
        fio = open(fid_prefix + "col_idx.csv",'w')
        fio.write( "\n".join("%s %s" %x for x in col_idx ) )
        fio.close()
        print("column index has been saved as fid_prefix" + "col_idx.csv")

    np.savez(fid_prefix + "A_matrix.npz", Amat)
    print("A matrix has been saved as %s A_matrix.npz"%fid_prefix)

def main(input_fid, save_prefix, verbose=False):
    Amat, metadata, col_idx = data_mat_parser(input_fid, verbose)
    save_npz(save_prefix, Amat, metadata, col_idx)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Utility for reading dosecalc data (h5) into a python sparse matrix")
    parser.add_argument('infile', type=str, help="path specifying a valid dosecalc h5 result file")
    parser.add_argument('out', '-o', type=str, help="prefix or absolute path to output of npz sparse matrix", default=None)
    parser.add_argument('--verbose', '-v', action='count', help="show verbose output")
    args = parser.parse_args()
    main(args.infile, args.out, args.verbose)
