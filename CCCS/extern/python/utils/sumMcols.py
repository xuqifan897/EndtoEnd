#!/usr/bin/env python
"""
Reconstructs final dose distribution from M-matrix by summing across the columns and outputting in .raw
format (linear packed bytes - float values) and in standard rs4pi compliant "Final Dose" file (h5)
"""
import sys
import os.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))
import argparse
import struct
import numpy as np
import read_dose_data
import h5py


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Beam angle selection and fluence map optimization algorithm for use in radiation therapy treatment planning',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('dosefile', type=str, help='path to dose calculation output file containing beamlet coefficients data')
    parser.add_argument('--out', '-o', nargs='?', default="fmo_results", type=str, help='results output file name')
    parser.add_argument('--verbose', '-v', nargs='?', type=int, default=1, const=1, help='enable verbose output')
    args = parser.parse_args()
    dosefile = args.dosefile
    outfile = os.path.splitext(args.out)[0] + '.raw'
    verbose = args.verbose

    if verbose: print(f'Reading dose data from \"{dosefile}\"')
    result = read_dose_data.data_mat_parser(dosefile)

    if isinstance(result, tuple):
        (A, meta_dict, col_idx) = result
    else:
        raise Exception("There was a problem reading the dose data file")

    # sum rows of A - get total dose
    fulldose = A.sum(1)

    # save to .raw
    doseshape = meta_dict['calc_specs']["full_dicom_size"]
    if verbose: print('Writing result to \"{}\" with size ({},{},{})'.format(outfile, *doseshape))
    with open(outfile, 'wb') as f:
        f.write(struct.pack('{}f'.format(A.shape[0]), *fulldose[:,0]))

    # write to hdf5
    h5_outfile = outfile.replace('.raw', '.h5')
    if verbose: print(f'Writing result to \"{h5_outfile}\"')
    with h5py.File(h5_outfile, 'w') as f:
        g = f.create_group('filetype')
        g.attrs.create('ftmagic', int('0x2A', 16), dtype=np.uint8)
        g.attrs.create('ftversionmajor', 0, dtype=np.uint8)
        g.attrs.create('ftversionminor', 8, dtype=np.uint8)

        d = f.create_dataset('dose', data=np.asarray(fulldose).reshape(doseshape), dtype=np.dtype('<f4'))
        d.attrs.create('dicom_start_cm', meta_dict['calc_specs']['dicom_start_cm'], dtype=np.dtype('<f4'))
        d.attrs.create('voxel_size_cm', meta_dict['calc_specs']['voxel_size_cm'], dtype=np.dtype('<f4'))
        d.attrs.create('mem_layout', "row_major_zyx", dtype="S14")
