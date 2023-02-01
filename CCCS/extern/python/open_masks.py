import sys, os
import numpy as np
import h5py

def pad_array(arr, size, crop_start, crop_size):
    """fill an array of 'size' with zeros then inject 'arr' into it according to 'crop_start' and 'crop_size'"""
    padarr = np.zeros(size, dtype=np.uint8)
    slices = tuple([slice(crop_start[ii], crop_start[ii]+crop_size[ii]) for ii in range(arr.ndim)])
    padarr[slices] = arr
    return padarr

def open_masks(h5file):
    masks = {}
    with h5py.File(h5file) as h5fd:
        for group in h5fd.values():
            props = {k: v for k,v in group['ArrayProps'].attrs.items()}
            arr = group['mask'][:].reshape(props['crop_size'][::-1])
            padmask = pad_array(arr, size=props['size'][::-1], crop_start=props['crop_start'][::-1], crop_size=props['crop_size'][::-1])
            masks[group.attrs['name'].decode('UTF-8')] = padmask
    return masks


if __name__ == '__main__':
    outdir = 'masks'
    os.makedirs(outdir, exist_ok=True)
    for name, mask in open_masks(sys.argv[1]).items():
        np.save(os.path.join(outdir, name+'.npy'), mask)
    print('mask files saved to "{}"'.format(outdir))

