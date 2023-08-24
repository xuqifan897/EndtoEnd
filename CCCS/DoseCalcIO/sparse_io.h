#ifndef __SPARSE_IO_H__
#define __SPARSE_IO_H__

#include <helper_cuda.h>
#include <helper_math.h>
#include <iostream>
#include <vector>
#include <string>

#include "./io_data_structs.h"
#include "./roi.h"

#include "H5Cpp.h" // Namespace H5::

//////////////// Helpers ////////////////
// low-level function for construction sparse data vectors
const float _default_sparsity_threshold = 1e-12;
const int _default_compress_lvl = -1;
const int _default_chunksize = 50;
void _sparsify(SparseData& sparsedata, const float* array, ArrayProps& props, float thresh=_default_sparsity_threshold);
void _reduce_by_roi(SparseData& reduceddata, SparseData& sparsedata, ROIMaskList& roi_list, ArrayProps& props);
void _sparsify_and_reduce(SparseData& data, const float* array, ArrayProps& props, ROIMaskList* roi_list=nullptr, float thresh=_default_sparsity_threshold);
int _write_beamlet_to_hdf5(H5::Group&, const SparseData&, HEADER_BEAMLET&, int compress_lvl=_default_compress_lvl, uint64_t chunksize=_default_chunksize);
int _write_beamlet_metadata(H5::Group&, const HEADER_BEAMLET&, const SparseData&);
int _write_beam_metadata(H5::Group&, const HEADER_BEAM&);
int _write_patient_metadata(H5::Group&, const HEADER_PATIENT&);
/////////////////////////////////////////

//////////////// Write-to-file methods ////////////////

// write to plain ascii text file
int sparse_to_text(const std::string& filename, const SparseData& data, const void* _);

// write to binary file
int sparse_to_binary(const std::string& filename, const SparseData& data, const void* _);

// // write to hdf5 file
int sparse_to_hdf5(const std::string& filename, const SparseData& data, const void* header);
///////////////////////////////////////////////////////


// creates DOK sparse vector from a dense linearized 3D array of values strictly greater than "thresh" using
// the write-to-file method defined by "tofile" function pointer
int write_sparse_beamlet_to_file(const std::string& filename, const float* array, ArrayProps& props, int(*tofile)(const std::string&, const SparseData&, const void* header)=sparse_to_hdf5, void* header=nullptr, ROIMaskList* roi_list=nullptr, float thresh=_default_sparsity_threshold);
int write_beam_metadata(const std::string& filename, const HEADER_BEAM&);
int write_patient_metadata(const std::string&, const HEADER_PATIENT&);
// int test_read_binary(const std::string& filename, KeyValPairs& in_sparse);

#endif // __SPARSE_IO_H__
