#ifndef __FMAPS_IO_H__
#define __FMAPS_IO_H__

#include <vector>
#include <string>

#include "./io_data_structs.h"
#include "./beam.h"
#include "./volume.h"

#include "H5Cpp.h" // Namespace H5::

// Methods for reading FMO result file
int _read_patient_metadata(H5::Group& grp, HEADER_PATIENT& patient_header);
int _read_beams(H5::Group& h5group, std::vector<HEADER_BEAM>& beam_headers);
int read_fmaps(const std::string& filename, HEADER_PATIENT* patient_header, std::vector<HEADER_BEAM>& beam_headers );
int read_fmaps(const std::string& filename, HEADER_PATIENT* patient_header, std::vector<BEAM>& beams );

int _write_file_version(H5::Group&, uint ftmagic, uint ftversionmajor, uint ftversionminor);

#endif //__FMAPS_IO_H__
