#ifndef __READ_SPECFILE_H__
#define __READ_SPECFILE_H__

// save on include by declaring struct name only
struct MONO_KERNELS;

int read_spectrum_file(MONO_KERNELS *mono, bool verbose);

#endif // __READ_SPECFILE_H__
