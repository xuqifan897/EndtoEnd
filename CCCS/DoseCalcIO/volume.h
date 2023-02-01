#ifndef __VOLUME_H__
#define __VOLUME_H__

#include <type_traits>
#include <iostream>
#include <exception>
#include <memory>
#include <vector>
#include <vector_types.h>
#include <vector_functions.h>
#include "H5Cpp.h"

#include "./binary_io.h"
#include "./macros.h"

typedef unsigned int uint;
int _write_file_version(H5::Group&, uint, uint, uint);

// Data volume properties
struct FrameOfReference {
    uint3  size;    // data array dimensions
    float3 start;   // scanner coords of first voxel in data array [unit: mm]
    float3 spacing; // voxel size in [unit: mm]
    uint nvoxels() const { return size.x*size.y*size.z; }
    float3 end() const { return make_float3(start.x + spacing.x*size.x, start.y + spacing.y*size.y, start.z + spacing.z*size.z); }
};
std::ostream& operator<<(std::ostream& out, const FrameOfReference& frame);

template <typename T>
class Volume {
    public:
        Volume() {}
        // copy constructor
        Volume(const FrameOfReference& frame) :
            start{frame.start}, size{frame.size}, voxsize{frame.spacing}, _vect{std::vector<T>(frame.nvoxels())}
            {}
        Volume(T* data, const FrameOfReference& frame) :
            start{frame.start}, size{frame.size}, voxsize{frame.spacing}, _vect{std::vector<T>(data, data+frame.nvoxels())}
            {}

        float3 start;                // coords of first voxel in scanner coord system (GCS)
        float3 voxsize;              // voxel dimensions in mm
        uint3  size;                 // volume shape (units of integer voxels)
        std::vector<T> _vect;        // linear array flattened in C-Major order (Depth->Row->Column)

        T* data()  {return _vect.data();};
        const T* data() const {return _vect.data();};
        // copy count elements from ptr to vector
        void set_data(T* ptr, uint count) { _vect = std::vector<T>(ptr, ptr+count); }
        // initialize to size=count with zeros
        void set_data(uint count) { _vect = std::vector<T>(count); }

        inline uint nvoxels() { return size.x * size.y * size.z; }
        inline uint mem_size() { return nvoxels() * sizeof(T); }

        const FrameOfReference get_frame() const {
          return FrameOfReference{size, start, voxsize};
        }

        // OP OVERLOADS
        T& operator [](uint idx) { return _vect[idx]; }
        const T& operator [](uint idx) const { return _vect[idx]; }
        T& operator [](int idx) { return _vect[idx]; }
        const T& operator [](int idx) const { return _vect[idx]; }
        T& operator [](uint3 coords) { return at(coords); }
        const T& operator [](uint3 coords) const { return at(coords); }
        T& operator [](int3 coords) { return at(coords); }
        const T& operator [](int3 coords) const { return at(coords); }
        T& at(uint idx) { return _vect.at(idx); }
        const T& at(uint idx) const { return _vect.at(idx); }
        T& at(uint x, uint y, uint z) { return _vect.at(x + size.x*(y + size.y*z)); }
        const T& at(uint x, uint y, uint z) const { return _vect.at(x + size.x*(y + size.y*z)); }
        T& at(uint3 coords) { return at(coords.x, coords.y, coords.z); }
        const T& at(uint3 coords) const { return at(coords.x, coords.y, coords.z); }
        T& at(int3 coords) { return at(coords.x, coords.y, coords.z); }
        const T& at(int3 coords) const { return at(coords.x, coords.y, coords.z); }
        bool check_bounds(uint idx) { return (idx < nvoxels()); }
        bool check_bounds(uint x, uint y, uint z) { return (x < size.x && y < size.y && z < size.z); }
        bool check_bounds(uint3 coords) { return (coords.x < size.x && coords.y < size.y && coords.z < size.z); }
        bool check_bounds(int3 coords) { return (coords.x>=0 && coords.y>=0 && coords.z>=0 &&
                coords.x < size.x && coords.y < size.y && coords.z < size.z); }

        // FILE I/O
        static int _readFromHDF5(Volume<T>& vol, H5::Group& h5group, const std::string& dset_name="data") {
            // read data array
            auto dset = h5group.openDataSet(dset_name);
            _readTypedH5Dataset(vol._vect, vol.size, dset);

            { // read start coords
                auto att = dset.openAttribute("dicom_start_cm");
                float temp[3];
                att.read(H5::PredType::NATIVE_FLOAT, &temp);
                ARR3VECT(vol.start, temp);
            } { // read voxelsize
                auto att = dset.openAttribute("voxel_size_cm");
                float temp[3];
                att.read(H5::PredType::NATIVE_FLOAT, &temp);
                ARR3VECT(vol.voxsize, temp);
            }
            return true;
        }
        static int readFromFile(Volume<T>& vol, const std::string& infile, const std::string& dset_name="data", int verbose=false) {
            try {
                vol = Volume<T>{};
                auto h5file = H5::H5File(infile, H5F_ACC_RDONLY);
                H5::Group rootgroup = h5file.openGroup("/");
                if (!Volume<T>::_readFromHDF5(vol, rootgroup, dset_name)) {
                    if (verbose) { std::cerr << "Failed to read Volume from \""<<infile<<"\""<<std::endl; }
                    return false;
                }
            } catch (H5::FileIException &file_exists_error) {
                if (verbose) { std::cerr << "Failed to read Volume from \""<<infile<<"\""<<std::endl; }
                return false;
            }
            return true;
        }

        static int _writeToHDF5(H5::Group& h5group, const T* data, const FrameOfReference& frame, const std::string& dset_name) {

            // write dose float-volume to dataset
            H5::DataSet dset = _createTypedH5Dataset(h5group, data, frame, dset_name);

            // add attributes
            {
                hsize_t dims[] = {3};
                auto dspace = H5::DataSpace(1, dims);
                auto att = dset.createAttribute("dicom_start_cm", H5::PredType::IEEE_F32LE, dspace);
                float temp[3];
                VECT3ARR(temp, frame.start);
                att.write(H5::PredType::NATIVE_FLOAT, &temp);
            } {
                hsize_t dims[] = {3};
                auto dspace = H5::DataSpace(1, dims);
                auto att = dset.createAttribute("voxel_size_cm", H5::PredType::IEEE_F32LE, dspace);
                float temp[3];
                VECT3ARR(temp, frame.spacing);
                att.write(H5::PredType::NATIVE_FLOAT, &temp);
            }
            return true;
        }
        int writeToFile(const std::string& outfile, const std::string& dset_name="volume") {
            FrameOfReference frame = { size, start, voxsize };
            Volume<T>::writeToFile(outfile, _vect.data(), frame, dset_name);
        }
        static int writeToFile(const std::string& outfile, T* data, FrameOfReference& frame, const std::string& dset_name="volume") {
            auto h5file = H5::H5File(outfile, H5F_ACC_TRUNC);
            auto rootgroup = h5file.openGroup("/");
            _write_file_version(rootgroup, FTMAGIC, FTVERSIONMAJOR, FTVERSIONMINOR);
            return _writeToHDF5(rootgroup, data, frame, dset_name);
        }
        int writeToRawFile(const std::string& outfile, bool verbose=false) {
            return write_binary_data<T>(this->data(), size, outfile.c_str(), verbose);
        }

    protected:
        static const uint FTMAGIC = 0x2A;
        static const uint FTVERSIONMAJOR = 1;
        static const uint FTVERSIONMINOR = 0;

        static H5::DataSet typedH5DatasetFactory(H5::Group& h5group, const T* data, const FrameOfReference& frame, const std::string& dset_name, const H5::PredType& fileDataType, const H5::PredType& memDataType) {
            hsize_t dims_mem[] = {frame.nvoxels()};
            hsize_t dims_file[] = {frame.size.z, frame.size.y, frame.size.x};
            auto dspace_mem = H5::DataSpace(1, dims_mem);
            auto dspace_file = H5::DataSpace(3, dims_file);
            H5::DataSet dset = h5group.createDataSet(dset_name, fileDataType, dspace_file);
            dset.write(data, memDataType, dspace_mem, dspace_file);
            return dset;
        }

        static void typedH5DatasetLoader(std::vector<T>& outvect, uint3& outsize, H5::DataSet& h5dset, const H5::PredType& memDataType) {
            auto file_space = h5dset.getSpace();
            hsize_t N = file_space.getSimpleExtentNpoints();
            std::unique_ptr<T[]> temp(new T[N]);
            hsize_t dims[] = {N};
            auto mem_space = H5::DataSpace(1, dims);
            h5dset.read(temp.get(), memDataType, mem_space, file_space);
            outvect.clear();
            outvect.assign(&temp[0], &temp[N]);
            hsize_t file_dims[3];
            file_space.getSimpleExtentDims(file_dims);
            outsize.x = file_dims[2];
            outsize.y = file_dims[1];
            outsize.z = file_dims[0];
        }

        // These need to be specialized for each templated Volume type that will be used
        static H5::DataSet _createTypedH5Dataset(H5::Group& h5group, const T* data, const FrameOfReference& frame, const std::string& dset_name);
        static void _readTypedH5Dataset(std::vector<T>&, uint3&, H5::DataSet&);

};
typedef Volume<float>   FloatVolume;
typedef Volume<short>   ShortVolume;
typedef Volume<char>    CharVolume;
typedef Volume<uint8_t> BoolVolume;

template <typename T>
std::ostream& operator<<(std::ostream& out, const Volume<T>& vol) {
  out << vol.get_frame();
  return out;
}

////////////////

#endif // __VOLUME_H__
