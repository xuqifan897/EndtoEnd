#include "volume.h"

std::ostream& operator<<(std::ostream& out, const FrameOfReference& frame) {
  out << "Frame:" << std::endl
      << "  Size:    ("<<frame.size.x    <<", "<< frame.size.y    <<", "<< frame.size.z <<")"<< std::endl
      << "  Start:   ("<<frame.start.x   <<", "<< frame.start.y   <<", "<< frame.start.z <<")"<< std::endl
      << "  Spacing: ("<<frame.spacing.x <<", "<< frame.spacing.y <<", "<< frame.spacing.z <<")"<< std::endl;
  return out;
}


template<>
H5::DataSet Volume<uint8_t>::_createTypedH5Dataset(H5::Group& h5group, const uint8_t* data, const FrameOfReference& frame, const std::string& dset_name) {
    return typedH5DatasetFactory(h5group, data, frame, dset_name, H5::PredType::STD_U8LE, H5::PredType::NATIVE_UINT8);
}
template<>
void Volume<uint8_t>::_readTypedH5Dataset(std::vector<uint8_t>& vol, uint3& size, H5::DataSet& h5dset) {
    typedH5DatasetLoader(vol, size, h5dset, H5::PredType::NATIVE_UINT8);
}

template<>
H5::DataSet Volume<short>::_createTypedH5Dataset(H5::Group& h5group, const short* data, const FrameOfReference& frame, const std::string& dset_name) {
    return typedH5DatasetFactory(h5group, data, frame, dset_name, H5::PredType::STD_I16LE, H5::PredType::NATIVE_SHORT);
}
template<>
void Volume<short>::_readTypedH5Dataset(std::vector<short>& vol, uint3& size, H5::DataSet& h5dset) {
    typedH5DatasetLoader(vol, size, h5dset, H5::PredType::NATIVE_SHORT);
}

template<>
H5::DataSet Volume<char>::_createTypedH5Dataset(H5::Group& h5group, const char* data, const FrameOfReference& frame, const std::string& dset_name) {
    return typedH5DatasetFactory(h5group, data, frame, dset_name, H5::PredType::STD_I8LE, H5::PredType::NATIVE_CHAR);
}
template<>
void Volume<char>::_readTypedH5Dataset(std::vector<char>& vol, uint3& size, H5::DataSet& h5dset) {
    typedH5DatasetLoader(vol, size, h5dset, H5::PredType::NATIVE_CHAR);
}

template<>
H5::DataSet Volume<float>::_createTypedH5Dataset(H5::Group& h5group, const float* data, const FrameOfReference& frame, const std::string& dset_name) {
    return typedH5DatasetFactory(h5group, data, frame, dset_name, H5::PredType::IEEE_F32LE, H5::PredType::NATIVE_FLOAT);
}
template<>
void Volume<float>::_readTypedH5Dataset(std::vector<float>& vol, uint3& size, H5::DataSet& h5dset) {
    typedH5DatasetLoader(vol, size, h5dset, H5::PredType::NATIVE_FLOAT);
}
