#include "roi.h"

#include <cstdio>
#include <iostream>
#include <cmath>
#include <list>
#include <limits>

#include "helper_cuda.h"
#include "helper_math.h"

/* helpers at a higher level of abstraction to low-level rtstruct lib for matching names, constructing encapsulated objects... */

///// DATA STRUCTURES /////

BaseROIMask::~BaseROIMask() {}
int BaseROIMask::writeToFile(std::string fname, bool verbose) {
    // open hdf5 file
    H5::H5File h5file = H5::H5File(fname, H5F_ACC_TRUNC);
    H5::Group rootgroup = h5file.openGroup("/");
    if (!_writeToHDF5(rootgroup)) {
        if (verbose){ std::cout << "Failed to write ROIMask to \""<<fname<<"\""<<std::endl; }
        return false;
    }
    return true;
}

uint8_t DenseROIMask::operator[](uint64_t idx) {
    if (idx >= nvoxels()) { return false; }
    return mask[idx];
}
uint DenseROIMask::nones() {
    // cache if not done yet, otherwise get from cache
    if (_nones_cached == 0) {
        uint n = 0;
        for (uint i=0; i<nvoxels(); ++i) { if (mask[i] == 1) { ++n; } }
        _nones_cached = n;
    }
    return _nones_cached;
}
// write/read HDF5
int DenseROIMask::_writeToHDF5(H5::Group& h5group) {
    H5::DataSpace scalarspace;

    // create attribute - name
    {
        H5::StrType str_t{H5::PredType::C_S1, name.length()+1};
        auto att = h5group.createAttribute("name", str_t, scalarspace);
        att.write(str_t, name);
    }

    // write array props to group
    auto array_props_group = h5group.createGroup("ArrayProps");
    props._writeToHDF5(array_props_group);

    // create dataset - dense mask
    if (nvoxels() > 0) {
        hsize_t dims[] = { nvoxels() };
        H5::DataSpace simplespace(1, dims);
        // no bool in HDF5, use int 0,1 instead
        auto dset = h5group.createDataSet("mask", H5::PredType::STD_U8LE, simplespace);
        dset.write(mask.data(), H5::PredType::NATIVE_UINT8);
    }

    return true;
}
int DenseROIMask::_readFromHDF5(DenseROIMask& mask, H5::Group& h5group) {
    // read attribute - name
    {
        auto att = h5group.openAttribute("name");
        H5::DataType str_t = att.getDataType();
        H5std_string buf("");
        att.read(str_t, buf);
        mask.name = buf;
    }

    // read group - ArrayProps
    {
        H5::Group props_group = h5group.openGroup("ArrayProps");
        ArrayProps props {};
        ArrayProps::_readFromHDF5(props, props_group);
        mask.props = props;
    }

    // read dataset - dense mask
    {
        H5::DataSet dset = h5group.openDataSet("mask");
        H5::DataSpace dspace = dset.getSpace();
        hsize_t N;
        dspace.getSimpleExtentDims(&N, NULL);
        mask.mask = std::vector<uint8_t>(N);
        dset.read(mask.mask.data(), H5::PredType::NATIVE_UINT8, dspace, dspace);
    }

    return true;
}



std::vector<std::string> ROIMaskList::getROINames() {
    std::vector<std::string> names{};
    for (const auto& roi : _coll) {
        names.push_back( roi->name );
    }
    return names;
}
std::vector<uint64_t> ROIMaskList::getROICapacities() {
    std::vector<uint64_t> capacities{};
    for (const auto& roi : _coll) {
        capacities.push_back( roi->nones() );
    }
    return capacities;
}
// write/read HDF5
int ROIMaskList::writeToFile(std::string fname, bool verbose) {
    // open hdf5 file
    H5::H5File h5file = H5::H5File(fname, H5F_ACC_TRUNC);
    H5::Group rootgroup = h5file.openGroup("/");
    if (!_writeToHDF5(rootgroup)) {
        if (verbose){ std::cout << "Failed to write ROIMaskList to \""<<fname<<"\""<<std::endl; }
        return false;
    }
    return true;
}
int ROIMaskList::readFromFile(ROIMaskList& masklist, std::string fname, bool verbose) {
    H5::Exception::dontPrint();
    try {
        masklist = ROIMaskList();
        H5::H5File h5file = H5::H5File(fname, H5F_ACC_RDONLY);
        H5::Group rootgroup = h5file.openGroup("/");
        if (!ROIMaskList::_readFromHDF5(masklist, rootgroup)) {
            if (verbose){ std::cout << "Failed to read ROIMaskList from \""<<fname<<"\""<<std::endl; }
            return false;
        }
    }
    catch (H5::FileIException &file_exists_error) {
        if (verbose){ std::cout << "Failed to read ROIMaskList from \""<<fname<<"\""<<std::endl; }
        return false;
    }
    return true;
}
int ROIMaskList::_writeToHDF5(H5::Group& h5group) {
    H5::DataSpace scalarspace;

    // store each ROIMask to its own group
    uint index = 0;
    for (auto&& roi : _coll) {
        ++index;
        auto roi_group = h5group.createGroup(roi->name);
        if (!roi->_writeToHDF5(roi_group)) { return false; }

        // also store 1-based index for each indicating order
        H5::Attribute att = roi_group.createAttribute("index", H5::PredType::STD_U16LE, scalarspace);
        att.write(H5::PredType::NATIVE_UINT, &index);
    }
    return true;
}
int ROIMaskList::_readFromHDF5(ROIMaskList& masklist, H5::Group& h5group, bool verbose) {
    // read groupnames and indices
    struct IndexedString { uint16_t idx; std::string str; };
    struct OpData { OpData(H5::Group& g) : h5group{g} {}; std::list<IndexedString> groups={}; H5::Group& h5group; };
    OpData opdata(h5group);
    int iter_idx = 0; // iter_count is returned here
    h5group.iterateElems(".", &iter_idx,
        [](hid_t loc_id, const char* name, void* opdata) -> herr_t {
            // iterator body
            // construct an IndexedString for each group and add to "groups" list
            OpData* data = static_cast<OpData*>(opdata);
            H5::Group roi_group = data->h5group.openGroup(name);
            auto att = roi_group.openAttribute("index");
            uint16_t index; att.read(H5::PredType::NATIVE_UINT16, &index);
            data->groups.push_back( IndexedString{index, std::string(name)} );
            return 0;
        }, (void*)&opdata);

    // sort (index, groupname) list on index ascending
    opdata.groups.sort( [](IndexedString& a, IndexedString& b) -> bool {
            // true if a belongs before b
            return (a.idx <= b.idx);
            } );

    if (verbose) { std::cout << "Reading ROIMaskList ("<<iter_idx<<"):" << std::endl; }
    for (auto v : opdata.groups) {
        if (verbose) { std::cout << v.idx << ": " << v.str << std::endl; }
        // open group and load ROI data
        H5::Group roi_group = h5group.openGroup(v.str);
        std::shared_ptr<DenseROIMask> mask = std::make_shared<DenseROIMask>();
        if (!DenseROIMask::_readFromHDF5(*mask, roi_group)) {
            if (verbose){ std::cout << "Failed to read ROIMaskList from H5 Group: \""<<v.str<<"\""<<std::endl; }
            return false;
        }
        masklist._coll.emplace_back(mask);
    }
    return true;
}
///////////////////////////


///// FUNCTIONS /////
// TODO: Add case insensitivity option
int getROIIndex(RTStruct& rtstruct, const std::string& search, bool exact, bool verbose) {
    // find the ptv contour and load the data points
    if (verbose) {
        if (exact) { printf("Exact ROI matching requested for: \"%s\"\n", search.c_str()); }
        else { printf("Substring ROI matching requested for sub: \"%s\"\n", search.c_str()); }
    }

    for (unsigned int r=0; r<rtstruct.getNumberOfROIs(); r++) {
        if (exact && strcmp(rtstruct.getROIName(r), search.c_str()) == 0) {
            return r;
        } else if (!exact && strstr(rtstruct.getROIName(r), search.c_str()) != NULL) {
            return r;
        }
    }
    return -1; // failed to match
}


int loadStructureSet(StructureSet& roi, RTStruct& rtstruct, uint roi_idx, bool verbose) {
    //TODO: dont repeat if already loaded
    if ( !rtstruct.loadRTStructData( roi_idx, verbose)) {
        if (verbose) printf("Failed to load data for structure: \"%s\"\n", rtstruct.getROIName(roi_idx));
        return false;
    }
    // contour data are stored as a series of sub-contours
    roi.sub_cntr_count = rtstruct.getROISubCntrCount(roi_idx);
    // each sub-contour contains a list of points, organized by slice
    roi.sub_cntr_points_count.resize(roi.sub_cntr_count);
    // total number of points in the contour
    roi.total_points_count = rtstruct.getROITotalPointsCount(roi_idx);
    // allocate points array
    roi.points.resize(3*roi.total_points_count);

    uint points_copied = 0;
    for (uint s=0; s < roi.sub_cntr_count; s++ ) {
        roi.sub_cntr_points_count[s] = rtstruct.getROISubCntrPointCount(roi_idx,s);
        float *sub_cntr_points = rtstruct.getROISubCntrPoints( roi_idx, s );
        memcpy( roi.points.data() + 3 * points_copied, sub_cntr_points, 3 * roi.sub_cntr_points_count[s] * sizeof(float) );
        points_copied += roi.sub_cntr_points_count[s];
    }
    if (verbose) printf("\n %d ROI points copied. \n",points_copied);
    return true;
}

ArrayProps getROIExtents(const StructureSet& roi, const FrameOfReference& frame, bool verbose) {
    float3 min_coord{ std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(),  std::numeric_limits<float>::max()};
    float3 max_coord{std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min()};

    // slice iterator
    uint mark = 0;
    for (int c=0; c<(int)roi.sub_cntr_count; ++c) {
        // point iterator
        for (int i=0; i<(int)roi.sub_cntr_points_count[c]; ++i) {
            // get coord
            float x = roi.points[3*(i+mark)  ];
            float y = roi.points[3*(i+mark)+1];
            float z = roi.points[3*(i+mark)+2];
            // if (verbose) { printf("  c:%d, i:%d, x:%f, y:%f, z:%f\n", c, i, x, y, z); }

            // update limits
            if (x < min_coord.x) { min_coord.x = x; }
            if (y < min_coord.y) { min_coord.y = y; }
            if (z < min_coord.z) { min_coord.z = z; }
            if (x > max_coord.x) { max_coord.x = x; }
            if (y > max_coord.y) { max_coord.y = y; }
            if (z > max_coord.z) { max_coord.z = z; }
        }
        mark += roi.sub_cntr_points_count[c];
    }

    // convert to indices
    int3 start, end;
    start.x = floor( (0.1f*min_coord.x - frame.start.x) / (frame.spacing.x) );
    start.y = floor( (0.1f*min_coord.y - frame.start.y) / (frame.spacing.y) );
    start.z = floor( (0.1f*min_coord.z - frame.start.z) / (frame.spacing.z) );
    end.x   = ceil(  (0.1f*max_coord.x - frame.start.x) / (frame.spacing.x) );
    end.y   = ceil(  (0.1f*max_coord.y - frame.start.y) / (frame.spacing.y) );
    end.z   = ceil(  (0.1f*max_coord.z - frame.start.z) / (frame.spacing.z) );

    // bounds checking
    if (start.x < 0) { start.x = 0; }
    if (start.y < 0) { start.y = 0; }
    if (start.z < 0) { start.z = 0; }
    if (start.x > (int)frame.size.x) { start.x = frame.size.x; }
    if (start.y > (int)frame.size.y) { start.y = frame.size.y; }
    if (start.z > (int)frame.size.z) { start.z = frame.size.z; }
    if (end.x < 0) { end.x = 0; }
    if (end.y < 0) { end.y = 0; }
    if (end.z < 0) { end.z = 0; }
    if (end.x > (int)frame.size.x) { end.x = frame.size.x; }
    if (end.y > (int)frame.size.y) { end.y = frame.size.y; }
    if (end.z > (int)frame.size.z) { end.z = frame.size.z; }

    if (verbose) {
        printf("Frame.size: (%d, %d, %d)  ||  Frame.start: (%6.2f, %6.2f, %6.2f)  ||  Frame.spacing: (%6.2f, %6.2f, %6.2f)\n",
                frame.size.x, frame.size.y, frame.size.z,
                frame.start.x, frame.start.y, frame.start.z,
                frame.spacing.x, frame.spacing.y, frame.spacing.z
                );
        printf("Min Coords: (%6.2f, %6.2f, %6.2f)  ||  Max Coords: (%6.2f, %6.2f, %6.2f)\n",
                min_coord.x, min_coord.y, min_coord.z, max_coord.x, max_coord.y, max_coord.z);
        printf("Min Indices: (%d, %d, %d)  ||  Max Indices: (%d, %d, %d)\n",
                start.x, start.y, start.z, end.x, end.y, end.z);
    }

    // calculate FrameOfReference
    ArrayProps props {};
    props.size = frame.size;
    props.crop_start = make_uint3(start);
    props.crop_size = make_uint3(end-start);
    if (verbose) {
        printf("crop_start: (%d, %d, %d)  ||  crop_size: (%d, %d, %d)\n",
                props.crop_start.x, props.crop_start.y, props.crop_start.z,
                props.crop_size.x, props.crop_size.y, props.crop_size.z);
    }
    return props;
}

float3 getROICentroid(Volume<uint8_t>& mask, const FrameOfReference& frame) {
    float3 centroid = make_float3(0, 0, 0);
    // find the center of mass of the identified contour volume
    float com_count = 0;
    for (unsigned int v=0; v<frame.size.x*frame.size.y*frame.size.z; v++)
    {
        if (mask[v] > 0)
        {
            int X = (v % (frame.size.x*frame.size.y)) % frame.size.x;
            int Y = (v % (frame.size.x*frame.size.y)) / frame.size.x;
            int Z = (v / (frame.size.x*frame.size.y));

            float x = frame.start.x + ((float)X + 0.5f) * frame.spacing.x;
            float y = frame.start.y + ((float)Y + 0.5f) * frame.spacing.y;
            float z = frame.start.z + ((float)Z + 0.5f) * frame.spacing.z;

            centroid += make_float3(x,y,z);
            com_count += 1.f;
        }
    }
    centroid.x /= com_count;
    centroid.y /= com_count;
    centroid.z /= com_count;

    return centroid;
}
/////////////////////
