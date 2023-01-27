#ifndef __ROI_H__
#define __ROI_H__

#include <memory>

#include "../RTClasses/rtstruct.h"
#include "Utilities/logging.h"
#include "./volume.h"
#include "./macros.h"
#include "./io_data_structs.h"

#include "H5Cpp.h" // Namespace H5::

///// DATA STRUCTURES /////
struct StructureSet {
    unsigned int       sub_cntr_count = 0; // nslices in "CLOSED_PLANAR" DICOM storage mode
    std::vector<uint>  sub_cntr_points_count;  // list of coord. counts for each sub_cntr
    unsigned int       total_points_count = 0;
    std::vector<float> points;
};

class BaseROIMask {
    public:
        BaseROIMask() {}
        BaseROIMask(std::string name, ArrayProps props) : name{name}, props(props) {}
        virtual ~BaseROIMask() = 0;

        std::string name;       // name of ROI structure
        ArrayProps props;       // custom bbox containing ROI for more efficient iterating

        virtual uint8_t operator[](uint64_t idx) = 0;
        uint nvoxels() { return props.nvoxels(); }
        virtual uint nones() = 0;
        uint nzeros() { return nvoxels() - nones(); }

        int writeToFile(std::string fname, bool verbose=false);
        virtual int _writeToHDF5(H5::Group& h5group) = 0;
};

class DenseROIMask : public BaseROIMask {
    public:
        DenseROIMask() {}
        DenseROIMask(std::string name, std::vector<uint8_t> mask, ArrayProps props) : BaseROIMask(name, props), mask{mask} {}
        virtual ~DenseROIMask() {}
        std::vector<uint8_t> mask; // linearized binary array (1: in ROI, 0: out ROI)
        virtual uint8_t operator[](uint64_t idx);
        virtual uint nones();

        int _writeToHDF5(H5::Group& h5group);
        static int _readFromHDF5(DenseROIMask& mask, H5::Group& h5group);

    protected:
        uint _nones_cached = 0;
};

class ROIMaskList {
    public:
        ROIMaskList() {}
        ~ROIMaskList() {}
        std::vector<std::shared_ptr<DenseROIMask> > _coll;
        void push_back(DenseROIMask* mask) { _coll.emplace_back(mask); }
        std::vector<std::string> getROINames();
        std::vector<uint64_t>    getROICapacities();

        uint64_t size() { return _coll.size(); }

        int writeToFile(std::string fname, bool verbose=false);
        static int readFromFile(ROIMaskList& masklist, std::string fname, bool verbose=false);
        int _writeToHDF5(H5::Group& h5group);
        static int _readFromHDF5(ROIMaskList& masklist, H5::Group& h5group, bool verbose=false);
};
///////////////////////////


// Search all ROIs in rtstruct by name and return the index of the matched entry
//   Supports exact matching or substring matching
int getROIIndex(RTStruct& rtstruct, const std::string& search, bool exact=false, bool verbose=false);

// Load an roi at entry "roi_idx" from "rtstruct" and package coordinates into StructureSet object
int loadStructureSet(StructureSet& roi, RTStruct& rtstruct, uint roi_idx, bool verbose=false);

// Fit a 3d rectangular box to be as small as possible while confining the "roi" within
ArrayProps getROIExtents(const StructureSet& roi, const FrameOfReference& frame, bool verbose=false);

// fill "mask" with linearized dense binary matrix (0/1) according to contours in StructureSet "roi"
// int getDenseROIMask(std::vector<uint8_t>& mask, StructureSet& roi, const FrameOfReference& frame, const ArrayProps& bbox);

float3 getROICentroid(Volume<uint8_t>& mask, const FrameOfReference& frame);


#endif // __ROI_H__
