#ifndef __RTIMAGE_H__
#define __RTIMAGE_H__

#include <string>
#include <array>

#include <helper_cuda.h>
#include <helper_math.h>
#include "dcmtk/dcmdata/dctk.h"

#define RTIMAGE_SOP_CLASS_UID "1.2.840.10008.5.1.4.1.1.481.1"
#define CTIMAGE_SOP_CLASS_UID "1.2.840.10008.5.1.4.1.1.2"
#define MRIMAGE_SOP_CLASS_UID "1.2.840.10008.5.1.4.1.1.4"

#ifndef UINT16_MAX
#define UINT16_MAX 65535
#endif

class RTImage
{
public:
    RTImage() {}
    RTImage(const std::string& dname, bool verbose=false) : dicom_dir{dname} {
        loadDicomInfo(verbose);
        loadRTImageData(verbose);
    }
    ~RTImage();

    class SLICE_DATA {
      public:
        std::string filename;
        std::string sop_instance_uid;
        std::string reference_frame_uid;
        int instance_number;
        float slice_location;
        std::string patient_position;
        float3 image_position_patient;
        std::array<float, 6> image_orientation_patient;

        SLICE_DATA() {}

        // copy constructor
        SLICE_DATA(const SLICE_DATA& copyfrom) {
            filename = copyfrom.filename;
            sop_instance_uid = copyfrom.sop_instance_uid;
            reference_frame_uid = copyfrom.reference_frame_uid;
            instance_number = copyfrom.instance_number;
            slice_location = copyfrom.slice_location;
            image_position_patient = copyfrom.image_position_patient;
            patient_position = std::string(copyfrom.patient_position);
            image_orientation_patient = copyfrom.image_orientation_patient;
        }
    };

    bool loadDicomInfo(bool verbose=false);
    int  loadRTImageData(bool verbose=false, bool flipXY=false);
    void saveRTImageData (const std::string& outpath, bool anonymize_switch );
    int  saveRTImageData (const std::string& outpath, float *newData, bool anonymize_switch );
    bool importSOPClassUID(const std::string& fname);
    void importPatientPosition(unsigned int i);
    void importImagePositionPatient( unsigned int i );
    void importInstanceNumber( unsigned int i );
    void importPatientInfo();
    void anonymize( DcmDataset *dataset );

    /* GETTERS */
    std::string  getDicomDirectory()                          { return dicom_dir; };
    float        getArrayVoxel(int i, int j, int k)           { return data_array[i + data_size.x*(j + data_size.y*k)]; };
    unsigned int getImageCount()                              { return image_count; };
    int3         getDataSize()                                { return data_size; };
    float3       getVoxelSize()                               { return voxel_size; };
    float        getRescaleSlope()                            { return rescale_slope; };
    float        getRescaleIntercept()                        { return rescale_intercept; };
    float        getDataMin()                                 { return data_min; };
    float        getDataMax()                                 { return data_max; };
    float*       getDataArray()                               { return data_array; };
    float        getWindowCenter()                            { return window_center; };
    float        getWindowWidth()                             { return window_width; };
    int          getSliceInstanceNumber(unsigned int i)       { return slice[i].instance_number; };
    std::string  getSlicePatientPosition(unsigned int i)      { return slice[i].patient_position; }
    float3       getSliceImagePositionPatient(unsigned int i) { return slice[i].image_position_patient; };
    std::array<float, 6> getSliceImageOrientationPatient(unsigned int i) { return slice[i].image_orientation_patient; }
    std::string  getSliceSOPInstanceUID(unsigned int i)       { return slice[i].sop_instance_uid; };
    std::string  getSliceReferenceFrameUID(unsigned int i)    { return slice[i].reference_frame_uid; };
    std::string  getSliceFilename(unsigned int i)             { return slice[i].filename; };

    /* SETTERS */
    void setDicomDirectory(const std::string& dname)                                   { dicom_dir = dname; }
    void setSliceSOPInstanceUID( unsigned int i, const std::string& uid)               { slice[i].sop_instance_uid = uid; }
    void setSliceFilename(unsigned int i, const std::string& fname)                    { slice[i].filename = fname; }
    void setArrayVoxel(int x, int y, int z, float v)                                   { data_array[x + data_size.x*(y + data_size.y*z)] = v; };
    void setDataSize(int x, int y, int z)                                              { data_size = int3{x,y,z}; }
    void setVoxelSize(float v_x, float v_y, float v_z)                                 { voxel_size = float3{v_x, v_y, v_z}; }
    void setDataOrigin(float v_x, float v_y, float v_z)                                { data_origin = float3{v_x, v_y, v_z}; }
    void setRescaleSlope(float v)                                                      { rescale_slope = v; };
    void setRescaleIntercept(float v)                                                  { rescale_intercept = v; };
    void setDataMin(float v)                                                           { data_min = v; };
    void setDataMax(float v)                                                           { data_max = v; };
    void setWindowCenter(float v)                                                      { window_center = v; };
    void setWindowWidth(float v)                                                       { window_width = v; };
    void setSliceInstanceNumber(unsigned int i, int n)                                 { slice[i].instance_number = n; };
    void setSliceImagePositionPatient(unsigned int i, float v_x, float v_y, float v_z) { slice[i].image_position_patient = float3{v_x, v_y, v_z}; }

    // Type of image
    enum class SOPTYPE { RTIMAGE, CT, MR };
    SOPTYPE sop_type;
    std::string getSOPTypeName() { return getSOPTypeName(sop_type); }
    std::string getSOPTypeName(SOPTYPE t) {
        switch(t) {
            case SOPTYPE::RTIMAGE  :  return "RTIMAGE";
            case SOPTYPE::CT       :  return "CT";
            case SOPTYPE::MR       :  return "MR";
            default                :  return "unknown";
        }
    }

    // data ordering - see https://public.kitware.com/IGSTKWIKI/index.php/DICOM_data_orientation
    enum class DATAORDER { LPS, RAS };
    DATAORDER data_order;

protected:
    std::string dicom_dir;
    std::string dicom_date;
    std::string pt_series_description;
    std::string pt_name;
    std::string pt_id;
    std::string pt_study_id;
    std::string pt_study_instance_uid;
    std::string pt_series_instance_uid;

    unsigned int image_count;
    SLICE_DATA *slice;

    int3 data_size;
    float3 voxel_size;
    float3 data_origin;
    float slice_thickness;

    std::string orient{};
    float3 orient_x;
    float3 orient_y;

    float window_center;
    float window_width;
    float rescale_slope;
    float rescale_intercept;

    float data_min;
    float data_max;
    float *data_array;
};

#endif // __RTIMAGE_H__

