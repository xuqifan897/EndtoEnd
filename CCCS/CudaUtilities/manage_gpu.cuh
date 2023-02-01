/* GPU management function definitions */
#ifndef __MANAGE_GPU_H__
#define __MANAGE_GPU_H__

#include <helper_cuda.h>

// general management
void setActiveGPU(unsigned int dev=0);
/* void freeGPUMemory(unsigned int dev=0); */

// requirements checks
inline bool IsGPUCapableP2P(cudaDeviceProp *pProp);
inline bool IsAppBuiltAs64();

class cudaMemInfo {
    public:
        float total;
        float free;
        float used;
        std::string unit = "GB";
        int gpuid;

        void print_available() {
            formatted_memory_print(gpuid, "available", free, total, "GB");
        }
        void print_consumed() {
            formatted_memory_print(gpuid, "consumed", free, total, "GB");
        }

    private:
        static void formatted_memory_print(int gpuid, const char* memo, float q1, float q2, const char* unit);
};
cudaMemInfo query_device_memory();
cudaMemInfo query_device_memory(int gpuid);

/* Check if uva possible for P2P memory access */
/* return number of gpu's that are uva enabled, and array of device ids */
void init_devices(int &nDevices, int *gpuid, const int maximum_device_count, const int nrequested=-1, const int first_device=0, bool verbose=false);
void init_devices_uva(int &nDevices, int *gpuid, const int maximum_device_count, const int nrequested=-1, const int first_device=0, bool verbose=false);

#endif // __MANAGE_GPU_H__
