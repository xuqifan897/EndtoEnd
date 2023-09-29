#include "manage_gpu.cuh"

#include <vector>
#include <helper_cuda.h>

// set CUDA device according to command line prompt or default to device 0
void setActiveGPU(unsigned int dev) {
    checkCudaErrors( cudaSetDevice(dev) );
    checkCudaErrors( cudaDeviceReset() );  // consider removing this - reset will kill any currently running GPU process
}

inline bool IsGPUCapableP2P(cudaDeviceProp *pProp)
{
#ifdef _WIN32
    return (bool)(pProp->tccDriver ? true : false);
#else
    return (bool)(pProp->major >= 2);
#endif
}
inline bool IsAppBuiltAs64()
{
    return sizeof(void*) == 8;
}

void cudaMemInfo::formatted_memory_print(int gpuid, const char* memo, float q1, float q2, const char* unit) {
    printf("||| (Device %d) - Memory %s: %0.3f %s of %0.3f %s |||\n", gpuid, memo, q1, unit, q2, unit);
}

cudaMemInfo query_device_memory() {
    cudaMemInfo meminfo;
    size_t free_byte, total_byte;
    checkCudaErrors( cudaMemGetInfo(&free_byte, &total_byte) );

    meminfo.free  = (float)free_byte/(1024.0*1024.0*1024.0);
    meminfo.total = (float)total_byte/(1024.0*1024.0*1024.0);
    meminfo.used  = meminfo.total - meminfo.free;
    meminfo.unit  = "GB";

    return meminfo;
}
cudaMemInfo query_device_memory(int gpuid) {
    checkCudaErrors( cudaSetDevice(gpuid) );
    cudaMemInfo meminfo = query_device_memory();
    meminfo.gpuid = gpuid;
    return meminfo;
}

void init_devices(int &nDevices, int *gpuid_arr, const int maximum_device_count, const int ndev_requested, int first_device, bool verbose) {
    // Skip UVA initiailization
    checkCudaErrors( cudaGetDeviceCount(&nDevices) );

    // check for valid first_device request
    if (first_device < 0 || first_device >= nDevices) {
        first_device = 0;
    }

    for (int dev=0; dev<nDevices; dev++) {
        gpuid_arr[dev] = dev + first_device;
    }
}
void init_devices_uva(int &nDevices, int *gpuid_arr, const int maximum_device_count, const int ndev_requested, int first_device, bool verbose) {
    /* if gpucount is 1, dont perform UVA init
     * otherwise, perform UVA checks and init on all devices including/after gpuid=first_device and return
     * gpuid_arr containing the gpuid for each initialized device
     */

    // check 64bit build
    if (!IsAppBuiltAs64()) {
        printf("UVA only supported on 64-bit applications");
        exit(1);
    }

    // get number of devices...
    checkCudaErrors(cudaGetDeviceCount(&nDevices));
    // ...then reduce if necessary

    // check for valid first_device request
    if (first_device < 0 || first_device >= nDevices) {
        first_device = 0;
    }

    if (verbose) {
        printf("Number of devices [requested: %d | identified: %d]\n", ndev_requested, nDevices);
        if (first_device>0) {
            printf("  (requesting device IDs >= %d)\n", first_device);
        }
    }

    // check each for uva, enable on those capable and return array of capable device ids
    std::vector<cudaDeviceProp> prop(maximum_device_count);
    std::vector<int> supports_id(maximum_device_count);
    int gpucount = 0;
    for (int i=first_device; i<nDevices && gpucount<maximum_device_count && (ndev_requested < 0 || gpucount<ndev_requested); i++) {
        checkCudaErrors(cudaGetDeviceProperties(&prop[i], i));
        // Check that supports uva
        if (IsGPUCapableP2P(&prop[i])) {
            supports_id[gpucount++] = i;
        }
    }

    if (gpucount > 1) {
        // enable uva on all supported devices
        int uva_count = 0;
        for (int i=0; i<gpucount; i++) {
            bool use_gpu = true;
            for (int j=0; j<gpucount; j++) {
                if (i == j) { continue; }

                // check that device supports_id[i] has peer access to all others
                int can_access_peer_ab;
                checkCudaErrors(cudaDeviceCanAccessPeer(&can_access_peer_ab, supports_id[i], supports_id[j]));
                int can_access_peer_ba;
                checkCudaErrors(cudaDeviceCanAccessPeer(&can_access_peer_ba, supports_id[j], supports_id[i]));
                if (!can_access_peer_ab || !can_access_peer_ba) {
                    use_gpu = false;
                    break;
                }
            }
            if (use_gpu) {
                gpuid_arr[uva_count++] = supports_id[i];
            }
        }
        printf("Number of UVA capable devices: %d of %d\n", uva_count, nDevices);
        for (int i=0; i<uva_count; i++) {
            for (int j=i+1; j<uva_count; j++) {
                // enable uva for this pair of devices
                checkCudaErrors(cudaSetDevice(gpuid_arr[i]));
                checkCudaErrors(cudaDeviceEnablePeerAccess(gpuid_arr[j], 0));
                checkCudaErrors(cudaSetDevice(gpuid_arr[j]));
                checkCudaErrors(cudaDeviceEnablePeerAccess(gpuid_arr[i], 0));
                if (verbose) {
                    //printf("> %s (GPU%d) supports UVA: %s\n", prop[gpuid_arr[i]].name, gpuid_arr[i], (prop[gpuid_arr[i]].unifiedAddressing ? "Yes" : "No"));
                    //printf("> %s (GPU%d) supports UVA: %s\n", prop[gpuid_arr[j]].name, gpuid_arr[j], (prop[gpuid_arr[j]].unifiedAddressing ? "Yes" : "No"));
                    printf("  -->P2P enabled between GPUs: %d <--> %d\n", gpuid_arr[i], gpuid_arr[j]);
                }
            }
        }
        // return true number of devices in use
        nDevices = uva_count;
    } else {
        // return true number of devices in use
        nDevices = 1;
        gpuid_arr[0] = supports_id[0];
    }
}
