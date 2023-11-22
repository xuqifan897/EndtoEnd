#include "binary_io.h"
#include "brain_defs.h"
#include "paths.h"
#include "kernel.h"
#include "debugLog.h"

#include "cuda_runtime.h"
#include "helper_cuda.h"
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int old::load_density(SHM_DATA* data)
{
    std::string denseFile = (fs::path(Paths::Instance()->temp_dir()) / fs::path("density.raw")).string();
    std::ifstream file(denseFile);
    if (! file)
    {
        std::cerr << "Failed to open density file: " << denseFile << std::endl;
        return false;
    }
    data->density.resize(data->size_data);
    file.read((char*)(data->density.data()), data->size_data*sizeof(float));
    file.close();
    return true;
}

int old::load_data(CONSTANTS* host, SHM_DATA* data)
{
    data->size_data = host->nvoxels();

    if(load_density(data)<0)
    {
        std::cout << "Failed to load density data." << std::endl;
        return false;
    }

    KERNEL kern;
    std::string kern_fname = Paths::Instance()->cum_kernel_file();
    KERNEL::readFromFile(kern, kern_fname);
    data->radial_boundary.resize(kern.nradii);
    memcpy(data->radial_boundary.data(), kern.radial_boundary, kern.nradii*sizeof(float));
    data->kernel_array.resize(kern.ntheta*kern.nradii);
    memcpy(data->kernel_array.data(), kern.total_matrix.data(), kern.ntheta*kern.nradii*sizeof(float));
    host->ntheta = kern.ntheta;
    host->nradii = kern.nradii;

    int numangles = host->nphi * kern.ntheta/2;
    load_convolution_theta_angles(host->conv_theta_deg, numangles);
    load_convolution_phi_angles(host->conv_phi_deg, numangles);

    return 0;
}

int old::writeREVTerma(const float* d_revTerma, 
    const float* d_revDens, CONSTANTS* constants,
    const fs::path& debugDir)
{
    int volume = constants->max_rev_size.x * constants->max_rev_size.y 
        * constants->max_rev_size.z;
    std::vector<float> h_revTerma(volume);
    checkCudaErrors(cudaMemcpy((void*)(h_revTerma.data()), d_revTerma, volume*sizeof(float), cudaMemcpyDeviceToHost));

    if (! fs::is_directory(debugDir))
        fs::create_directories(debugDir);

    fs::path TermaFile = debugDir / fs::path("REVTerma.bin");
    std::ofstream f(TermaFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << TermaFile.string() << std::endl;
        return 1;
    }
    f.write((char*)(h_revTerma.data()), volume*sizeof(float));
    f.close();

    checkCudaErrors(cudaMemcpy((void*)(h_revTerma.data()), d_revDens, volume*sizeof(float), cudaMemcpyDeviceToHost));
    fs::path DensFile = debugDir / fs::path("REVDens.bin");
    f.open(DensFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << DensFile.string() << std::endl;
        return 1;
    }
    f.write((char*)(h_revTerma.data()), volume*sizeof(float));
    f.close();
    return 0;
}

int old::writeREVDose(const cudaTextureObject_t revDoseTex, 
    CONSTANTS* constants, const fs::path& debugDir)
{
    int volume = constants->max_rev_size.x * constants->max_rev_size.y 
        * constants->max_rev_size.z;
    std::vector<float> h_revDose(volume);
    float* d_revDose;
    checkCudaErrors(cudaMalloc((void**)(&d_revDose), volume*sizeof(float)));

    dim3 blockSize(8, 8, 8);
    dim3 gridSize;
    gridSize.x = int(std::ceil((float)constants->max_rev_size.x / blockSize.x));
    gridSize.y = int(std::ceil((float)constants->max_rev_size.y / blockSize.y));
    gridSize.z = int(std::ceil((float)constants->max_rev_size.z / blockSize.z));
    readTexture3D<<<gridSize, blockSize>>>(d_revDose, revDoseTex, constants->max_rev_size.x, 
        constants->max_rev_size.y, constants->max_rev_size.z);
    checkCudaErrors(cudaMemcpy(h_revDose.data(), d_revDose, volume*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(d_revDose));

    if (! fs::is_directory(debugDir))
        fs::create_directories(debugDir);
    fs::path DoseFile = debugDir / fs::path("REVDose.bin");
    std::ofstream f(DoseFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << DoseFile.string();
        return 1;
    }
    f.write((char*)(h_revDose.data()), volume*(sizeof(float)));
    f.close();
    return 0;
}

int old::writeREVDebug(float* d_debugLog, CONSTANTS* constants, const fs::path& debugDir)
{
    int volume = constants->max_rev_size.x * constants->max_rev_size.y 
        * constants->max_rev_size.z;
    std::vector<float> h_debugLog(volume);
    checkCudaErrors(cudaMemcpy(h_debugLog.data(), d_debugLog, volume*sizeof(float), cudaMemcpyDeviceToHost));
    fs::path DebugFile = debugDir / fs::path("REVDebug.bin");
    std::ofstream f(DebugFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << DebugFile.string();
        return 1;
    }
    f.write((char*)(h_debugLog.data()), volume*sizeof(float));
    f.close();
    return 0;
}

__global__ void readSurface3D(float* dest, cudaSurfaceObject_t surf,
    int width, int hight, int depth)
{
    int idx_x = threadIdx.x + blockIdx.x * blockDim.x;
    int idx_y = threadIdx.y + blockIdx.y * blockDim.y;
    int idx_z = threadIdx.z + blockIdx.z * blockDim.z;
    if (idx_x >= width || idx_y >= hight || idx_z >= depth)
        return;
    int mem_idx = idx_x + width * (idx_y + hight * idx_z);
    dest[mem_idx] = surf3Dread<float>(surf, idx_x*sizeof(float), idx_y, idx_z);
}

int old::writeREVSurf(const cudaSurfaceObject_t surfDose, 
    CONSTANTS* constants, const boost::filesystem::path& debugDir)
{
    int volume = constants->max_rev_size.x * constants->max_rev_size.y * constants->max_rev_size.z;
    std::vector<float> h_surfDose(volume);
    float* d_surfDose;
    checkCudaErrors(cudaMalloc((void**)(&d_surfDose), volume*sizeof(float)));
    dim3 blockSize(8, 8, 8);
    dim3 gridSize;
    gridSize.x = int(std::ceil(static_cast<float>(constants->max_rev_size.x)/blockSize.x));
    gridSize.y = int(std::ceil(static_cast<float>(constants->max_rev_size.y)/blockSize.y));
    gridSize.z = int(std::ceil(static_cast<float>(constants->max_rev_size.z)/blockSize.z));
    readSurface3D<<<gridSize, blockSize>>>(d_surfDose, surfDose, 
        constants->max_rev_size.x, constants->max_rev_size.y, constants->max_rev_size.z);
    checkCudaErrors(cudaMemcpy(h_surfDose.data(), d_surfDose, volume*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(d_surfDose));

    fs::path file = debugDir / fs::path("REVDoseSurf.bin");
    std::ofstream f(file.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << file.string() << std::endl;
        return 1;
    }
    f.write((char*)(h_surfDose.data()), volume*sizeof(float));
    f.close();
    return 0;
}

int old::writeBEVDose(const float* d_bevDose, const PILLAR_GRID& hPG, 
    const fs::path& debugDir)
{
    int volume = hPG.pillar_grid_nvoxels();
    std::vector<float> h_bevDose(volume);
    checkCudaErrors(cudaMemcpy((void*)(h_bevDose.data()), 
        (void*)d_bevDose, volume*sizeof(float), cudaMemcpyDeviceToHost));
    
    if (! fs::is_directory(debugDir))
        fs::create_directories(debugDir);
    fs::path DoseFile = debugDir / fs::path("BEVDose.bin");
    std::ofstream f(DoseFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << DoseFile.string();
        return 1;
    }
    f.write((char*)(h_bevDose.data()), volume*sizeof(float));
    f.close();
    return 0;
}

int old::writeResults(const RES_LOG& results)
{
    const auto& resultDir = Paths::Instance()->result_dir();
    for (int i=0; i<results.size(); i++)
    {
        fs::path BeamFolder = resultDir / fs::path(std::string("beam") + std::to_string(i));
        if (! fs::is_directory(BeamFolder))
            fs::create_directories(BeamFolder);
        const auto& beam = results[i];
        const auto& hPG = std::get<0>(beam);
        const auto& beamlets = std::get<1>(beam);
        for (int j=0; j<beamlets.size(); j++)
        {
            fs::path BeamletFile = BeamFolder / fs::path(std::string("beamlet") + std::to_string(j));
            std::ofstream f(BeamletFile.string());
            if (! f)
            {
                std::cerr << "Could not open file: " << BeamletFile.string() << std::endl;
                return 1;
            }
            f.write((char*)(beamlets[j].data()), beamlets[j].size()*sizeof(float));
            f.close();
        }
        // log metadata
        fs::path MetaDataFile = BeamFolder / fs::path("metadata.txt");
        std::ofstream f(MetaDataFile.string());
        if (! f)
        {
            std::cerr << "Could not open file: " << MetaDataFile.string() << std::endl;
            return 1;
        }
        f << "Number of beamlets: " << hPG.numBeamlets << std::endl;
        f << "Pillar size: " << hPG.pillarSize << std::endl;
        f << "Pillar Dimension: " << hPG.pillarDims << std::endl;
        f.close();
    }
    return 0;
}