#include "debugLog.h"
#include "cudaInit.h"
#include <iostream>
#include <iomanip>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

int old::datavols_log(CONSTANTS* constants, SHM_DATA* datavols)
{
    // kernel array
    std::stringstream output;
    int hight = constants->ntheta;
    int width = constants->nradii;
    for (int i=0; i<hight; i++)
    {
        for (int j=0; j<width; j++)
        {
            int idx = i * width + j;
            output << std::left << std::setw(12) << 
                std::setprecision(4) << datavols->kernel_array[idx];
        }
        output << std::endl;
    }
    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("kernel_array.txt");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f << output.str();
    f.close();

    outputFile = fs::path(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("density.bin");
    f.open(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    int volume = constants->size.x * constants->size.y * constants->size.z * sizeof(float);
    f.write((char*)(datavols->density.data()), volume);
    f.close();

    output.clear();
    for (int i=0; i<datavols->radial_boundary.size(); i++)
        output << std::left << std::setw(12) << std::setprecision(4) << datavols->radial_boundary[i];
    outputFile = fs::path(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("radial_boundary.txt");
    f.open(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f << output.str();
    f.close();
    return 0;
}

int old::mono_kernel_log(MONO_KERNELS* mono)
{
    std::stringstream output;
    output << "mono->spectrum_file = " << mono->spectrum_file << std::endl;
    output << "mono->nkernels = " << mono->nkernels << std::endl;

    int width = 10;
    output << std::left << std::setw(width) << "energy" << 
        std::left << std::setw(width) << "fluence" <<
        std::left << std::setw(width) << "mu" << 
        std::left << std::setw(width) << "mu_en" << std::endl;
    for (int i=0; i<mono->nkernels; i++)
    {
        output << std::fixed << std::setprecision(2) << std::setw(width) << mono->energy[i] << 
            std::fixed << std::setprecision(4) << std::setw(width) << mono->fluence[i] << 
            std::fixed << std::setprecision(4) << std::setw(width) << mono->mu[i] << 
            std::fixed << std::setprecision(4) << std::setw(width) << mono->mu_en[i] << std::endl;
    }

    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("spectrum_log.txt");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f << output.str();
    f.close();
    return 0;
}

int old::constants_log(CONSTANTS* constants)
{
    std::stringstream output;
    output << "nphi: " << constants->nphi << std::endl;
    output << "ntheta: " << constants->ntheta << std::endl;
    output << "nradii: " << constants->nradii;

    output << std::endl << "conv_theta_deg: (";
    for (int i=0; i<constants->conv_theta_deg.size(); i++)
    {
        output << constants->conv_theta_deg[i];
        if (i < constants->conv_theta_deg.size()-1)
            output << ", ";
    }
    output << ")" << std::endl;

    output << std::endl << "conv_phi_deg: (";
    for (int i=0; i<constants->conv_phi_deg.size(); i++)
    {
        output << constants->conv_phi_deg[i];
        if (i < constants->conv_phi_deg.size()-1)
            output << ", ";
    }
    output << ")" << std::endl << std::endl;

    output << "ss_factor: " << constants->ss_factor << std::endl;
    output << "rev_latspacing: " << constants->rev_latspacing << std::endl;
    output << "rev_longspacing: " << constants->rev_longspacing << std::endl;
    output << "kernel_extent: " << constants->kernel_extent << std::endl;
    output << "penumbra: " << constants->penumbra << std::endl;
    output << "beamhard_correct: " << constants->beamhard_correct << std::endl;
    output << "beam_count: " << constants->beam_count << std::endl;
    output << "start: (" << constants->start.x << ", " 
        << constants->start.y << ", " << constants->start.z << ")" << std::endl;
    output << "voxel: (" << constants->voxel.x << ", " 
        << constants->voxel.y << ", " << constants->voxel.z << ")" << std::endl;
    output << "size: (" << constants->size.x << ", " 
        << constants->size.y << ", " << constants->size.z << ")" << std::endl;
    output << "max_rev_size: (" << constants->max_rev_size.x << ", " 
        << constants->max_rev_size.y << ", " << constants->max_rev_size.z << ")" << std::endl;
    output << "calc_bbox_start: (" << constants->calc_bbox_start.x << ", " 
        << constants->calc_bbox_start.y << ", " << constants->calc_bbox_start.z << ")" << std::endl;
    output << "calc_bbox_size: (" << constants->calc_bbox_size.x << ", " 
        << constants->calc_bbox_size.y << ", " << constants->calc_bbox_size.z << ")" << std::endl;
    
    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("constants_log.txt");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f << output.str();
    f.close();
    return 0;
}

int old::beams_log(std::vector<BEAM>& beams)
{
    std::stringstream output;
    output << "Number of beams: " << beams.size() << std::endl;
    output << "Take the first beam for example." << std::endl;
    BEAM& firstBeam = beams[0];
    output << "uid: " << firstBeam.uid << std::endl;
    output << "azimuth: " << firstBeam.azimuth << std::endl;
    output << "zenith: " << firstBeam.zenith << std::endl;
    output << "coll: " << firstBeam.coll << std::endl;
    output << "source: (" << firstBeam.source.x << ", " << 
        firstBeam.source.y << ", " << firstBeam.source.z << ")" << std::endl;
    output << "direction: (" << firstBeam.direction.x << ", " << 
        firstBeam.direction.y << ", " << firstBeam.direction.z << ")" << std::endl;
    output << "isocenter: (" << firstBeam.isocenter.x << ", " << 
        firstBeam.isocenter.y << ", " << firstBeam.isocenter.z << ")" << std::endl;
    output << "sad: " << firstBeam.sad << std::endl;
    output << "fmap_size: (" << firstBeam.fmap_size.x << ", " << 
        firstBeam.fmap_size.y << ")" << std::endl;
    output << "beamlet_size: (" << firstBeam.beamlet_size.x << ", " << 
        firstBeam.beamlet_size.y << ")" << std::endl;
    
    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("beams_log.txt");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f << output.str();
    f.close();
    return 0;
}

__global__ void old::readTexture2D(float* output, cudaTextureObject_t texObj, int width, int hight)
{
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
    if (x>=width || y>=hight)
        return;
    uint idx = (x + width * y);
    output[idx] = tex2D<float>(texObj, (float)x+0.5, (float)y+0.5);
}

__global__ void old::readTexture3D(float* output, cudaTextureObject_t texObj, int width, int hight, int depth)
{
    uint x = threadIdx.x + blockIdx.x * blockDim.x;
    uint y = threadIdx.y + blockIdx.y * blockDim.y;
    uint z = threadIdx.z + blockIdx.z * blockDim.z;
    if (x>=width || y>=hight || z>=depth)
        return;
    uint idx = x + (y + z * hight) * width;
    output[idx] = tex3D<float>(texObj, x+0.5, y+0.5, z+0.5);
}

int old::texKern_log(CONSTANTS* constants)
{
    int width = constants->nradii;
    int hight = constants->ntheta;
    float* d_result;
    checkCudaErrors(cudaMalloc((void**)(&d_result), width*hight*sizeof(float)));
    dim3 gridSize(1, 1, 1);
    dim3 blockSize(width, hight, 1);
    readTexture2D<<<gridSize, blockSize>>>(d_result, texKern, width, hight);

    // output result
    float* h_result = new float[width*hight];
    checkCudaErrors(cudaMemcpy(h_result, d_result, width*hight*sizeof(float), cudaMemcpyDeviceToHost));
    // convert binary result to text result
    std::stringstream result;
    for (int i=0; i<hight; i++)
    {
        for (int j=0; j<width; j++)
            result << std::left << std::setw(12) << std::setprecision(4) << h_result[i * width + j];
        result << std::endl;
    }
    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("texKern_log.txt");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile << std::endl;
        return 1;
    }
    f << result.str();
    f.close();

    checkCudaErrors(cudaFree((void*)d_result));
    delete[] h_result;
    return 0;
}

int old::texSpectrum_log(MONO_KERNELS* mono_kernels)
{
    int width = mono_kernels->nkernels;
    int hight = 4;
    float* d_result;
    checkCudaErrors(cudaMalloc((void**)(&d_result), width*hight*sizeof(float)));
    dim3 gridSize(1, 1, 1);
    dim3 blockSize(width, hight, 1);
    readTexture2D<<<gridSize, blockSize>>>(d_result, texSpectrum, width, hight);

    // output result
    float* h_result = new float[width*hight];
    checkCudaErrors(cudaMemcpy(h_result, d_result, width*hight*sizeof(float), cudaMemcpyDeviceToHost));
    std::stringstream result;
    for (int i=0; i<hight; i++)
    {
        for (int j=0; j<width; j++)
            result << std::left << std::setw(12) << std::setprecision(4) << h_result[i * width + j];
        result << std::endl;
    }
    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("texSpectrum_log.txt");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f << result.str();
    f.close();

    checkCudaErrors(cudaFree((void*)d_result));
    delete[] h_result;
    return 0;
}

int old::texDens_log(CONSTANTS* constants)
{
    auto& size = constants->size;
    uint volume=size.x*size.y*size.z*sizeof(float);
    float* d_result;
    checkCudaErrors(cudaMalloc((void**)(&d_result), volume));
    dim3 blockSize(8, 8, 8);
    dim3 gridSize;
    gridSize.x = int(std::ceil((float)size.x / blockSize.x));
    gridSize.y = int(std::ceil((float)size.y / blockSize.y));
    gridSize.z = int(std::ceil((float)size.z / blockSize.z));
    readTexture3D<<<gridSize, blockSize>>>(d_result, texDens, size.x, size.y, size.z);

    // output result
    float* h_result = new float[size.x*size.y*size.z];
    checkCudaErrors(cudaMemcpy(h_result, d_result, volume, cudaMemcpyDeviceToHost));
    fs::path outputFile(Paths::Instance()->result_dir());
    outputFile = outputFile / fs::path("texDens_log.bin");
    std::ofstream f(outputFile.string());
    if (! f)
    {
        std::cerr << "Could not open file: " << outputFile.string() << std::endl;
        return 1;
    }
    f.write((char*)h_result, volume);
    f.close();
    checkCudaErrors(cudaFree((void*)d_result));
    delete[] h_result;
    return 0;
}

int old::debugLog(SHM_DATA* datavols, MONO_KERNELS *mono, CONSTANTS *constants, std::vector<BEAM>& beams)
{
    if (datavols_log(constants, datavols))
        return 1;
    if (mono_kernel_log(mono))
        return 1;
    if (constants_log(constants))
        return 1;
    if (beams_log(beams))
        return 1;
    if (texKern_log(constants))
        return 1;
    if (texSpectrum_log(mono))
        return 1;
    if (texDens_log(constants))
        return 1;
    return 0;
}