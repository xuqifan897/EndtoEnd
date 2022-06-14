#include "args.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

using namespace E2E;
using namespace std;

CCCSkernel::CCCSkernel(int NA): num_angles(NA)
{
    this->angles = (float*)malloc(NA*sizeof(float));
    this->Atheta = (float*)malloc(NA*sizeof(float));
    this->Btheta = (float*)malloc(NA*sizeof(float));
    this->atheta = (float*)malloc(NA*sizeof(float));
    this->btheta = (float*)malloc(NA*sizeof(float));

    for(int i=0; i<this->num_angles; i++)
    {
        this->angles[i] = 0;
        this->Atheta[i] = 0;
        this->Btheta[i] = 0;
        this->atheta[i] = 0;
        this->btheta[i] = 0;
    }
}

CCCSkernel::~CCCSkernel()
{
    if (this->angles != nullptr)
        free(this->angles);
    if (this->Atheta != nullptr)
        free(this->Atheta);
    if (this->Btheta != nullptr)
        free(this->Btheta);
    if (this->atheta != nullptr)
        free(this->atheta);
    if (this->btheta != nullptr)
        free(this->btheta);
}

CCCSkernel::CCCSkernel(CCCSkernel& old): num_angles{old.num_angles}
{
    if (old.angles != nullptr)
    {
        this->angles = (float*)malloc(this->num_angles*sizeof(float));
        for (int i=0; i<this->num_angles; i++)
            this->angles[i] = old.angles[i];
    }
    else
        this->angles = nullptr;
    
    if (old.Atheta != nullptr)
    {
        this->Atheta = (float*)malloc(this->num_angles*sizeof(float));
        for (int i=0; i<this->num_angles; i++)
            this->Atheta[i] = old.Atheta[i];
    }
    else
        this->Atheta = nullptr;
    
    if (old.Btheta != nullptr)
    {
        this->Btheta = (float*)malloc(this->num_angles*sizeof(float));
        for (int i=0; i<this->num_angles; i++)
            this->Btheta[i] = old.Btheta[i];
    }
    else
        this->Btheta = nullptr;

    if (old.atheta != nullptr)
    {
        this->atheta = (float*)malloc(this->num_angles*sizeof(float));
        for (int i=0; i<this->num_angles; i++)
            this->atheta[i] = old.atheta[i];
    }
    else
        this->atheta = nullptr;
    
    if (old.btheta != nullptr)
    {
        this->btheta = (float*)malloc(this->num_angles*sizeof(float));
        for (int i=0; i<this->num_angles; i++)
            this->btheta[i] = old.btheta[i];
    }
    else
        this->btheta = nullptr;
}

CCCSkernel::CCCSkernel(CCCSkernel&& old): num_angles(old.num_angles)
{
    this->angles = exchange(old.angles, nullptr);
    this->Atheta = exchange(old.Atheta, nullptr);
    this->Btheta = exchange(old.Btheta, nullptr);
    this->atheta = exchange(old.atheta, nullptr);
    this->btheta = exchange(old.btheta, nullptr);
}

FCBBkernel::FCBBkernel(int ND): num_depths(ND)
{
    this->depths = (float*)malloc(this->num_depths*sizeof(float));
    this->doses = (float*)malloc(this->num_depths*sizeof(float));
    for (int i=0; i<this->num_depths; i++)
    {
        this->depths[i] = 0;
        this->doses[i] = 0;
    }

    this->A = 0;
    this->B = 0;
    this->a = 0;
    this->b = 0;

    this->min_depth = 0;
    this->max_depth = 0;
    this->d_doses = 0;
}

FCBBkernel::~FCBBkernel()
{
    if (this->depths != nullptr)
        free(this->depths);
    if (this->doses != nullptr)
        free(this->doses);
    texDecon();
}

FCBBkernel::FCBBkernel(FCBBkernel& old): \
num_depths(old.num_depths), A(old.A), B(old.B), a(old.a), b(old.b)
{
    if (old.depths != nullptr)
    {
        this->depths = (float*)malloc(this->num_depths*sizeof(float));
        for (int i=0; i<this->num_depths; i++)
            this->depths[i] = old.depths[i];
    }
    else
        this->depths = nullptr;

    if (old.doses != nullptr)
    {
        this->doses = (float*)malloc(this->num_depths*sizeof(float));
        for (int i=0; i<this->num_depths; i++)
            this->doses[i] = old.doses[i];
    }
    else
        this->doses = nullptr;
}

FCBBkernel::FCBBkernel(FCBBkernel&& old): \
num_depths(old.num_depths), A(old.A), B(old.B), a(old.a), b(old.b)
{
    this->depths = exchange(old.depths, nullptr);
    this->doses = exchange(old.doses, nullptr);
}

void FCBBkernel::texInit()
{
    if (this->d_doses != 0)
        checkCudaErrors(cudaFreeArray(this->d_doses));
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( \
        32, 0, 0, 0, cudaChannelFormatKindFloat);
    checkCudaErrors(cudaMallocArray( \
        &(this->d_doses), &channelDesc, this->num_depths));
    checkCudaErrors(cudaMemcpyToArray(this->d_doses, 0, 0, \
        this->doses, this->num_depths*sizeof(float), cudaMemcpyHostToDevice));

    cudaResourceDesc texRes;
    memset(&texRes, 0, sizeof(cudaResourceDesc));
    texRes.resType = cudaResourceTypeArray;
    texRes.res.array.array = this->d_doses;

    cudaTextureDesc texDescr;
    memset(&texDescr, 0, sizeof(cudaTextureDesc));
    texDescr.normalizedCoords = true;
    texDescr.filterMode = cudaFilterModeLinear;
    texDescr.addressMode[0] = cudaAddressModeBorder;
    texDescr.readMode = cudaReadModeElementType;

    checkCudaErrors(cudaCreateTextureObject(&(this->tex), &texRes, &texDescr, NULL));
}

void E2E::FCBBkernel::texDecon()
{
    if (this->tex)
        checkCudaErrors(cudaDestroyTextureObject(this->tex));
    if (this->d_doses != 0)
        checkCudaErrors(cudaFreeArray(this->d_doses));
}

extern "C" void depthDose(dim3 gridSize, dim3 blockSize, float* output, \
    uint nPoints, cudaTextureObject_t& texObj);

void E2E::testDepthDose(FCBBkernel* kernel)
{
    (*kernel).texInit();
    int nPoints = 512;
    dim3 gridSize(8);
    dim3 blockSize(64);
    float* h_output = (float*)malloc(nPoints*sizeof(float));
    float* d_output = nullptr;
    checkCudaErrors(cudaMalloc((void**)&d_output, nPoints*sizeof(float)));
    depthDose(gridSize, blockSize, d_output, nPoints, (*kernel).tex);

    checkCudaErrors(cudaMemcpy(h_output, d_output, nPoints*sizeof(float), \
        cudaMemcpyDeviceToHost));
    string outputFile{"/data/qifan/projects_qlyu/EndtoEnd3/data/patient1_out/depthDose.dat"};
    ofstream outFile(outputFile);
    if (! outFile.is_open())
    {
        cout << "Could not open this file: " << outputFile << endl;
        exit;
    }
    outFile.write((char*)h_output, nPoints*sizeof(float));
    outFile.close();
}