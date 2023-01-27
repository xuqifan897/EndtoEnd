#ifndef __DUMMY_DEFS_H__
#define __DUMMY_DEFS_H__

// CUDA utilities and system includes
#include <helper_math.h>
#include "DoseCalcIO/beam.h"

#define MAXIMUM_SERVERS_PER_NODE    4

#define DEFAULT_DEVICE_COUNT    6
#define DEFAULT_SERVER_COUNT    2
#define DEFAULT_STREAM_COUNT    1
#define DEFAULT_CONVO_COUNT     32

#define DEFAULT_NPHI            8
#define DEFAULT_NTHETA          8
#define DEFAULT_DELR            0.2
#define DEFAULT_VOXEL           0.2, 0.2, 0.2
#define DEFAULT_SIZE            128, 128, 128
#define DEFAULT_BEAM_ORIG       0., 0., 0.
#define DEFAULT_BEAM_SIZE       5., 0., 5.
#define DEFAULT_CTBASE          "slab"
#define DEFAULT_SPECFILE        "spec_6mv"


// this class was created for passing pertinent information between machines over sockets
// I'm pretty sure it isn't used anywhere in the mgcs-precomp or mgcs-omnidose programs
class MGCS_PARAMS
{
  public:
    unsigned int nphi;      // azimuthal angles
    unsigned int ntheta;    // zenithal angles
    float delr;             // radial step size
    uint3 size;             // dimensions of data in voxels
    float3 voxel;           // size of voxel in mm
    BEAM beam;
    char CTbase[32];
    char SPECfile[32];

    int BID;               // ID of the branching client node // PARENT
    int SID;               // ID of the server node // CHILD

    int sconvos[MAXIMUM_SERVERS_PER_NODE];           // number of convolutions for each server on this node

    int sservers;          // number of servers to launch on this node
    int sbranches;         // number of branches off node BID, if nbranches = 0, server is a leaf

    int nnodes;            // maximum nodes in mgcs
    int nconvos;           // total number of convolutions

    int nservers;          // maximum servers per node
    int ndevices;          // GPUs per server

    int nstreams;          // number of streams per GPU

    bool verbose;
    bool timing;

    void calculate_servers()
    {
        int node_convos = nconvos / nnodes;
        if (node_convos < 1) node_convos = 1;
        int server_convos = node_convos / nservers;
        if (server_convos < 1) server_convos = 1;
        printf("\n %d convolutions / %d nodes = %d convolutions per node + %d extra", nconvos, nnodes,
                                                                                      node_convos, nconvos % nnodes );
        printf("\n %d convolutions / %d servers = %d convolutions per server + %d extra", node_convos, nservers,
                                                                                          server_convos, node_convos % nservers );

        int convos_so_far = (SID-1) * node_convos;
        int convos_left = nconvos - convos_so_far;
        printf("\n Convolutions so far = %d",convos_so_far);
        printf("\n Convolutions left = %d",convos_left);

        if (convos_left > 0)
        {
            int s = 0;
            while(convos_left >= server_convos && s < nservers)
            {
                sconvos[s] = server_convos;
                convos_left -= server_convos;
                s++;
            }
            sservers = s;
            if (convos_left > 0 && s < nservers)
            {
                sconvos[s] = convos_left;
                convos_left = 0;
                sservers++;
            }
        }
        else
            sservers = 0;
        printf("\n %d servers for node %d",sservers,SID);
        for (int s=0; s<sservers; s++)
            printf("\n     %d convos for server %d",sconvos[s],s);
    }

    void calculate_branches()
    {
        int nextID = 2*SID + 1;

        int node_convos = nconvos / nnodes;
        if (node_convos < 1) node_convos = 1;

        int convos_so_far = (nextID-1) * node_convos;
        int convos_left = nconvos - convos_so_far;
        printf("\n %d ID: Convolutions so far = %d", nextID, convos_so_far);
        if (convos_left > 0)
            sbranches = 2;
        else
        {
            nextID = 2*SID;
            convos_so_far = (nextID-1) * node_convos;
            convos_left = nconvos - convos_so_far;
            printf("\n %d ID: Convolutions so far = %d", nextID, convos_so_far);
            if (convos_left > 0)
                sbranches = 1;
        }
        printf("\n %d branches for node %d",sbranches,SID);
    }

    void create_options_string( char * options )
    {
        sprintf(options,"%d %d\n%f\n%d %d %d\n%f %f %f\n%f %f %f\n%s\n%s\n%d %d\n%d %d\n%d %d %d %d %d\n%d %d\n%d %d %d %d",
                        nphi, ntheta,
                        delr,
                        size.x, size.y, size.z,
                        voxel.x, voxel.y, voxel.z,
                        beam.source.x, beam.source.y, beam.source.z,
                        CTbase,
                        SPECfile,
                        BID, SID,
                        sservers, sbranches,
                        nconvos, nservers, ndevices, nstreams, nnodes,
                        verbose, timing,
                        sconvos[0], sconvos[1], sconvos[2], sconvos[3] );
    }

    void parse_options_string( char * options )
    {
        // parse options
        sscanf(options,"%d %d\n%f\n%d %d %d\n%f %f %f\n%f %f %f\n%s\n%s\n%d %d\n%d %d\n%d %d %d %d %d\n%d %d\n%d %d %d %d",
                       &nphi, &ntheta,
                       &delr,
                       &size.x, &size.y, &size.z,
                       &voxel.x, &voxel.y, &voxel.z,
                       &beam.source.x, &beam.source.y, &beam.source.z,
                       &CTbase,
                       &SPECfile,
                       &BID, &SID,
                       &sservers, &sbranches,
                       &nconvos, &nservers, &ndevices, &nstreams, &nnodes,
                       &verbose, &timing,
                       &sconvos[0], &sconvos[1], &sconvos[2], &sconvos[3] );
    }
};

//macro to access grid values
#define ARRAY_VALUE(data, i, j, k, count)\
    (data[i + count.x * (j + k * count.y)])

#endif // __DUMMY_DEFS_H__






