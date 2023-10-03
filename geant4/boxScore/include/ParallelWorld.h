#ifndef ParallelDetectorConstruction_h
#define ParallelDetectorConstruction_h 1

#include "G4VUserParallelWorld.hh"
#include "G4LogicalVolume.hh"

namespace bs
{
    class ParallelWorld : public G4VUserParallelWorld
    {
    public:
        ParallelWorld(G4String& parallelWorldName, float offset, int dimz);
        virtual ~ParallelWorld() = default;
        float& getOffset() {return this->Offset;}
        int& getDimZ() {return this->DimZ;}

    protected:
        virtual void Construct();
        virtual void ConstructSD();
    
    private:
        float Offset;  // the offset of the starting voxel
        int DimZ;  // the dimension along the z direction
        int DimXY;  // the dimension along the x and y axes
        float voxelSize;  // in half dimension
        std::vector<G4LogicalVolume*> logicals;
    };
}

#endif