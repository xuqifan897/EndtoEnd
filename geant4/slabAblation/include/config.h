#ifndef PHANTOM
#define PHANTOM 4

// WORLD means that the sensitive detector is only associated to the WORLD logical volume
    // It turns out that the sensitive detector does not produce hits
// SLABS means that each slab logical volume is associated to a separate sensitive detector
    // It turns out that only the sensitive detector of the last layer works
// SHARED means that all slab logical volumes are associated to the same sensitive detector
#define WORLD 0
#define SENSDET 2

// whether to deposit energy in the pre-step or the post step
// 0 means pre-step, 1 means post-step, 2 means middle point
#define STEP_ORDER 0

// to sparsify the material density and enlarge voxel size
#define SPARSE 100

#define RES 

#endif
