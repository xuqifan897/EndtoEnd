#ifndef PHANTOM
#define PHANTOM 3

// WORLD means that the sensitive detector is only associated to the WORLD logical volume
    // It turns out that the sensitive detector does not produce hits
// SLABS means that each slab logical volume is associated to a separate sensitive detector
    // It turns out that only the sensitive detector of the last layer works
// SHARED means that all slab logical volumes are associated to the same sensitive detector
#define WORLD 0
#define SENSDET SHARED

#define RES 

#endif
