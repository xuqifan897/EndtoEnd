#ifndef __PROFILE_H__
#define __PROFILE_H__

#ifdef USE_NVTX

#include "nvToolsExt.h"
#include "nvToolsExtCuda.h"
#include "nvToolsExtCudaRt.h"

const uint32_t colors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff, 0xff00ffff, 0xffff0000, 0xffffffff };
const int num_colors = sizeof(colors)/sizeof(uint32_t);

#define NVTX_PUSH_RANGE(name,cid) { \
    int color_id = cid; \
    color_id = color_id%num_colors;\
    nvtxEventAttributes_t eventAttrib = {0}; \
    eventAttrib.version = NVTX_VERSION; \
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
    eventAttrib.colorType = NVTX_COLOR_ARGB; \
    eventAttrib.color = colors[color_id]; \
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
    eventAttrib.message.ascii = name; \
    nvtxRangePushEx(&eventAttrib); \
}
#define NVTX_POP_RANGE nvtxRangePop()
#define NVTX_NAME_STREAM(stream, name) nvtxNameCudaStreamA(stream, name)

#else

#define NVTX_PUSH_RANGE(name,cid)
#define NVTX_POP_RANGE
#define NVTX_NAME_STREAM(stream, name)

#endif

#endif // __PROFILE_H__
