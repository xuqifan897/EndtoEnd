#ifndef __PRECOMP_DEBUG_H__
#define __PRECOMP_DEBUG_H__

/************* SWITCHABLES **************/
/****************************************/
//            DEBUG SWITCHES            //
/****************************************/
// #define DEBUG_PRINT_FINDCNTRVOLUME

/****************************************/
//              IMPLEMENT               //
/****************************************/
#ifdef DEBUG_PRINT_FINDCNTRVOLUME
    #define DEBUG(cmd) cmd;
#else
    #define DEBUG(cmd) {}
#endif

/****************************************/
/****************************************/

#endif // __PRECOMP_DEBUG_H__
