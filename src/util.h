#ifndef UTIL_H

#define UTIL_H

#include "common.h"

struct singleAtom{
    int atomid;  
    string resname; 
    string symbol; 
    int resid; 
    float x; 
    float y; 
    float z;
}; 

struct MolTracker{
    string resname; 
    int nAtoms ;
    int sIDx;  // starting index of molecule  
    int eIDx;  // end index of molecule 
    bool remove; // logical varible, if need to be removed 
}; 

// Center-of-Mass 
struct COM{
    int comID; 
    float max_dist; 
}; 

extern int str2int(string str);
extern float str2float(string str); 
extern singleAtom line2coord(string line); 
extern void applyPBC(float &dx, float &dy, float &dz,float* box); 
    


#endif /* end of include guard COMMON_H */

