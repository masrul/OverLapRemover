#ifndef UTIL_H

#define UTIL_H

#include "common.h"

struct singleAtom{
    int atomid;  
    std::string resname; 
    std::string symbol; 
    int resid; 
    float x; 
    float y; 
    float z;
}; 

struct MolTracker{
    std::string resname; 
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

extern int str2int(std::string str);
extern float str2float(std::string str); 
extern singleAtom line2coord(std::string line); 
extern void applyPBC(float &dx, float &dy, float &dz,float* box); 
    


#endif /* end of include guard COMMON_H */

