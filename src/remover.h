#ifndef REMOVER_H

#define REMOVER_H

#include "common.h"
#include "util.h"

class Remover{
  private:
    string coord_file; 
    float rcut; 
    int target;
    
  public: 
    int nAtoms;
    int nMols; 
    int nOverlaps; 
    vector<string>symbols; 
    vector<string>resnames;     
    vector<float>x;
    vector<float>y;
    vector<float>z;
    vector<int>resids; 
    float box[3]; 
    vector<MolTracker> molTracker; 
    vector<COM> comTracker; 

    Remover(string,int,float);  
    void read(); 
    void genMolTracker(); 
    void genComTracker(); 
    void checkOverlap(); 
    bool isOverlapped(int,int); 
    void write(string);

}; 

#endif /* end of include guard REMOVER_H */

