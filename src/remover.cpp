#include "remover.h"

// Constructor 
Remover::Remover(string coord_file,int target, float rcut){
    this->coord_file = coord_file; 
    this->target = target;
    this->rcut = rcut; 
}


// Read gromacs coordinate file 
void Remover::read(){ 
    ifstream file(this->coord_file); 
    string line; 
    
    getline(file,line); 
    getline(file,line); 
    
    this->nAtoms=str2int(line); 
    
    for (int i=0; i<this->nAtoms;++i){
        getline(file,line); 
        singleAtom atom = line2coord(line); 
        
        this->resnames.push_back(atom.resname); 
        this->symbols.push_back(atom.symbol); 
        this->resids.push_back(atom.resid); 
        this->x.push_back(atom.x);
        this->y.push_back(atom.y);
        this->z.push_back(atom.z); 
    }
    getline(file,line); 
    istringstream input(line); 
    input >> box[0] >> box[1] >> box[2]; 
}

// track molecule's starting index and nAtoms etc. 
void Remover::genMolTracker(){
    this->nMols = *max_element(this->resids.begin(),this->resids.end()); 
    
    int nAtoms_mol[this->nMols]; 

    for (int i=0; i<this->nMols; ++i){ 
        nAtoms_mol[i] =0;
        MolTracker m; 
        this->molTracker.push_back(m);
    }

    for (int i=0; i<this->nAtoms; ++i){
        int resid=this->resids[i];
        nAtoms_mol[resid-1] +=1; 
    }
    
    
    int nAtoms_sofar=0;   
    for (int i=0; i<this->nMols;++i){ 
        this->molTracker[i].sIDx = nAtoms_sofar; 
        this->molTracker[i].eIDx = this->molTracker[i].sIDx + nAtoms_mol[i] -1;
        this->molTracker[i].nAtoms= nAtoms_mol[i];
        this->molTracker[i].resname = this->resnames[nAtoms_sofar]; 
        this->molTracker[i].remove = false; 
        nAtoms_sofar +=nAtoms_mol[i]; 
    }

}

// track center-of-mass (COM) of molecules and max_dist from COM 
void Remover::genComTracker(){
    
    for (int i=0; i<this->nMols; ++i){ 
        COM com; 
        this->comTracker.push_back(com);
    }
    
    int startIDx; 
    int endIDx; 
    for (int i = 0; i< this->nMols; ++i){ 

        // find com 
        startIDx = this->molTracker[i].sIDx; 
        endIDx = startIDx + this->molTracker[i].nAtoms; 
        this->comTracker[i].comID = startIDx; // does not matter 

        // find max distance from comID to any ID
        vector<float> distances; 
        distances.clear();
        float xcom= this->x[this->comTracker[i].comID];
        float ycom= this->y[this->comTracker[i].comID];
        float zcom= this->z[this->comTracker[i].comID];
        for ( int j=startIDx;j<endIDx;++j){
            float dx = this->x[j] - xcom;     
            float dy = this->y[j] - ycom;     
            float dz = this->z[j] - zcom;     
            applyPBC(dx,dy,dz,box);
            float dist = dx*dx + dy*dy + dz*dz;
            distances.push_back(dist); 
        } 
        float max_dist = *max_element(distances.begin(),distances.end()); 
        max_dist = sqrt(max_dist); 
        this->comTracker[i].max_dist=max_dist; 

    }
}

// Check if iMoly and jMoly is overlapped 
bool Remover::isOverlapped(int iMoly, int jMoly){
    
    float com_dist; 
    float dist; 
    float dx;
    float dy;
    float dz; 
    float imax_dist; 
    float jmax_dist; 

    dx= this->x[this->comTracker[iMoly].comID] - this->x[this->comTracker[jMoly].comID];
    dy= this->y[this->comTracker[iMoly].comID] - this->y[this->comTracker[jMoly].comID];
    dz= this->z[this->comTracker[iMoly].comID] - this->z[this->comTracker[jMoly].comID]; 
    applyPBC(dx,dy,dz,box);
    
    com_dist = sqrt(dx*dx + dy*dy + dz*dz);   
    imax_dist = this->comTracker[iMoly].max_dist; 
    jmax_dist = this->comTracker[jMoly].max_dist; 
   
    // A quick check of overlapping using COM before jump into O(N^2) search !
    if (com_dist - imax_dist - jmax_dist > this->rcut){
        return false;
    }
    else{ 
        
        for (int i=this->molTracker[iMoly].sIDx; i <= this->molTracker[iMoly].eIDx;++i){
            for (int j=this->molTracker[jMoly].sIDx; j<= this->molTracker[jMoly].eIDx;++j){
                dx = this->x[i] -this->x[j]; 
                dy = this->y[i] -this->y[j]; 
                dz = this->z[i] -this->z[j]; 
                applyPBC(dx,dy,dz,box);
                
                dist = sqrt(dx*dx + dy*dy + dz*dz);  

                if (dist < this->rcut){ 
                    return true;
                }
            }
        }
    }
    return false; 
}


// This function loop over molecules and call isOverlapped() for all pairs
void Remover:: checkOverlap(){ 
    
    int iMoly;
        
    iMoly = this->target;  

    this->nOverlaps = 0 ; 
    for (int jMoly=0;jMoly<this->nMols;++jMoly){
        if (iMoly == jMoly)continue;
        if (this->isOverlapped(iMoly,jMoly)) { 
            this ->nOverlaps +=1; 
            this->molTracker[jMoly].remove =true; 
        }
    } 
    cout <<endl; 

    cout << "Total overlapped found: "<<this->nOverlaps<<endl;
    
} 


// Write back without overlapping with target molecule 
void Remover:: write(string outFile){ 
    
    int nAtoms_write=0; 
    for (int i=0;i<this->nMols;++i){
        if (this->molTracker[i].remove == false){
            nAtoms_write += this->molTracker[i].nAtoms; 
        }
    } 
    
    ofstream outCoord;
    outCoord.open(outFile); 

    outCoord << "Created by Masrul Huda"<<endl;
    outCoord << setw(10) <<left<<nAtoms_write<<endl; 
    
    int resID = 0; 
    int atomID = 0;
    for (int iMoly=0; iMoly<this->nMols; ++iMoly){
        if (this->molTracker[iMoly].remove)continue;
        
        int sIDx = this->molTracker[iMoly].sIDx; 
        int eIDx = this->molTracker[iMoly].eIDx; 
        ++resID; 

        for (int i=sIDx; i<=eIDx; ++i){
            ++atomID; 
            outCoord << setw(5) << right << resID; 
            outCoord << setw(5) << right << this->resnames[i]; 
            outCoord << setw(5) << right << this->symbols[i];
            outCoord << setw(5) << right << atomID;
            outCoord << fixed<< setw(8) << setprecision(3)<<right << this->x[i];
            outCoord << fixed<< setw(8) << setprecision(3)<<right << this->y[i];
            outCoord << fixed<< setw(8) << setprecision(3)<<right << this->z[i];
            outCoord << endl; 
        }
    }
    
    outCoord << fixed<< setw(10) << setprecision(5)<<right << this->box[0];
    outCoord << fixed<< setw(10) << setprecision(5)<<right << this->box[1];
    outCoord << fixed<< setw(10) << setprecision(5)<<right << this->box[2];
    outCoord << endl; 
}
