#include "remover.h"


int main(int argc, char* argv[])
{  
    
    string inputFile;
    string outFile;
    float rcut=-1.0;
    int target; 
    
    // Parese Command Line Argument 
    int i=1;
    while(i<argc && argv[i][0] == '-') {
        string opt = string(argv[i]) ;
        if(opt == "-i") inputFile=argv[++i] ;
        if(opt == "-o") outFile=argv[++i] ;
        if(opt == "-rcut") rcut = atof(argv[++i]) ;
        if(opt == "-target") target = atoi(argv[++i]) ;
        ++i ;
    }
    
    // Inistantiate Remover object and process 
    if (rcut < 0.0)rcut =0.35;  
    Remover remover(inputFile,target,rcut);     

    remover.read(); 
    remover.genMolTracker();
    remover.genComTracker(); 
    remover.checkOverlap();
    remover.write(outFile); 
    
    return 0;
}
