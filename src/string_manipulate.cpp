#include "util.h" 


int str2int(string str){
    int intVal; 
    stringstream intStr(str); 
    intStr >> intVal;
    return intVal; 
} 

float str2float(string str){
    float floatVal; 
    stringstream floatStr(str); 
    floatStr >> floatVal; 
    return floatVal; 
}

// It takes a line from gromacs file, and return symbol,x,y,z etc. 
singleAtom line2coord(string line){ 
    // "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" 
    singleAtom atom; 
    string sub_str;  

    // atom id 
    sub_str = line.substr(0,5); 
    atom.resid=str2int(sub_str);  

    // resname  
    atom.resname = line.substr(5,5); 
    
    // symbol 
    atom.symbol = line.substr(10,5); 

    // resid 
    sub_str = line.substr(15,5); 
    atom.atomid=str2int(sub_str);  
    
    // x 
    sub_str = line.substr(20,8); 
    atom.x=str2float(sub_str); 

    // y 
    sub_str = line.substr(28,8); 
    atom.y=str2float(sub_str);  

    // z 
    sub_str = line.substr(36,8); 
    atom.z=str2float(sub_str);  
    
    return atom; 

}
