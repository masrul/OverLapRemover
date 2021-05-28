void applyPBC(float &dx, float &dy, float &dz,float* box){
    
    if (dx > 0.5*box[0]){
        dx -=box[0]; 
    }
    else if(dx < -0.5*box[0]){
        dx +=box[0]; 
    }

    if (dy > 0.5*box[1]){
        dy -=box[1]; 
    }
    else if(dy < -0.5*box[1]){
        dy +=box[1]; 
    }

    if (dz > 0.5*box[2]){
        dz -=box[2]; 
    }
    else if(dz < -0.5*box[2]){
        dz +=box[2]; 
    }
}
