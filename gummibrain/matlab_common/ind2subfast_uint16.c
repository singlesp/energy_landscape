#include "mex.h"
#include "math.h"

/*
 * function [ET_table,EV_table,ETV_index]=edge_tangents(V,Ne)
 */

/* Coordinates to index */
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}
int mindex2(int x, int y, int sizx) { return y*sizx+x;}

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* Vertex list input */
    unsigned short *gridsize;

    /* Outputs */
    unsigned short *sub;
        
    /* Table Size */
    mwSize table_Dims[2]={1,3};
    const mwSize *GridDims;
        
    /* Loop variables */
    unsigned int gridN;
    unsigned int ijk;
    unsigned short i,j,k;    

    
    /* Check for proper number of arguments. */
    if(nrhs!=1) {
        /*mexErrMsgTxt("At least 2 inputs are required.");*/
        mexErrMsgTxt("1 input is required.");
    } else if(nlhs!=1) {
        mexErrMsgTxt("1 output is required");
    }
    
    /* Connect Inputs */
    gridsize=(unsigned short *)mxGetPr(prhs[0]);
    GridDims = mxGetDimensions(prhs[0]);
    
    gridN=1;
    for(i=0;i<(unsigned short)(GridDims[0]*GridDims[1]);i++)
        gridN=gridN*gridsize[i];
    
    /* Reserve memory */
    table_Dims[0]=gridN; table_Dims[1]=3;
    plhs[0]= mxCreateNumericArray(2, table_Dims, mxUINT16_CLASS, mxREAL);
    /* Connect Outputs */

    sub = (unsigned short *)mxGetPr(plhs[0]);
 
    ijk=0;
    for (k=0; k<gridsize[2]; k++) {
        for (j=0; j<gridsize[1]; j++) {
            for (i=0; i<gridsize[0]; i++) {
                sub[ijk] = i+1;
                sub[ijk+gridN] = j+1;
                sub[ijk+gridN*2] = k+1;
                ijk++;
            }
        }
    }

}

