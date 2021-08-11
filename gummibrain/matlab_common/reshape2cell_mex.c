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
    double *List;

    /* Neighbour list input */
    const mxArray *Ne;
    mxArray *PneigMatlab;
    mwSize *PneigDims;
    int PneigLength=0;

    double *Pneig;
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;

    const mwSize *ListDims;
    int ListN=0;

    /* Loop variables */
    int i, j, k, n;
    
    /* Check for proper number of arguments. */
    if(nrhs!=2) {
        mexErrMsgTxt("2 inputs are required.");
    } else if(nlhs!=0) {
        mexErrMsgTxt("0 output is required");
    }
    
    
    /* Connect Inputs */
    List=(double *)mxGetPr(prhs[0]);
    Ne=prhs[1];
    
    
    ListDims = mxGetDimensions(prhs[0]);
    ListN=ListDims[0]*ListDims[1];
    
    /* Get number of VertexN */
    VertexDims = mxGetDimensions(prhs[1]);
    VertexN=VertexDims[0]*VertexDims[1];

    j=0;
    for (i=0; i<VertexN; i++) {
        PneigMatlab=mxGetCell(Ne, i);
        if( PneigMatlab == NULL)
        {
            PneigLength=-1;
        }
        else
        {
            PneigDims=(mwSize *)mxGetDimensions(PneigMatlab);
            PneigLength=(PneigDims[0]*PneigDims[1]);
            Pneig=(double *)mxGetPr(PneigMatlab);
        }

        if(j > ListN-PneigLength) {
            mexErrMsgTxt("Not enough elements in list");
            break;
        }
        
        for (k=0; k<PneigLength; k++) {
            Pneig[k]=List[j++];
        }
    }
    
    /* Remove temporary memory */

}

