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
    double *V;
    /* Neighbour list input */
    const mxArray *Ne;
    mxArray *PneigMatlab;
    mwSize *PneigDims;
    int PneigLength=0;
    int PneigCount=0;
    double t=0;
    double tmpV=0;
    
    const double *startidx_tmp=NULL;
    int startidx=0;
    
    const double *itertmp=NULL;
    int iterN=1000;
    
    int verbose=1;
    
    mxArray *V_idx1;
    double *Pidx1;
    mxArray *V_idx2;
    double *Pidx2;    
    int idx1N=0;
    int idx2N=0;
    
    int patchN=0;
    
    /* Outputs */

    double *V_fill;
    
    double *Pneig;

    /* Table Size */
    mwSize table_Dims[2]={1,3};
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;


    /* Loop variables */
    int i, j, k, n;
    
    /* Check for proper number of arguments. */
    if(nrhs<3) {
        /*mexErrMsgTxt("At least 2 inputs are required.");*/
        mexErrMsgTxt("2 inputs are required.");
    } else if(nlhs!=1) {
        mexErrMsgTxt("1 output is required");
    }
    
    
    /* Connect Inputs */
    V=(double *)mxGetPr(prhs[0]);
    Ne=prhs[1];


    /* check for value := 4 */
    startidx_tmp = (const double*) mxGetData(prhs[2]);
    if (!mxIsInf(*startidx_tmp) &&
        !mxIsNaN(*startidx_tmp))
        /* overwrite default tps value */
        startidx = (int)*startidx_tmp;
    else
        mexErrMsgTxt("Invalid startidx.");

    
    /* iterN argument given? */
    if ((nrhs > 3) &&
        (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[3]) == 1)) {
            
        /* check for value := 4 */
        itertmp = (const double*) mxGetData(prhs[3]);
        iterN = (int)*itertmp;
        if (!mxIsInf(*itertmp) &&
            !mxIsNaN(*itertmp))
            
            /* overwrite default tps value */
            iterN = (int)*itertmp;
    }
    
        /* iterN argument given? */
    if (nrhs > 4) {
        verbose=0;
    }
    
    /* Get number of VertexN */
    VertexDims = mxGetDimensions(prhs[0]);
    VertexN=VertexDims[0];

    /* Reserve memory */
    table_Dims[0]=VertexN; table_Dims[1]=1;
    plhs[0]= mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    /* Connect Outputs */
    
    V_idx1 = mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    Pidx1 = (double *)mxGetPr(V_idx1);
    V_idx2 = mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    Pidx2 = (double *)mxGetPr(V_idx2);    
    
    /* make index 0-based */
    startidx--; 
    
    Pidx1[0] = startidx;
    idx1N = 1;
    patchN = 1;
    
    V_fill = (double *)mxGetPr(plhs[0]);
    memcpy(V_fill, V, VertexDims[0]*VertexDims[1]*sizeof(double));
    
    V_fill[startidx] = 2;
    
    for (n=0; n<iterN; n++) {
        idx2N = 0;
        for (j=0; j<idx1N; j++) {
            
            i = (int)Pidx1[j];
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
            
            for (k=0; k<PneigLength; k++) {
                if(V_fill[(int)Pneig[k]-1] == 1)
                {
                    V_fill[(int)Pneig[k]-1] = 2;
                    Pidx2[idx2N++] = (int)Pneig[k]-1;
                }

            }

        }
        if(idx2N == 0){
            break;
        }
        patchN += idx2N;
        memcpy(Pidx1,Pidx2, idx2N*sizeof(double));
        idx1N = idx2N;
    }

    /* Remove temporary memory */
    mxDestroyArray(V_idx1);
    mxDestroyArray(V_idx2);
    
    if(verbose > 0){
        if(n == iterN){
            mexPrintf("Maximum iterations (%d) reached.  Quitting.\n",iterN);
        } else {
           mexPrintf("Filled %d vertices in %d iterations.\n",patchN,n);
        }    
    }
}

