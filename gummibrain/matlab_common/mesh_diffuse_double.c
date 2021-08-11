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
    
    const double *itertmp=NULL;
    int iterN=1;
    
    /* Outputs */

    double *V_diffuse;
    
    double *Pneig;

    /* Table Size */
    mwSize table_Dims[2]={1,3};
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;


    /* Loop variables */
    int i, j, k, n;
    
    /* Check for proper number of arguments. */
    if(!(nrhs==2 || nrhs == 3)) {
        mexErrMsgTxt("At least 2 inputs are required.");
        /*mexErrMsgTxt("2 inputs are required.");*/
    } else if(nlhs!=1) {
        mexErrMsgTxt("1 output is required");
    }
    
    
    /* Connect Inputs */
    V=(double *)mxGetPr(prhs[0]);
    Ne=prhs[1];
    
    
    /* iterN argument given? */
    if ((nrhs > 2) &&
        (mxGetClassID(prhs[2]) == mxDOUBLE_CLASS) &&
        (mxGetNumberOfElements(prhs[2]) == 1)) {
            
        /* check for value := 4 */
        itertmp = (const double*) mxGetData(prhs[2]);
        iterN = *itertmp;
        if (!mxIsInf(*itertmp) &&
            !mxIsNaN(*itertmp))
            
            /* overwrite default tps value */
            iterN = (int)*itertmp;
    }
    
    /* Get number of VertexN */
    VertexDims = mxGetDimensions(prhs[0]);
    VertexN=VertexDims[0];

    /* Reserve memory */
    table_Dims[0]=VertexN; table_Dims[1]=VertexDims[1];
    plhs[0]= mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    /* Connect Outputs */
    
    V_diffuse = (double *)mxGetPr(plhs[0]);

    for (n=0; n<iterN; n++) {

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
            
            for(j=0;j<VertexDims[1];j++)
            {
                tmpV=0;
                PneigCount=0;    
                for (k=0; k<PneigLength; k++) {
                    t = V[j*VertexN + (int)Pneig[k]-1];
                    if(mxIsFinite(t)) {
                        tmpV += t;
                        PneigCount ++;
                    }
                }
                t = V[j*VertexN + i];
                if(mxIsFinite(t))
                {
                    tmpV += t;
                    PneigCount ++;
                }   
                if(PneigCount > 0)
                {
                    V_diffuse[j*VertexN + i] = tmpV/PneigCount;
                } else {
                    V_diffuse[j*VertexN + i] = V[j*VertexN + i];

                    /* If the center vertex and ALL of its neighbors are NaN, put NaN there
                     * .... oh wait... that would be the same as putting V[i] back */

                    /* V_diffuse[i] = mxGetNaN(); */
                }
            }
            /* V_diffuse[i] = (V_diffuse[i] + V[i])/2; */
        }
        if(iterN > 1 && n < iterN-1)
        {
            memcpy(V, V_diffuse, VertexDims[0]*VertexDims[1]*mxGetElementSize(prhs[0]));
        }
    }

    /* Remove temporary memory */

}

