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
    mxArray *V_tmp;
    double *W;
    
    mxArray *PneigMatlab;
    mwSize *PneigDims;
    int PneigLength=0;
    int PneigCount=0;
    double t=0;
    double w,wsum;
    double nx,ny,nz;
    double x,y,z;
    double a,b;
    double *ap, *bp;
    
    const double *itertmp=NULL;
    int iterN=1;
    
    /* Outputs */

    double *V_smooth;
    double *V_smooth_tmp;
    
    double *Pneig;

    /* Table Size */
    mwSize table_Dims[2]={1,3};
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;


    /* Loop variables */
    int i, j, k, n;
    
    /* Check for proper number of arguments. */
    if(nrhs != 6) {
        mexErrMsgTxt("6 inputs are required.");
        /*mexErrMsgTxt("2 inputs are required.");*/
    } else if(nlhs!=1) {
        mexErrMsgTxt("1 output is required");
    }
    
    
    /* Connect Inputs */
    V=(double *)mxGetPr(prhs[0]);
    Ne=prhs[1];
    W=(double *)mxGetPr(prhs[2]);
    
    /* iterN argument given? */

    itertmp = (const double*) mxGetData(prhs[3]);
    iterN = *itertmp;
    if (!mxIsInf(*itertmp) &&
        !mxIsNaN(*itertmp))
        iterN = (int)*itertmp;
    
    ap = (const double*)mxGetData(prhs[4]);
    bp = (const double*)mxGetData(prhs[5]);
    a = ap[0];
    b = bp[0];
    
    /* Get number of VertexN */
    VertexDims = mxGetDimensions(prhs[0]);
   /* mexPrintf("%d %d\n",VertexDims[0],VertexDims[1]);*/

    VertexN=VertexDims[0];

    /* Reserve memory */
    table_Dims[0]=VertexN; table_Dims[1]=3;
    plhs[0]= mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    V_tmp = mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);
    /* Connect Outputs */
    
    V_smooth = (double *)mxGetPr(plhs[0]);
    V_smooth_tmp = (double *)mxGetPr(V_tmp);

    memcpy(V_smooth,V, VertexDims[0]*VertexDims[1]*sizeof(double));
    
    for (n=0; n<iterN; n++) {
        memcpy(V_smooth_tmp,V_smooth, VertexDims[0]*VertexDims[1]*sizeof(double));
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
            
            x=0;
            y=0;
            z=0;
            nx=0;
            ny=0;
            nz=0;
            wsum=0;
            PneigCount=0;    
            for (k=0; k<PneigLength; k++) {
                t = V_smooth_tmp[(int)(Pneig[k]-1)];
                if(mxIsFinite(t)) {
                    w = W[(int)(Pneig[k]-1)];
                    nx += w*t;
                    ny += w*V_smooth_tmp[(int)((Pneig[k]-1) + VertexN*1)];
                    nz += w*V_smooth_tmp[(int)((Pneig[k]-1) + VertexN*2)];
                    wsum += w;
                    PneigCount ++;
                }
            }
            
            /*
            if(PneigCount > 0)
            {
                nx/=PneigCount;
                ny/=PneigCount;
                nz/=PneigCount;
            }*/
            if(wsum > 0)
            {
                nx/=wsum;
                ny/=wsum;
                nz/=wsum;
            }
                    
            t = V_smooth_tmp[i];
            if(mxIsFinite(t))
            {
                w = W[i];
                x = w*t;
                y = w*V_smooth_tmp[(int)(i + VertexN*1)];
                z = w*V_smooth_tmp[(int)(i + VertexN*2)];                
                /*PneigCount ++;*/
            }
            if(PneigCount > 0)
            {
                V_smooth[i] = (1-a)*x + a*nx;
                V_smooth[(int)(i + VertexN*1)] = (1-a)*y + a*ny;
                V_smooth[(int)(i + VertexN*2)] = (1-a)*z + a*nz;
            } else {
                V_smooth[i] = V_smooth_tmp[i];
                V_smooth[(int)(i + VertexN*1)] = V_smooth_tmp[(int)(i + VertexN*1)];
                V_smooth[(int)(i + VertexN*2)] = V_smooth_tmp[(int)(i + VertexN*2)];

                /* If the center vertex and ALL of its neighbors are NaN, put NaN there
                 * .... oh wait... that would be the same as putting V[i] back */

                /* V_diffuse[i] = mxGetNaN(); */
            }
        }
        
    }
    mxDestroyArray(V_tmp);
    /* Remove temporary memory */

}

