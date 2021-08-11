#if !defined(_WIN32)
#define dgels dgels_
#endif

#include "mex.h"
#include "lapack.h"

#ifndef PI
#  define PI 3.141592653589793238
#endif


/* Coordinates to index */
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizx*sizy+y*sizx+x;}
int mindex2(int x, int y, int sizx) { return y*sizx+x;}

double wrapTo2Pi(double lambda) { return lambda-(2*PI)*floor(lambda/(2*PI)); }
double wrapToPi(double lambda) { return wrapTo2Pi(lambda+PI)-PI; };

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /* Vertex list input */
    double *Apolar, *Aeccen, *Vsph;
    /* Neighbour list input */
    const mxArray *Ne;
    mxArray *PneigMatlab;
    double *Pneig;
    mwSize *PneigDims;
    int PneigLength=0;

    /* Matrix division variables*/

		size_t m,n,p,lwork;
    char *trans="N";
    double *A2, *B2, *W;
    mxArray *Awork, *Bwork, *work;    
    mwSignedIndex info;

    /* Outputs */

    double *Vfieldsign;

    /* Table Size */
    mwSize table_Dims[2]={1,3};
    
    /* Number of vertices */
    const mwSize *VertexDims;
    int VertexN=0;

    /* Loop variables */
    int i, k;
    
    /* Check for proper number of arguments. */
    if(nrhs!=4) {
        /*mexErrMsgTxt("At least 4 inputs are required.");*/
        mexErrMsgTxt("4 inputs are required.");
    } else if(nlhs!=1) {
        mexErrMsgTxt("1 output is required");
    }
    
    /* */
    
    /* Connect Inputs */
    Aeccen=(double *)mxGetPr(prhs[0]);
    Apolar=(double *)mxGetPr(prhs[1]);
    Vsph=(double *)mxGetPr(prhs[2]);
    Ne=prhs[3];

    /* Get number of VertexN */
 
    VertexDims = mxGetDimensions(prhs[0]);
    VertexN=VertexDims[0];

    /* Reserve memory */
    table_Dims[0]=VertexN; table_Dims[1]=1;
    plhs[0] = mxCreateDoubleMatrix(table_Dims[0], table_Dims[1], mxREAL);
    /*plhs[0]= mxCreateNumericArray(2, table_Dims, mxDOUBLE_CLASS, mxREAL);*/
    /* Connect Outputs */
    
    Vfieldsign = (double *)mxGetPr(plhs[0]);

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

				m = PneigLength;
				p = 3;
				n = 2;
				lwork=m*p*2;
				
				work = mxCreateDoubleMatrix(lwork, 1, mxREAL);
				W = mxGetPr(work);
				
				Awork = mxCreateDoubleMatrix(m, p, mxREAL);
				A2 = mxGetPr(Awork);

				Bwork = mxCreateDoubleMatrix(m, n, mxREAL);   
				B2 = mxGetPr(Bwork);

        for (k=0; k<PneigLength; k++) {
            A2[k]=Vsph[i]-Vsph[(int)Pneig[k]-1]; /*du*/
            A2[k+m]=Vsph[i+VertexN]-Vsph[(int)Pneig[k]-1+VertexN]; /*dv*/
            A2[k+2*m]=1;
            
            B2[k]=wrapToPi(Aeccen[i]-Aeccen[(int)Pneig[k]-1]); /*dR*/
            B2[k+m]=wrapToPi(Apolar[i]-Apolar[(int)Pneig[k]-1]); /*dTh*/
        }

        dgels(trans,&m,&p,&n,A2,&m,B2,&m,W,&lwork,&info);
        Vfieldsign[i] = B2[0]*B2[m+1]-B2[m+0]*B2[1];

				mxDestroyArray(Awork);
				mxDestroyArray(Bwork);
				mxDestroyArray(work);        
    }

    /* Remove temporary memory */
    
}

