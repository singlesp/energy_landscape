/*
    Core MEX routine implementing a fast version of the classical 1D running
    median filter of size W, where W is odd.

    Usage:
    m = fastmedfilt1d_core(x, xic, xfc, W2)

    Input arguments:
    - x          Input signal
    - xic        Initial boundary condition
    - xfc        End boundary condition
    - W2         Half window size (so that window size is W=2*W2+1)

    Output arguments:
    - m          Median filtered signal

    (c) Max Little, 2010. If you use this code, please cite:
    Little, M.A. and Jones, N.S. (2010),
    "Sparse Bayesian Step-Filtering for High-Throughput Analysis of Molecular
    Machine Dynamics"
    in Proceedings of ICASSP 2010, IEEE Publishers: Dallas, USA.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "matrix.h"

/* Real variable type */
#define  REAL        double

/* Input parameters */
#define  X_IN        prhs[0]
#define  B_EXACT     prhs[1]
#define CELL_IN      prhs[2]

/* Output parameters */
#define  M_OUT       plhs[0]

/* Function definition */
#define  SYNTAX   "m = fastmed1d_core(X,exact)"

/* Fast, destructive median selection, credit to Nicolas Devillard for
   the original C implementation. */
#define SWAP(a,b) { register REAL t=(a);(a)=(b);(b)=t; }

REAL quickSelect(REAL arr[], long low, long high, bool exact)
{

    long high0 = high;
    REAL mednext;
        
    long median;
    long middle, ll, hh;
    median = (low + high) / 2;

    bool oddlength = ((high-low+1) % 2) == 1;
    
    for (;;)
    {
        /* One element only */
        if (high <= low)
            break;

        /* Two elements only */
        if (high == low + 1)
        {
            if (arr[low] > arr[high])
                SWAP(arr[low], arr[high]);
            break;
        }

        /* Find median of low, middle and high items; swap to low position */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            SWAP(arr[middle], arr[high]);
        if (arr[low] > arr[high])
            SWAP(arr[low], arr[high]);
        if (arr[middle] > arr[low])
            SWAP(arr[middle], arr[low]);

        /* Swap low item (now in position middle) into position (low+1) */
        SWAP(arr[middle], arr[low+1]);

        /* Work from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;)
        {
            do
                ll++;
            while (arr[low] > arr[ll]);
            do
                hh--;
            while (arr[hh] > arr[low]);

            if (hh < ll)
                break;

            SWAP(arr[ll], arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        SWAP(arr[low], arr[hh]);

        /* Reset active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
    if(oddlength || !exact)
        return arr[median];
    else
    {
        mednext=arr[median+1];
        for(ll=median+2;ll<=high0;ll++)
        {
            if(arr[ll] < mednext)
                mednext = arr[ll];
        }
        return (arr[median]+mednext)/2;
    }

}




/* Main entry point */
/* lhs - output parameters */
/* rhs - input parameters */
void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
#if !defined(V4_COMPAT)
    const mxArray *prhs[]         /* array of pointers to input arguments */
#else
    mxArray *prhs[]         /* array of pointers to input arguments */
#endif
)
{
    long samples, dimensions, elements, i;
    bool exact = false;
    REAL            *scalar;            /* Reused pointer to get access to scalar input parameters */
    REAL            *xInput;            /* Input signal vector */
    REAL            *mOutput;           /* Median filtered output signal vector */
    mxClassID category;
    const mxArray *cell_element_ptr;
    int total_num_of_cells;
    
    /* Check for proper number of arguments */
    if (!((nrhs == 1 || nrhs == 2 || nrhs == 3) && (nlhs == 1)))
    {
        mexErrMsgTxt("Incorrect number of parameters.\nSyntax: "SYNTAX);
    }

    samples    = mxGetM(X_IN);
    dimensions = mxGetN(X_IN);
    elements   = samples * dimensions;

    xInput = mxGetPr(X_IN);
    if(nrhs >= 2)
        exact = mxGetPr(B_EXACT) > 0;
    
    /*if(nrhs >= 3)
    {
        category = mxGetClassID(CELL_IN);
        if(category == mxCELL_CLASS){
            total_num_of_cells = mxGetNumberOfElements(CELL_IN); 
            for(i=0;i<total_num_of_cells;i++){
                cell_element_ptr = mxGetCell(CELL_IN, i);
                mexPrintf("%d elements\n",mxGetNumberOfElements(cell_element_ptr));
            }
        }
    } else {
        */
        /* Create output signal, get pointer access */
        M_OUT = mxCreateDoubleMatrix(1,dimensions, mxREAL);
        mOutput = mxGetPr(M_OUT);
        
        
        for(i = 0; i < dimensions; i++){
            mOutput[i] = quickSelect(xInput, i*samples,(i+1)*samples-1,exact);
        }
    /*}*/


    /* Release allocated memory */

    return;
}
