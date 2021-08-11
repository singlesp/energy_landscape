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

#include "rt_nonfinite.h"
#include "ordfilt3_med.h"
#include "fastmedc.h"

/* Fast, destructive median selection, credit to Nicolas Devillard for
   the original C implementation. */
#define SWAP(a,b) {register real_T t=(a);(a)=(b);(b)=t; }
/*#define SWAP(x, y) do { typeof(x) temp##x##y = x; x = y; y = temp##x##y; } while (0)*/
/*#define SWAP(a,x,y) {register real_T t=a->data[x]; a->data[x]=a->data[y]; a->data[y]=t; }*/

static real_T quickSelect(real_T *arr, int32_T low, int32_T high, boolean_T exact);

static real_T quickSelect(real_T *arr, int32_T low, int32_T high, boolean_T exact)
{

    int32_T high0 = high;
    real_T mednext;
    real_T tmp;
    boolean_T oddlength = ((high-low+1) % 2) == 1;

    int32_T median;
    int32_T middle, ll, hh;
    median = (low + high) / 2;

    
    for (;;)
    {
        /* One element only */
        if (high <= low)
            break;

        /* Two elements only */
        if (high == low + 1)
        {
            if (arr[low] > arr[high])
		SWAP(arr[low],arr[high]);
            break;
        }

        /* Find median of low, middle and high items; swap to low position */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])
            SWAP(arr[middle],arr[high]);
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



void fastmedc(const emxArray_real_T *x, emxArray_real_T *mOutput)
{
    int32_T samples, dimensions, elements, i;
    boolean_T exact = false;

    samples    = x->size[0];
    dimensions = x->size[1];
    elements   = samples * dimensions;

    exact = true;

    /* Create output signal, get pointer access */
    /*M_OUT = mxCreateDoubleMatrix(1,dimensions, mxREAL);*/
    /*mOutput = mxGetPr(M_OUT);*/


    for(i = 0; i < dimensions; i++){
        mOutput->data[i] = quickSelect(x->data, i*samples,(i+1)*samples-1,exact);
    }


    /* Release allocated memory */

    return;
}
