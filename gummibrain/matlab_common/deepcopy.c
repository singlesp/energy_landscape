#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
     const mxArray *prhs[])
{
    plhs[0] = mxDuplicateArray(prhs[0]);
}
