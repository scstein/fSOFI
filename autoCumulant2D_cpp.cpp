// Author: Simon Christoph Stein
// E-Mail: scstein@phys.uni-goettingen.de
// Date: 2017

// Use OpenMP if available to parallelize computation
#ifdef _OPENMP 
    #include <omp.h>
#endif

// mex
#include "mex.h"

// stdlib
#include <cmath>
#include <vector>

// Our headers
#include "mexUtil.h"

using namespace std;


// -- Type definitions --
typedef Array1D_t<double> Array1D;
typedef Array1D_t<int> Array1Dint;
typedef Array2D_t<double> Array2D;
typedef Array3D_t<double> Array3D;

// -- Prototypes --
void computeAvgImg(Array3D& imgstack, vector< vector<double> >& avg_img);

void cumulant_ord1(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img);
void cumulant_ord2(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img, Array1Dint& tlag);
void cumulant_ord3(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img, Array1Dint& tlag);
void cumulant_ord4(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img, Array1Dint& tlag);

// --  Global Variables  -- //
mxArray** plhs;  // Access to left hand side is global, so we can output to it in functions


// NOTE: Can be optimized by utilizing redundancy in computations (cache computed values of stack-avg and some multiplications)
// NOTE: Can possibly be optimized by using contiguous memory instead of vector< vector ...> constructs (better caching)

/* % USAGE: [ cum_img ] = autoCumulant2D_cpp( imgstack, order, time_lags  )
% Computes the cumulant 'cum' of order 'order' (max 4) for the 3D (x,y,time)
% image stack 'imgstack'  with time lags 'time_lags'. The required 
% number of specified time lags is #(order-1). */
void mexFunction(int nlhs, mxArray* plhs_arg[], int nrhs, const mxArray* prhs[])
{
//     printf("Thread: %i", omp_get_num_threads()); // Test if openMP is included
    
    // Access to left hand side is global, so we can output to it in functions
    plhs = plhs_arg;
    
    /* Check for proper number of arguments. */
    if(nrhs<3)
    {
        mexErrMsgTxt("Invalid number of input arguments! \n(Image stack, cumulant order, time lags) required.");
    }
    
    /* Map access to input data */
    Array3D imgstack( prhs[0] );
    int order = int(*mxGetPr(prhs[1]) + 0.5);
    
    // Get time lags
    Array1D tmp (prhs[2]); // Get double values from Matlab
    Array1Dint tlag(tmp.nElements);
    for(int i=0; i<tlag.nElements; ++i)
        tlag[i] = static_cast<int>(tmp[i]);
    
    int max_tlag = 0;
    for( int iEl = 0; iEl < tlag.nElements; ++iEl)
    {
        if( tlag[iEl]>max_tlag)
            max_tlag = tlag[iEl];           
    }
    
    if(tlag.nElements != order-1)
    {
        mexErrMsgTxt("Need #(order-1) time lags.");
    }
    
    if(max_tlag >= imgstack.nDepth)
    {
        mexErrMsgTxt("Max. time lag needs to be smaller than the number of frames in the movie!");
    }
    
    
        
    // -- Allocate output memory -- //
    int nDims = 2; // Number of dimensions for output
    mwSize dim_out[2] = { imgstack.nRows, imgstack.nCols };
    plhs[0] = mxCreateNumericArray( nDims, dim_out , mxDOUBLE_CLASS, mxREAL);
    
    Array2D cum_img ( plhs[0] ); // Map access to output
    
    // Compute average image
    vector< vector<double> > avg_img(imgstack.nCols, vector<double>(imgstack.nRows,0) );    
    computeAvgImg(imgstack, avg_img);
    
    // Compute the cumulant
    switch(order)
    {
        case 1: cumulant_ord1(imgstack, avg_img, cum_img); break;
        case 2: cumulant_ord2(imgstack, avg_img, cum_img, tlag); break;
        case 3: cumulant_ord3(imgstack, avg_img, cum_img, tlag); break;
        case 4: cumulant_ord4(imgstack, avg_img, cum_img, tlag); break;
        default: mexErrMsgTxt("Only cumulants up to order 4 supported!");
    }
}


void computeAvgImg(Array3D& imgstack, vector< vector<double> >& avg_img)
{
    // Note: can not be parallelized by frame -> multiple writes on same output pixel
    for(int iFrame = 0; iFrame<imgstack.nDepth; ++iFrame)
        #pragma omp parallel for
        for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
            for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
            {
                avg_img[iCol][iRow] += imgstack(iRow,iCol,iFrame);
            }
    
    #pragma omp parallel for
    for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
        for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
        {
            avg_img[iCol][iRow] /= imgstack.nDepth;
        }
}



void cumulant_ord1(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img)
{
    // <F>  (For joint cumulants this is actually 0)
    #pragma omp parallel for
    for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
        for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
        {
            cum_img(iRow,iCol) = avg_img[iCol][iRow];
        }
}

void cumulant_ord2(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img,  Array1Dint& tlag)
{
    int max_tlag = 0;
    for( int iEl = 0; iEl < tlag.nElements; ++iEl)
    {
        if( tlag[iEl]>max_tlag)
            max_tlag = tlag[iEl];           
    }
    
    // <dF0 * dF1>
    for(int iFrame = 0; iFrame<imgstack.nDepth-max_tlag; ++iFrame)
        #pragma omp parallel for
        for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
            for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
            {
                cum_img(iRow,iCol) += (imgstack(iRow,iCol,iFrame)-avg_img[iCol][iRow]) * (imgstack(iRow,iCol,iFrame+tlag[0])-avg_img[iCol][iRow]);
            }
    
    #pragma omp parallel for
    for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
        for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
        {
            cum_img(iRow,iCol) /= (imgstack.nDepth-max_tlag);
        }
}


void cumulant_ord3(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img,  Array1Dint& tlag)
{
    int max_tlag = 0;
    for( int iEl = 0; iEl < tlag.nElements; ++iEl)
    {
        if( tlag[iEl]>max_tlag)
            max_tlag = tlag[iEl];           
    }
    
    // <dF0 * dF1 * dF2>
    for(int iFrame = 0; iFrame<imgstack.nDepth-max_tlag; ++iFrame)
        #pragma omp parallel for
        for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
            for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
            {
                cum_img(iRow,iCol) += (imgstack(iRow,iCol,iFrame)-avg_img[iCol][iRow]) * (imgstack(iRow,iCol,iFrame+tlag[0])-avg_img[iCol][iRow])  * (imgstack(iRow,iCol,iFrame+tlag[1])-avg_img[iCol][iRow]);
            }
    
    // Divide by number of samples
    #pragma omp parallel for
    for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
        for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
        {
            cum_img(iRow,iCol) /= (imgstack.nDepth-max_tlag);
        }
}


void cumulant_ord4(Array3D& imgstack, vector< vector<double> >& avg_img, Array2D& cum_img,  Array1Dint& tlag)
{
    int max_tlag = 0;
    for( int iEl = 0; iEl < tlag.nElements; ++iEl)
    {
        if( tlag[iEl]>max_tlag)
            max_tlag = tlag[iEl];           
    }
    
    // Reserve memory for intermediate results
    vector< vector<double> > dF0dF3(imgstack.nCols, vector<double>(imgstack.nRows,0) );
    vector< vector<double> > dF1dF2(imgstack.nCols, vector<double>(imgstack.nRows,0) );
    vector< vector<double> > dF0dF2(imgstack.nCols, vector<double>(imgstack.nRows,0) );
    vector< vector<double> > dF1dF3(imgstack.nCols, vector<double>(imgstack.nRows,0) );
    vector< vector<double> > dF0dF1(imgstack.nCols, vector<double>(imgstack.nRows,0) );
    vector< vector<double> > dF2dF3(imgstack.nCols, vector<double>(imgstack.nRows,0) );
    
    // Define helper function computing  Sum( dF(t+tau1)*dF(t+tau2) )
    auto cum2 = [&imgstack, &avg_img, max_tlag] (const unsigned int tau1, const unsigned int tau2, vector< vector<double> >& output)
    {
        for(int iFrame = 0; iFrame<imgstack.nDepth-max_tlag; ++iFrame)
            #pragma omp parallel for
            for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
                for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
                {
                    output[iCol][iRow] += (imgstack(iRow,iCol,iFrame + tau1)-avg_img[iCol][iRow]) * (imgstack(iRow,iCol,iFrame+tau2)-avg_img[iCol][iRow]);
                }
    };
    
    // Compute Sum( dF0 * dF1 * dF2 * dF3 )
    for(int iFrame = 0; iFrame<imgstack.nDepth-max_tlag; ++iFrame)
        #pragma omp parallel for
        for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
            for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
            {
                cum_img(iRow,iCol) += (imgstack(iRow,iCol,iFrame)-avg_img[iCol][iRow]) * (imgstack(iRow,iCol,iFrame+tlag[0])-avg_img[iCol][iRow])  * (imgstack(iRow,iCol,iFrame+tlag[1])-avg_img[iCol][iRow])* (imgstack(iRow,iCol,iFrame+tlag[2])-avg_img[iCol][iRow]);
            }
    
    // Compute pairwise 2nd order cumulants
    cum2(0,tlag[2], dF0dF3); // dF0 * dF3
    cum2(tlag[0],tlag[1], dF1dF2); // dF1 * dF2
    cum2(0,tlag[1], dF0dF2); // dF0 * dF2
    cum2(tlag[0],tlag[2], dF1dF3); // dF1 * dF3
    cum2(0,tlag[0], dF0dF1); // dF0 * dF1
    cum2(tlag[1],tlag[2], dF2dF3); // dF2 * dF3
    
    
    // Divide by number of samples to get the mean
    #pragma omp parallel for
    for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
        for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
        {
            cum_img(iRow,iCol) /= (imgstack.nDepth-max_tlag);
            dF0dF3[iCol][iRow] /= (imgstack.nDepth-max_tlag);
            dF1dF2[iCol][iRow] /= (imgstack.nDepth-max_tlag);
            dF0dF2[iCol][iRow] /= (imgstack.nDepth-max_tlag);
            dF1dF3[iCol][iRow] /= (imgstack.nDepth-max_tlag);
            dF0dF1[iCol][iRow] /= (imgstack.nDepth-max_tlag);
            dF2dF3[iCol][iRow] /= (imgstack.nDepth-max_tlag);
        }
    
    // Compute 4th order joint cumulant <dF0*dF1*dF2*dF3> - <dF0*dF3>*<dF1*dF2> - <dF0*dF2>*<dF1*dF3> - <dF0*dF1>*<dF2*dF3>
    #pragma omp parallel for
    for(int iCol = 0; iCol<imgstack.nCols; ++iCol)
        for(int iRow = 0; iRow<imgstack.nRows; ++iRow)
        {
            cum_img(iRow,iCol) -=  ( dF0dF3[iCol][iRow]*dF1dF2[iCol][iRow]
                    + dF0dF2[iCol][iRow]*dF1dF3[iCol][iRow]
                    + dF0dF1[iCol][iRow]*dF2dF3[iCol][iRow] );
        }
    
}




