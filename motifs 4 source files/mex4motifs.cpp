// mex4motifs.cpp : Defines the entry point for the DLL application.

/****************************************************************************/
#include "stdafx.h"
#include "Evaluate4motifs.h"
#include <mex.h>

void mexFunction(
		int nlhs,              // Number of left hand side (output) arguments
		mxArray *plhs[],       // Array of left hand side arguments
		int nrhs,              // Number of right hand side (input) arguments
		const mxArray *prhs[]  // Array of right hand side arguments
	)
	{
		SetInputParameters(nlhs,plhs,nrhs,prhs);		
		PrepareToRun(nlhs,plhs,nrhs,prhs);
		PerformCalculations(nlhs,plhs,nrhs,prhs);
		FreeResourses(nlhs,plhs,nrhs,prhs);
	}