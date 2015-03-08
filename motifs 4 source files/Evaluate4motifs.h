#ifndef EVALUATE_4MOTIFS
#define EVALUATE_4MOTIFS
	#include <mex.h>	
#include "header.h"
// External Functions:
	#ifdef __cplusplus
		extern "C" {
	#endif			

	// functions, called from mexFunction
	void	SetInputParameters(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
	void	PrepareToRun(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
	void	PerformCalculations(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);
	void	FreeResourses(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]);

	#ifdef __cplusplus
		}
	#endif	
#endif	// EVALUATE_4MOTIFS