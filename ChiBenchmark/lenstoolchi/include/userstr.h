//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//            Bayesian Inference / Massive Inference
// 
// Filename:  userstr.h
// 
// Purpose:   Define user structures.
//
// History:   JS    25 Apr 2001 - 4 Feb 2003
//-----------------------------------------------------------------------------
//
#ifndef USERSTRH
#define USERSTRH

typedef struct          // COMMON PARAMETERS AND STATISTICS
{
  long int  Nsample;    //   O # output ensembles
  double    atoms;      //   O <# atoms>
  long int  Nchi2;	    // # of calls to o_chi()
  double    *err;       // cumulated sum(x^2) for each parameter
  double    *avg;       // cumulated sum(x) for each parameter
  double    sum;        // cumulated sum(var/mean/mean) for all parameters
  long int  Ndim;       // total number of free parameter
} UserCommonStr;

/* 
typedef struct          // INDIVIDUAL PARAMETERS
{
  .....
} UserObjectStr;
*/

#endif
