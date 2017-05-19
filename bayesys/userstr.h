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
  int       Ncell;      // I   # object cells
  int       Nsample;    //   O # output ensembles
  double    atoms;      //   O <# atoms>
  double*   Mockbar;    //   O <mock data>                        [Ndata]
  double*   PrPos;      //   O Prob(object +ve)                   [Ncell]
  double*   Objbar;     //   O Object average                     [Ncell]
  double*   Objdev;     //   O Object std.dev.                    [Ncell]
} UserCommonStr;

/* 
typedef struct          // INDIVIDUAL PARAMETERS
{
  .....
} UserObjectStr;
*/

#endif
