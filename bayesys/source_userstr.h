/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 *            Bayesian Inference / Massive Inference
 *
 * Filename:  userstr.h
 *
 * Purpose:   Define structures.
 *            Structure members flagged *SYS* are used by the BayeSys3 system.
 *            Other "USER" members can be added for individual applications.
 *
 * History:   JS    25 Apr 2001 - 31 Oct 2002
 *-----------------------------------------------------------------------------
 */
#ifndef USERSTRH
#define USERSTRH

typedef struct         // GENERAL INFORMATION AND WORKSPACE
{
// USER likelihood and monitoring information will usually follow ....
  int       Ndata1;    // I    # data
  int       Nx;     // I    # object cells
  int       Ny;     // I    # object cells
  int       Ncell;     // I    # object cells
  int       Nsample;   //   O  # output ensembles
  int       Nparms;
  int       Npeaks;
  int		MaxPoints;
  int		MinPoints;
  double    Pbeam;
  double    Ddsum;
  double	Sigma;
  double    Xscale;     // I    X-scale
  double    Natoms;    //   O  <# atoms>
  double    Natoms2;   //   O  <# atoms squared>
  double*   Mockbar;   //   O  <mock data>                        [Ndata]
  double*   Object;    //   O  Prob(object +ve)                   [Ncell]
  double*   Mock1;     //   O  Dirty Mock average                 [Ncell]
  double*   Objbar;    //   O  Object average                     [Ncell]
  double*   Objdev;    //   O  Object std.dev.                    [Ncell]
  double*	Parmbar;
  double*	Parmdev;
  double*	Dmap;
  double*	Dbeam;
    double* Power; //NEW VARIABLE
  int		Display1;   // I    First display dimension
  int		Display2;   // I    Second display dimension
  char	  PGDevice[200]; // I	   PGPplot graphics device
  char	  XLabel[200];  // I	   XLabel for display
  char	  YLabel[200];  // I	   YLabel for display
  char	  Title[200];   // I	   Plot title
  char	  Resfil[200];  // I	   results file
} UserCommonStr;

typedef struct         // MEMBER OF ENSEMBLE
{
/// USER likelihood information for individual members can follow ....
  double *Mock1;
} UserObjectStr;

#endif
