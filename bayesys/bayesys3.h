//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//            Bayesian Inference
// 
// Filename:  bayesys3.h
// 
// Purpose:   Header for bayesys3.c
// 
// History:   JS    28 Jan 2002 - 18 Sep 2003
//-----------------------------------------------------------------------------
#ifndef BAYESYS3
#define BAYESYS3

/**************/
/* Structures */
/**************/

typedef struct          // GENERAL INFORMATION AND WORKSPACE
{
// BayeSys3 external 
  int       Ndim;       // I   Hypercube dimension = # variables in atom
  int       MinAtoms;   // I   Min #atoms >= 1
  int       MaxAtoms;   // I   Max #atoms >= MinAtoms (or 0 for unlimited)
  double    Alpha;      // I   Hyperparameter for #atoms (+ve or 0 or -ve)
  int       ENSEMBLE;   // I   # sample objects in ensemble
  int       Method;     // I   Method of algorithm
  double    Rate;       // I   Annealing rate
  int       Iseed;      // I O Random seed, +ve = fixed, -ve = time seed
  double    cool;       //(I)O Annealing (burnin<1, Bayes=1, maxLhood=inf)
  double    Evidence;   //   O log[e] Pr(Data)
  double    Information;//   O -log"Volume" = INT post log(post/prior) dx
  int       Nsystem;    //   O Elapsed iterates
  double    CPU;        //   O CPU in units of 1-atom Del|Try|Ins changes
  double    Success;    //   O # successes (as internally defined)
                        //       (100% efficiency would be Success ~ CPU)
  int       Valency;    // I   0, or MassInf valency

// MassInf external (referenced by BayeSys3 and used if Valency > 0)
  int       MassInf;    // I   Method of analysing Gaussian or Poisson data
                        //     0=monkeys, 1=positive, 2=pos/neg, 3=gauss,
                        //     add 10 to allow flux=0: add 100 if Poisson data
  int       Ndata;      // I   # data
  double*   Data;       // I   Data                                     [Ndata]
  double*   Acc;        // I   1/Gaussian_sigma, or Poisson background  [Ndata]
  double    ProbON;     // I   Prob(individual flux != 0), > 0
  double    FluxUnit0;  // I O Hyperparameter for unit of flux

// BayeSys3 internal
  int       Nbits;      //   O # bits per word
  unsigned* offset;     //  (O)Cube-to-Label mapping                     [Ndim]
  int*      permute;    //  (O)Cube-to-Label mapping                     [Ndim]
  unsigned  Rand[4];    //   O Random generator

// MassInf internal (referenced by BayeSys3 and used if Poisson data)
  int*      Counts;     //  (O)Annealed data counts                     [Ndata]

// USER monitoring information
  void*     UserCommon;
} CommonStr;

typedef struct          // MEMBER OBJECT OF ENSEMBLE
{
// BayeSys3 external 
  double    Lhood;      //   O log Prob(Data|object)
  int       Natoms;     //   O # atoms in object (used dimension of Cube)
  double**  Cubes;      //   O Atomic output           [Natoms][Ndim+Valency+1]
                        //     Ndim cube position, Valency fluxes, log(width)
                        //     De-allocated before exiting BayeSys3

// MassInf external (referenced by BayeSys3 and used there iff Valency > 0)
  double*   Mock;       //  (O)Workspace for Mock data                  [Ndata]
                        //     De-allocated before exiting BayeSys3
  double    FluxUnit;   //   O Unit of flux

// BayeSys3 internal
  int       Nstore;         // # atoms allocated (physical dim of Cube)
  unsigned* xLabel;         // Label control                             [Ndim]
  unsigned* yLabel;         // Label control                             [Ndim]
  unsigned* xOrigin;        // Label control                             [Ndim]
  unsigned* yOrigin;        // Label control                             [Ndim]
  unsigned* xTry;           // Label control                             [Ndim]
  unsigned* yTry;           // Label control                             [Ndim]
  unsigned* work;           // Cube-to-Label workspace                   [Ndim]
  int       clock;          // parallel operation
  int       flowcheck;      // parallel operation
  unsigned  Rand[4];        // Random generator

// MassInf internal (referenced by BayeSys3 and used there if Valency > 0)
  int       reset;          // Reset inserted MassInf fluxes?
  double*   g1;             // grad chisquared                        [Valency]
  double*   g2;             // grad chisquared                        [Valency]
  double*   A11;            // grad grad chisquared                   [Valency]
  double*   A12;            // grad grad chisquared                   [Valency]
  double*   A22;            // grad grad chisquared                   [Valency]
  int*      Xindex;         // Index to cross-terms                   [Valency]
  int*      nbits;          // Fragment numbers                       [Valency]
  int*      nbitx;          // Fragment numbers                       [Valency]
  int*      ibits;          // Fragment identifiers                   [<=Ndata]
  int*      ibitx;          // Fragment identifiers                   [<=Ndata]
  double*   zbits;          // Fragment quantities                    [<=Ndata]
  double*   zbitx;          // Fragment quantities                    [<=Ndata]
  double*   Foot;           // Footprint of unit flux                   [Ndata]
  int*      flags;          // Valency identifiers                      [Ndata]

// USER information for individual objects
  void*     UserObject;
} ObjectStr;

/**************/
/* Procedures */
/**************/
#ifdef __cplusplus
extern "C" {
#endif

                           // Full Bayesian calculation ended by UserMonitor
extern int BayeSys3(       //   O  +ve = UserMonitor return code, -ve = abort
      CommonStr* Common,   // I O  General information
      ObjectStr* Objects); //   O  ENSEMBLE of sample objects

                           
extern void BayeShape(     // Transform cube position to other shape
      double*    Coord,    //   O  Coordinates
      double*    Cube,     // I    Cube (or part of)
      int        N,        // I    Dimension of (part of) Cube
      int        Shape);   // I    Choice of coordinate shape

                           // Set empty object with 0 atoms
extern int UserEmpty(      //   O  >=0 is OK, -ve is error abort
      double*    Lhood,    //   O  loglikelihood
      CommonStr* Common,   // I    General information
      ObjectStr* Object);  //   O  Sample object

                           // Try one new atom
extern int UserTry1(       //   O  +ve = OK, 0 = DO NOT USE, -ve = error abort
      double*    Ltry,     //   O  DELTA(logLikelihood) from 1 trial atom
      CommonStr* Common,   // I    General information
      ObjectStr* Object);  // I    Sample object (DO NOT UPDATE Lhood)

                           // Try two new atoms
extern int UserTry2(       //   O  +ve = OK, 0 = DO NOT USE, -ve = error abort
      double*    Ltry1,    //   O  DELTA(logLikelihood) from 1st trial atom
      double*    Ltry2,    //   O  DELTA(logLikelihood) from both trial atoms
      CommonStr* Common,   // I    General information
      ObjectStr* Object);  // I    Sample object (DO NOT UPDATE Lhood)

                           // Insert one atom and update
extern int UserInsert1(    //   O  >=0 is OK, -ve is error abort
      double*    Lhood,    //   O  DELTA(loglikelihood)
      CommonStr* Common,   // I    General information
      ObjectStr* Object);  // I O  Sample object

                           // Insert two atoms and update
extern int UserInsert2(    //   O  >=0 is OK, -ve is error abort
      double*    Lhood,    //   O  DELTA(loglikelihood)
      CommonStr* Common,   // I    General information
      ObjectStr* Object);  // I O  Sample object

                           // Delete one atom and update
extern int UserDelete1(    //   O  >=0 is OK, -ve is error abort
      double*    Lhood,    //   O  DELTA(loglikelihood)
      CommonStr* Common,   // I    General information
      ObjectStr* Object);  // I O  Sample object

                           // Mock data from unit flux for MassInf
extern int UserFoot(       //   O  +ve = OK, 0 = DO NOT USE, -ve = error abort
      double*    Cube,     // I    Atom coordinates
      CommonStr* Common,   // I    General information
      int*       ibits,    //   O  Fragment coords
      double*    zbits,    //   O  Fragment quantities
      int*       nbits);   //   O  # fragments >= 0, SUM(nbits) <= Ndata

                           // Inspect progress and collect statistics
extern int UserMonitor(    //   O  0 = continue, +ve = finish, -ve = abort
      CommonStr* Common,   // I    General information
      ObjectStr* Objects); // I    ENSEMBLE of sample objects

#ifdef __cplusplus
};
#endif

/*************/
/* Constants */
/*************/
#undef  E_BAYESYS_PARMS
#define E_BAYESYS_PARMS   -220  // Wrong input parameters
#undef  E_BAYESYS_SYSERR
#define E_BAYESYS_SYSERR  -221  // Global error (atom insertion failure)
#undef  E_MASSINF_PARMS
#define E_MASSINF_PARMS   -210  // Wrong input parameters
#undef  E_MASSINF_OVERLAP
#define E_MASSINF_OVERLAP -211  // Valency footprints overlap wrongly
#undef  E_MASSINF_NBITS
#define E_MASSINF_NBITS   -212  // Too many fragments
#undef  E_MASSINF_RANGE
#define E_MASSINF_RANGE   -213  // Fragment outside data range
#undef  E_MASSINF_DATA
#define E_MASSINF_DATA    -214  // Negative count data supplied
#undef  E_MASSINF_COUNTS
#define E_MASSINF_COUNTS  -215  // Count data too large threatening overflow
#undef  E_MALLOC
#define E_MALLOC          -130  // Can't allocate memory

#endif
