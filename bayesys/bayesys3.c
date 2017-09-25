//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//             Bayesian Inference
//
// Filename:   bayesys3.c
//
// Purpose:    Obtain sample objects from posterior atomic distribution.
//
// Dedication: To my intellectual ancestors  the late Edwin T Jaynes,  and
//             Steve Gull, to my descendant Sibusiso Sibisi,  and the many
//             colleagues and friends  over the  past quarter-century  who
//             have inspired and encouraged the development of these ideas.
//
//             John Skilling, Kenmare, Ireland, November 2003
//             email: skilling@eircom.net
//=============================================================================
/*
    Copyright (c) 1999-2003 Maximum Entropy Data Consultants Ltd,
                            114c Milton Road, Cambridge CB4 1XE, England

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "license.txt"
*/
//=============================================================================
//
//  An OBJECT is a combination of a-priori-equivalent "atoms":
//  statisticians call this construction a "mixture model".
//
//  ATOMS
//  The number of atoms N in an object can range from MIN to MAX.
//  MIN >= 1 is required to avoid the algorithmically special null object N=0.
//  Generally, MAX >= MIN is required (with MAX=MIN allowed), but
//  MAX=infinity (absence of limit) can be specified with MAX=0.
//  [ Technically, this coding allows at most 1 atom per identifiable
//    location, so "infinity" means #(locations), perhaps 2^32 . ]
//
//  There are three standard priors Pr(N), selected by sign(Alpha).
//  Alpha = 0.
//      The prior is UNIFORM in MIN <= N <= MAX (and MAX must be finite).
//  Alpha > 0.
//      The prior is POISSON or BINOMIAL.
//      N  ~  Alpha +- sqrt(Alpha), subject to MIN and MAX
//  Alpha < 0.
//      The prior is GEOMETRIC.
//      N  ~  |Alpha| +- |Alpha|, subject to MIN and MAX
//
//  COORDINATES
//  The BayeSys prior for locating an atom in space is flat over the hypercube
//  (0,1)^Ndim  of dimension  Ndim .  In fact, the hypercube is treated as
//  wraparound-continuous, so is more properly described as a "hyper-torus".
//  BayeSys gives you coordinate vectors "double Cube[Ndim]" with values in
//  (0,1), but you can transform these into other "Coord" if you wish.
//
//  The BayeShape procedure lets you modify some or all of your coordinates
//  into nine alternative configurations:
//
//   Shape  Description           Prob(Coord[i])                    Range of i
//
//     0   Permutation         uniform on integers Perm(0...N-1)     0...N-1
//     1   +ve orthant         exp(-x) in x > 0                      0...N-1
//     2   Simplex volume      uniform in SUM(x) < 1                 0...N-1
//     3   Simplex surface     uniform on SUM(x) = 1                  0...N
//     4   Ordered             uniform in 0 < x[0] < x[1] < ... < 1  0...N-1
//     5   Bell                Normal(0,1) in -inf < x < inf         0...N-1
//     6   Sphere volume       uniform in SUM(x^2) < 1               0...N-1
//     7   Hemisphere surface  uniform on SUM(x^2) = 1, x[N]>0        0...N
//     8   Sphere surface      uniform on SUM(x^2) = 1                0...N
//
//  For the interior of an awkward shape, allow a circumscribing volume,
//  but return the value 0 from UserTry1 and UserTry2 (instead of the
//  usual positive acceptance code) for any location outside the shape.
//  The ensemble will keep within the domain, because BayeSys only accepts
//  points for which you give strictly positive return codes.
//
//  Thus the unit disc  x^2 + y^2 < 1 could be programmed in several ways:
// (a) Use BayeShape with Shape=6 which gives the N=2 disc interior directly.
// (b) Use N=2 hypercube (= unit square) but expand it to (-1,1) by setting
//        x = 2*Cube[0] - 1,  y = 2*Cube[1] - 1
//     and reject any point outside the disc by returning 0 from UserTry1
//     and UserTry2.
// (c) Use hypercube to yield a unit square, but spread the prior measure
//     uniformly through the unit disc with your own transformation, e.g.
//        radius^2 = Cube[0],  angle = 2 PI Cube[1] .
//
//  FLUXES
//  The Ndim coordinates can be supplemented by  Valency  intensities or
//  "fluxes" z, being attributes for which the joint distribution
//                  Prior(z).Likelihood(D|z,...)
//  is integrable and can be sampled from.
//  The MassInf library incorporated in BayeSys provides Flux procedures for
//  its priors and linear data.
//
//  DISPLAY
//  The coordinates and intensities are supplemented by a guidance width, being
//      log(fraction of hypercube volume that atom might plausibly range over).
//  This may help you to produce smooth displays from the atomic objects.
//
//  The number of atoms in an object is supplied to you as Natoms and
//  each atom's attributes are supplied to you consecutively as
//    double Cube[ 0,1,...,Ndim-1,  Ndim,...,Ndim+Valency-1,  Ndim+Valency ]
//              coordinates in (0,1)     intensities            log(width)
//=============================================================================
//
//  ENSEMBLE
//  The program uses an ensemble of sample objects.  Using several objects is
//  usually recommended, partly because objects that seem to be getting stuck
//  are overwritten by more successful ones which reduces the risk of failure,
//  and partly because the geometrical exploration engines only work with
//  several objects.
//
//=============================================================================
//
//  METHOD
//  Internally in the program, hypercube coordinates are mapped to an
//  extended-integer label whose range fills the hypercube but which preserves
//  a degree of locality: small changes in this integer will necessarily
//  correspond to small changes in coordinates.
//  The user can use the Method parameter to control the style of this mapping,
//  and the operation of the various internal engines that control the
//  evolution of the ensemble.
//
//     Method = 0    is the simplest algorithm, mapping hypercube coordinates
//                   onto extended-integer position labels by simple raster,
//                   and using the "LifeStory1" diffusion engine to create,
//                   destroy, and move atoms.
//                   This is very likely to be enhanced by switching on
//                   various bits of the Method integer.
//
// if( Method & 1 ), hypercube coordinates are mapped to position labels along
//                   a space-filling Hilbert curve.
//                   This is generally recommended, but unnecessary if Ndim=1.
//
// if( Method & 2 ), the "LifeStory" diffusion engine will include interactions
//                   with atoms' neighbours, otherwise not.
//                   This is generally recommended for its extra power.
//
// if( Method & 4 ), the algorithm includes the Chameleon1 engine, which lets
//                   atoms jump from one ensemble object to another, without
//                   changing position.
//                   This is ineffective if there is only one ensemble object.
//
// if( Method & 8 ), the algorithm includes the Chameleon2 engine, which lets
//                   atoms from different ensemble objects exchange position.
//                   This is ineffective if there is only one ensemble object.
//
// if( Method & 16), the algorithm will include the Leapfrog1 engine, which
//                   inverts atom positions with respect to another atom from
//                   the same or a different ensemble object.
//                   This is designed to assist when the posterior distribution
//                   in more than one dimension is highly non-spherical.
//
// if( Method & 32), the algorithm will include the Leapfrog2 engine, which
//                   reflects atom positions with respect to two other atoms
//                   from the same or different ensemble objects.
//                   This is designed to assist when the posterior distribution
//                   in more than one dimension is highly non-spherical.
//
// if( Method & 64), the algorithm will include the GuidedWalk engine, which
//                   moves atoms along directions parallel to displacements
//                   between neighbours, thus following the local shape of
//                   the posterior distribution.
//                   This may supersede the Leapfrog engines.
//
//  ________________________________________________________________________
// |                                                                        |
// | GuidedWalk Leapfrog2 Leapfrog1 Chameleon2 Chameleon1 LifeStory2 Hilbert| 1
// |    off        off       off       off       off      LifeStory1 raster | 0
// |________________________________________________________________________|
//       64        32        16         8         4          2         1 Method
//
//  In the interest of overall efficiency, the algorithm tries to equalise the
//  computation time between its engines so that even if all the optional
//  engines did nothing, the computation cost of including them would be
//  limited to the number of engines (currently a factor of 6).
//
//  Author generally recommends switching everything in by
//                   Method = 127 (equivalently -1) .
//
//=============================================================================
//
// The "Massive Inference" (MassInf) option is provided for applications where
// atoms have intensities or "fluxes" about which the data are linear.
//
// The PRIOR for the number of atoms and for their location is as for BayeSys3.
//
// Prior on flux z of atom is Pr(z) = ProbON * P(z) + (1 - ProbON) * delta(z) ,
// i.e. each flux is expected to be distributed as P(z) with probability ProbON
//      otherwise it is switched off with zero value.
// Common->ProbON supplies ProbON, and the "units" decimal digit of the switch
// Common->MassInf defines the shape P(z) of the prior according to
//
// (0) "monkeys"              P(z) = delta(z-q),   i.e. z = q = constant
//
// (1) "positive"             P(z) = exp(-z/q) / q    in z > 0
//
// (2) "positive/negative"    P(z) = exp(-|z|/q) / 2q
//                                         2     2               2
// (3) "Gaussian"             P(z) = exp(-z / 2 q ) / sqrt(2 pi q )
//
// When flux can be positive or negative, the "positive/negative" prior is
// more tolerant of dynamic range than is the Gaussian prior.
// In each case the flux unit q is fixed at the positive value supplied
// in Common->FluxUnit0, or (if <= 0) kept close to its most probable value.
//
// An atom can have Valency (=1,2,...) independent fluxes, provided
// (a) all atoms have the same valency,
// (b) the flux unit is the same for each,
// (c) the mock-data footprints of the different valencies of an atom
//     do not overlap (no common cells), and
// (d) (if using LifeStory2 engine) for any pair of atoms, each footprint
//     of one overlaps no more than one valency footprint of the other.
//
//
// To specify the LIKELIHOOD (used as its log), define
//
//                    f     = object = sum of atoms
//
//                    x[j]  = positional coordinate of atom j
//
//                    z[j]  = flux(es) of atom j
//
//                    Footprint(x) = mock data from unit flux at x
//
//                    Mock  = SUM[j] z[j] Footprint(x[j]) = mock data of object
//
// Likeliood Pr(Data | f) can be one of
//
// (0xx) "chisquared", derived from
//
//                    Data  = signal,
//
//                    Acc   = 1/sigma, can be 0,   (sigma = standard deviation)
//
//                    M     = # measured data, for which Acc[k] > 0.
//                               M                         2
//                    Z     = PRODUCT[k] sqrt(2 PI sigma[k] )
//                                                      2          2
//                    Chisq = SUM[k] (Mock[k] - Data[k]) / sigma[k]
//                                    -1   -Chisq / 2
//                    Pr(Data | f) = Z    e
//
// (1xx) "Poisson", for positive problems with "monkey" or "positive" prior;
//
//                    Data  = counts, can be floating-point but must be >= 0.0
//
//                    Acc   = extraneous "background", which may be small but
//                            must be strictly +ve wherever Data > 0.0
//                                               -F[k]     D[k]
//                    Pr(Data | f) = PRODUCT[k] e      F[k]    / D[k]!
//
//                    where  F = Mock + Acc, D = Data + Acc,  x! = GAMMA(x+1) :
//                    this form of likelihood avoids singularity when active
//                    cells happen to be empty of mock data, and reduces to
//                    standard Poisson when the background -> 0
//
// according to the "hundreds" digit of the switch Common->MassInf.
//
// For example, Common->MassInf =  1  is Gaussian data with "positive" prior;
//              Common->MassInf = 101  is Poisson data with "positive" prior.
//=============================================================================
// History:
// MassInf1 v1.12-2.36  1998-2000
// BayeSys1 code        1999-2000
// BayeSys2 v1.02-1.05  1 May - 10 Sep 2001
// BayeSys3 v1.01-1.32  2 Jan 2002 - 4 Feb 2003
//          v2.00      10 Feb 2003
//          v2.01      11 Jul 2003  Fixed-q0 option to fix FluxUnit0 in MassInf
//          v2.02      19 Jul 2003  MassInf options extended to allow flux=0
//          v2.10      12 Aug 2003  MassInf needs FluxUnit0 from user
//          v2.20      14 Aug 2003  Remove any overlay in MassInf user's "bits"
//          v3.00      20 Aug 2003  Poisson data as MassInf option
//          v3.01      12 Sep 2003  MassInf control of flux=0 by ProbON
//          v3.02      18 Sep 2003  BayeShape coordinate options
//          v3.03      11 Oct 2003  "Hilbert" name used instead of Peano
//          v3.04      20 Oct 2003  Display width of atom roughy halved
//          v3.05      27 Oct 2003  BayeShape permutation option
//          v3.10      29 Nov 2003  Free software under GNU LGPL
//          v3.11       3 Feb 2004  Control uses <#copies>.  Better Evid & Info
//-----------------------------------------------------------------------------
//
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "bayesys3.h"
#include "random.h"
#include "hilbert.h"

#undef  PARALLEL
#define PARALLEL 0   // # slaved parallel processors [0 = none]
#undef  DEBUG
#define DEBUG    0   // check Cube/Link consistency  [0 = off]

/***********************/
/* Internal structures */
/***********************/
typedef struct
{
    int            i;        // 1st object
    int            j;        // [2nd object]
    int            k;        // [3rd object]
    signed char    engine;   // operation type
    signed char    iType;    // type of 1st object, -1=WRITE, 0=Unused, 1=read
    signed char    jType;    // type of 2nd object, -1=WRITE, 0=Unused, 1=read
    signed char    kType;    // type of 3rd object, -1=WRITE, 0=Unused, 1=read
    int            iEarly;   // min start time # for 1st object i
    int            jEarly;   // min start time # for 2nd object j
    int            kEarly;   // min start time # for 3rd object k
    int            iLate;    // max start time # for 1st object i
    int            jLate;    // max start time # for 2nd object j
    int            kLate;    // max start time # for 3rd object k
    int            slave;    // assigned processor
    int            CPU;      // CPU time in this operation
    int            Success;  // # successes in this operation
} OperStr;                   // Operation

/*********************************/
/* Internal constants and macros */
/*********************************/
//
//  CALLOC(p,n,t)   allocates vector p[0:n-1] of type t
//  REALLOC(p,n,t)  (re-)allocates vector p[0:n-1] of type t
//  FREE(p)         frees CALLOC or REALLOC or NULL vector p[0:*], sets p=NULL
//  CALL            catches negative error codes
//  PLUS(x,y)       log(exp(x)+exp(y))
//  SWAP(x,y)       exchange
//
#undef  CALLOC
#define CALLOC(p,n,t) {p=NULL;\
if((n)>0&&!(p=(t*)calloc((size_t)(n),sizeof(t))))\
{CALLvalue=E_MALLOC;goto Exit;}/*printf("%p %d\n",p,(size_t)(n)*sizeof(t));*/}
#undef  REALLOC
#define REALLOC(p,n,t) {/*printf("%p -1\n",p);*/\
if((n)>0&&!(p=(t*)realloc(p,(size_t)(n)*sizeof(t))))\
{CALLvalue=E_MALLOC;goto Exit;}/*printf("%p %d\n",p,(size_t)(n)*sizeof(t));*/}
#undef  FREE
#define FREE(p) {if(p){/*printf("%p -1\n",p);*/(void)free((void*)p);} p=NULL;}
#undef  CALL
#define CALL(x) {if( (CALLvalue = (x)) < 0 ) goto Exit;}
#undef  PLUS
#define PLUS(x,y) ((x)>(y)?(x)+log(1.+exp((y)-(x))):(y)+log(1.+exp((x)-(y))))
#undef  SWAP
#define SWAP(a,b) {swap = a; a = b; b = swap;}

/***********************/
/* Internal prototypes */
/***********************/
// BayeSys3 main procedures
static void   ShapeSolve  (int, double, double, double*, double*);
static double ShapeCumul  (int, double, double*);
static void   ShapeIndex  (int, double*, double*);
static int    BayesAlloc  (CommonStr*, ObjectStr*, Node*, Node*);
static int    BayesInit   (CommonStr*, ObjectStr*);
static int    PriorInit   (CommonStr*, ObjectStr*, Node*, int, double*);
static int    Control     (CommonStr*, ObjectStr*, int, double*, double*);
static double Copies      (double*, int);
static int    Anneal      (CommonStr*, ObjectStr*, int, double);
static void   Sort2       (int*, double*, int);
static int    FillOcean   (CommonStr*, ObjectStr*, Node*, int);
static int    MCMCengines (CommonStr*, ObjectStr*, Node*);
static void   BayesFree   (CommonStr*, ObjectStr*, Node*, Node*);

// Engines
static int    DoOperations(OperStr*, CommonStr*, ObjectStr*, Node*, int);
static int    Do1operation(OperStr*, CommonStr*, ObjectStr*, Node*, int);
static int    SetPrior    (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    CopyObject  (OperStr*, CommonStr*, ObjectStr*);
static int    MCMCsetup   (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    LifeStory1  (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    LifeStory2  (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    Chameleon1  (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    Chameleon2  (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    Leapfrog1   (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    Leapfrog2   (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    GuidedWalk  (OperStr*, CommonStr*, ObjectStr*, Node*);
static int    MCMCexit    (OperStr*, CommonStr*, ObjectStr*, Node*);

// Event library
static double Birth       (CommonStr*, int);
static double Death       (CommonStr*, int);
static int    SafeCube    (CommonStr*, ObjectStr*);

// Atom control
static int    Empty       (CommonStr*, ObjectStr*);
static int    Try1        (double*, CommonStr*, ObjectStr*);
static int    Try2        (double*, double*, CommonStr*, ObjectStr*);
static int    Insert1     (CommonStr*, ObjectStr*, Node*);
static int    Insert2     (CommonStr*, ObjectStr*, Node*);
static int    Delete1     (int,double*,double*, CommonStr*, ObjectStr*, Node*);

// Label library
static void   Topology    (CommonStr*);
static void   CubetoLabel (unsigned*, ObjectStr*, CommonStr*);
static void   LabeltoCube (ObjectStr*, unsigned*, CommonStr*);

static void   RanLabel    (unsigned*, int, unsigned*);
static void   CopyLabel   (unsigned*, unsigned*, int);
static void   NewLabel    (unsigned*, unsigned*, unsigned*, int,
                           int, int, unsigned*);
static int    Outside     (unsigned*, unsigned*, unsigned*, int);

// MassInf ancillary library
static int    FluxAlloc   (CommonStr*, ObjectStr*);
static int    FluxInit    (CommonStr*, ObjectStr*);
static int    FluxCalib0  (CommonStr*, ObjectStr*, Node*);
static void   FluxCalib   (CommonStr*, ObjectStr*);
static void   FluxFree    (CommonStr*, ObjectStr*);

// MassInf outer procedures
static void   FluxEmpty   (double*, CommonStr*, ObjectStr*);
static int    FluxTry1    (double*, CommonStr*, ObjectStr*);
static int    FluxTry2    (double*, double*, CommonStr*, ObjectStr*);
static int    FluxInsert1 (double*, CommonStr*, ObjectStr*);
static int    FluxInsert2 (double*, CommonStr*, ObjectStr*);
static int    FluxDelete1 (double*, CommonStr*, ObjectStr*);

// MassInf application procedures for user Footprints
static int    AtomBits    (double*, CommonStr*, int*, double*, int*, int*);
static int    Overlap1    (int*, int*, int*, int);
static int    Overlap2    (int*, int*, int*, int*, int*, int*, int);
static void   SetIndex    (int, int*, int*, int*, int*, int*, int*);
static void   InsBits     (double*, double*, int*, int*, double*, int);
static void   DelBits     (double*, double*, int*, int*, double*, int);

// MassInf application procedures for quadratic chisquared
static void   SetGrad1    (double*, double*, double*, int, int*, int*, double*,
                           double*, double*);
static void   SetGrad2    (double*, double*, double*, double*, int, int*, int*,
                           double*, int*, int*, double*, double*, double*,
                           double*, double*, double*, int*);

// MassInf application procedures for quadratic probabilities
static double GaussLhood  (int, double*, double*, double*);
static double GaussTry1   (int, double, double, double, int, double*, double*);
static void   GaussTry2   (int, double, double, double, int, double*, double*,
                           double*, double*, double*, int*, double*, double*);
static void   GaussInsert1(int, Rand_t, double, double, double, int, double*,
                           double*, double*);
static void   GaussInsert2(int, Rand_t, double, double, double, int, double*,
                           double*, double*, double*, double*, int*, double*,
                           double*);
static double GaussLhood1 (int, double*, double*, double*, int);
static double GaussLhood2 (int, double*, double*, double*, double*, double*,
                           int*, double*, double*);
// without Valency ....
static double Gauss1Marginal(int, double, double, double, double, double);
static double Gauss2Marginal(int, double, double, double, double, double,
                             double, double, double);
static double Gauss1Sample  (int, Rand_t, double, double, double, double,
                             double);
static void   Gauss2Sample  (int, Rand_t, double, double, double, double,
                             double, double, double, double, double*, double*);
// without ProbON ....
static double gauss1marginal(int, double, double, double, double);
static double gauss2marginal(int, double, double, double, double, double,
                             double, double);
static double gauss1sample  (int, Rand_t, double, double, double, double);
static void   gauss2sample  (int, Rand_t, double, double, double, double,
                             double, double, double, double*, double*);

// MassInf application procedures for Poisson probabilities
static double PoissLhood    (int, double*, double*, double*);
static int    PoissTry1     (int, double, double, double, double*, double*,
                             double*, int*, int, int*, int*, double*, double*);
static int    PoissTry2     (int, double, double, double, double*, double*,
                             double*, int*, double*, int, int*, int*, double*,
                             int*, int*, double*, int*, double*, double*);
static int    PoissInsert1  (int, Rand_t, double, double, double, double*,
                             double*, double*, int*, int, int*, int*, double*,
                             double*);
static int    PoissInsert2  (int, Rand_t, double, double, double, double*,
                             double*, double*, int*, double*, int, int*, int*,
                             double*, int*, int*, double*, int*, double*,
                             double*);
static double PoissLhood1   (double*, double*, double*, int, int*, int*,
                             double*, double*, int);
static double PoissLhood2   (double*, double*, double*, double*, int, int*,
                             int*, double*, int*, int*, double*, int*,
                             double*, double*);
// without Valency ....
static int    Poiss1Marginal(int, double, double, double, double*, double*,
                             double*, int*, int, int*, double*, double*);
static int    Poiss2Marginal(int, double, double, double, double*, double*,
                             double*, int*, double*, int, int*, double*, int,
                             int*, double*, double*, double*);
static int    Poiss1Sample  (int, Rand_t, double, double, double, double*,
                             double*, double*, int*, int, int*, double*,
                             double*);
static int    Poiss2Sample  (int, Rand_t, double, double, double, double*,
                             double*, double*, int*, double*, int, int*,
                             double*, int, int*, double*, double*, double*);
// without ProbON ....
static int    Poisson1      (int, Rand_t, double, double, double*, double*,
                             double*, int*, int, int*, double*, double*);
static int    Poisson2      (int, Rand_t, double, double, double*, double*,
                             double*, int*, double*, int, int*, double*, int,
                             int*, double*, double*, double*);
static double poiss1lhood   (double, double*, double*, double*, int, int*,
                             double*);
static void   poiss2lhood   (double, double, double*, double*, double*,
                             double*, int, int*, double*, int, int*, double*,
                             double*, double*);

//=============================================================================
//
//                      BayeSys3 main procedures
// main
//  |
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
// BayeSys3
//  |   \_____________________________________________________________
//  |           |       |         |           |            |          |
//  |       PriorInit Control MCMCengines  BayesAlloc   BayesFree  FillOcean
//  |            \    Anneal   /     \     BayesInit
//  |             \     |     /       \     /
//  |              DoOperations       Topology
//  |                   |
//  |              Do1operation
//  |      _____________|______________________________________________
//  |     |             |         |           |            |           |
//  |  CopyObject   SetPrior   MCMCsetup  LifeStory1   LifeStory2   MCMCexit
//  |                   |         |       Chameleon1      /|
//  |                   |         |       Chameleon2     / |
//  |                   |         |       Leapfrog1     /  |
//  |                   |         |       Leapfrog2    /   |
//  |                   |         |       GuidedWalk  /    |
//  |                   |         |             \    /     |
//  |                  Empty     Empty         Insert1   Insert2
//  |                 Insert1   Insert1        Delete1    Try2
//  |                  Try1                     Try1
//  |
//  |      Empty     Insert1     Try1     Delete1     Insert2     Try2
//  |        |         | \_________|_\______|_\_________|_\_________|_\___
//  |        |         |           |        |           |           |     |
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//  |        |         |           |        |           |           |     |
//  |    UserEmpty UserInsert1 UserTry1 UserDelete1 UserInsert2 UserTry2  |
//  |                                                                     |
//  |                                                                     |
// UserMonitor                                                        UserFoot
//=============================================================================

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  BayeSys3
//
// Purpose:   Perform Markov chain Monte Carlo algorithm for
//            Bayesian analysis or maximisation over unit hypercube,
//            using an ensemble of one or more objects.
//
// History:   JS BayeSys1      3 Mar 1999 -  3 Nov 2000
//               BayeSys2     25 Apr 2001 - 10 Sep 2001
//               BayeSys3     12 Jan 2002 -  3 Feb 2004
//-----------------------------------------------------------------------------
//
int BayeSys3(            //   O  +ve finish code from UserMonitor (or error)
CommonStr* Common,       // I O  General information
ObjectStr* Objects)      //   O  ENSEMBLE of objects
{
static const double BIG = DBL_MAX / 3.0000;  // for annealing protection
static const int    NOCEAN = 50000;  // max # atoms in Ocean, for atom-widths
    int             NTABLE = 1000;   // current size of evidence table
#undef  NEXTRA
#define NEXTRA 12              /* # likelihoods in prior settings */
    int        ENSEMBLE = Common->ENSEMBLE;
    int        Valency  = Common->Valency;
    double     Lextra[NEXTRA]; // Old likelihoods for stability
    Node*      Links   = NULL; // Linked lists of atom labels
    Node       Ocean[1];       // Linked list of accumulated atoms
    double     cold;           // Old reciprocal temperature
    double     dcool;          // Cooling increment
    double     Lbar;           // <logL>
    double*    Ltable = NULL;  // table of strictly increasing logL values
    double*    ctable = NULL;  // table of strictly increasing cool values
    int*       ntable = NULL;  // table of occupation numbers for averaging
    double     Evid;           // SUM[1..] Ltable[i] * (ctable[i]-ctable[i-1])
    int        j;              // table counter
    int        k;              // object counter
    int        finish;         // exit code from UserMonitor
    double     t;              // temporary
    int        CALLvalue  = 0;
#if DEBUG
  printf("WARNING: DEBUG switched on\n");
#endif

// Check input parameters
    if( Common->ENSEMBLE < 1 )        return E_BAYESYS_PARMS;
    if( Common->MinAtoms < 1 )        return E_BAYESYS_PARMS;
    if( Common->MaxAtoms && Common->MaxAtoms < Common->MinAtoms )
                                      return E_BAYESYS_PARMS;
    if( Common->MaxAtoms < Common->MinAtoms && Common->Alpha == 0.0 )
                                      return E_BAYESYS_PARMS;
    if( Common->Rate <= 0.0 )         return E_BAYESYS_PARMS;
    if( Common->Ndim < 1 )            return E_BAYESYS_PARMS;
    if( Valency )
    {
        if( Common->MassInf !=   0 && Common->MassInf !=   1
         && Common->MassInf !=   2 && Common->MassInf !=   3
         && Common->MassInf != 100 && Common->MassInf != 101 )
                                      return E_MASSINF_PARMS;
        if(   Common->Ndata < 0 )     return E_MASSINF_PARMS;
        if( ! Common->Data )          return E_MASSINF_PARMS;
        if( ! Common->Acc )           return E_MASSINF_PARMS;
        if( Common->ProbON <= 0.0 )   return E_MASSINF_PARMS;
        if( Common->MassInf >= 100 )
            for( k = 0; k < Common->Ndata; k++ )
            {
               if( Common->Acc[k] < 0. || Common->Data[k]+Common->Acc[k] < 0. )
                   return E_MASSINF_DATA;
               if( Common->Acc[k] == 0.0 && Common->Data[k] != 0.0 )
                   return E_MASSINF_DATA;
            }
    }
// Allocate memory
    CALLOC(Ltable, NTABLE, double)
    CALLOC(ctable, NTABLE, double)
    CALLOC(ntable, NTABLE, int)
    CALLOC(Links,  ENSEMBLE, Node)
    CALL( BayesAlloc(Common, Objects, Links, Ocean) )
// Initialise system
    CALL( BayesInit(Common, Objects) )
// Initialise prior objects
    CALL( PriorInit(Common, Objects, Links, NEXTRA, Lextra) )

// Enter at prior (infinite temperature)
    dcool = BIG;
    ntable[0] = j = 0;
    ctable[0] = Ltable[0] = Evid = 0.0;
    do
    {
        Common->Nsystem ++;
// Update Ocean
        CALL( FillOcean(Common, Objects, Ocean, NOCEAN) )
// Controlled cooling
        cold = Common->cool;                       // previous value
        CALL( Control(Common, Objects, NEXTRA, Lextra, &t) )
        Common->cool += t;
        if( Common->cool > cold + 2.0000 * dcool ) // MAX allowable cooling to
            Common->cool = cold + 2.0000 * dcool;  // prevent fast acceleration
// User-limited cooling, at posterior (= 1) or beyond (> 1)
// Reset any nuisance parameters
// Collect statistics and display diagnostics
        CALL( UserMonitor(Common, Objects) )       // should not increase dcool
        finish = CALLvalue;                        // enables user to exit
        if( Common->cool > BIG )                   // avoid overflow
            Common->cool = BIG;
        dcool = Common->cool - cold;
// Cool ensemble by copying (user's nuisance parameters won't be copied)
        CALL( Anneal(Common, Objects, NEXTRA, dcool) )
// Evolve
        CALL( MCMCengines(Common, Objects, Links) )
// Evidence and Information diagnostics
        Lbar = 0.0;
        for( k = 0; k < ENSEMBLE; k++ )
            Lbar += Objects[k].Lhood;
        Lbar /= ENSEMBLE;
        if( j == NTABLE - 1 )
        {
            NTABLE += NTABLE / 2;
            REALLOC(Ltable, NTABLE, double)
            REALLOC(ctable, NTABLE, double)
            REALLOC(ntable, NTABLE, int)
        }
        if( Common->cool > ctable[j] )   // new table value
        {
            j++;
            ctable[j] = Common->cool;
            Ltable[j] = Lbar;
            ntable[j] = 1;
            Evid += (ctable[j] - ctable[j-1]) * Ltable[j];
        }
        else                             // update existing top value
        {
            Evid -= (ctable[j] - ctable[j-1]) * Ltable[j];
            Ltable[j] = (ntable[j] * Ltable[j] + Lbar) / (ntable[j] + 1);
            ntable[j] ++;
            Evid += (ctable[j] - ctable[j-1]) * Ltable[j];
        }
// Prune tables to satisfy known constraint: logL increasing function of cool
// If logL tries to decrease, feed the deficit backwards
        while( j > 1 && Ltable[j-1] >= Ltable[j] )
        {
            Evid -= (ctable[j] - ctable[j-1]) * Ltable[j];
            j--;
            Evid -= (ctable[j] - ctable[j-1]) * Ltable[j];
            ctable[j] = ctable[j+1];
            Ltable[j] = (ntable[j] * Ltable[j] + ntable[j+1] * Ltable[j+1])
                       / (ntable[j] + ntable[j+1]);
            ntable[j] += ntable[j+1];
            Evid += (ctable[j] - ctable[j-1]) * Ltable[j];
        }
        Common->Evidence = Evid;
        Common->Information = ctable[j] * Ltable[j] - Evid;

// Exit?
    } while( ! finish );
    CALLvalue = finish;

Exit:
// Free memory
    BayesFree(Common, Objects, Links, Ocean);
    FREE(Links)
    FREE(ntable)
    FREE(ctable)
    FREE(Ltable)
    return CALLvalue;
#undef NEXTRA
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  BayeShape
//
// Purpose:   Transfer BayeSys hypercube to User's shape.
//
//   Shape  Description           Prob(Coord[i])                    Range of i
//
//     0   Permutation         uniform on integers Perm(0...N-1)     0...N-1
//     1   Positive orthant    exp(-x) in x > 0                      0...N-1
//     2   Simplex volume      uniform in SUM(x) < 1                 0...N-1
//     3   Simplex surface     uniform on SUM(x) = 1                  0...N
//     4   Ordered             uniform in 0 < x[0] < x[1] < ... < 1  0...N-1
//     5   Bell                Normal(0,1) in -inf < x < inf         0...N-1
//     6   Sphere volume       uniform in SUM(x^2) < 1               0...N-1
//     7   Hemisphere surface  uniform on SUM(x^2) = 1, x[N]>0        0...N
//     8   Sphere surface      uniform on SUM(x^2) = 1                0...N
//
// Notes:     The hypercube has 2^N corners, at which there may be singular
//            anisotropy in the exploration.  It is better to keep any such
//            singularities separated and diluted.
//
//     0   Permutation of integers 0,1,...,N-1.
//         No singularity, permutation given by ranked order of Cube[0...N-1].
//
//     1   Positive orthant
//         No singularity.
//         The origin remains unchanged.
//         Other corners are at infinity.
//
//     2   Simplex volume
//         The origin remains unchanged.
//         Its N adjacent corners go to the remaining N vertices x[i]=1
//         of the simplex.
//         The remaining 2^N-N-1 corners are isolated singularities
//         distributed over the simplex roof SUM(x)=1.
//
//     3   Simplex surface
//         The origin moves to the upper vertex x[N]=1.
//         Its N adjacent corners go to the other surface vertices x[i]=1.
//         The remaining 2^N-N-1 corners are isolated singularities
//         distributed over the floor x[N]=0.
//
//     4   Ordered
//         The origin does not move.
//         Its N adjacent corners go to the other vertices (0,..0,0,1,1,..,1).
//         The remaining 2^N-N-1 corners are isolated singularities
//         distributed over the roof x[N-1]=1.
//         This is better than just sorting the x[i] from a hypercube,
//         with its N!-to-1 overlapping.
//
//     5   Bell
//         No singularity.
//         The hypercube centre becomes the origin.
//         All corners project out to infinity.
//
//     6   Sphere volume
//         All corners are isolated singularities on the sphere surface
//         SUM(x^2)=1, distributed with their original cubic symmetry.
//         This is much better than concentrating all the singularity at
//         the centre, as with using polar coordinates.
//
//     7   Hemisphere surface
//         All corners are isolated singularities at the equator x[N]=0,
//         distributed with their original cubic symmetry.
//
//     8   Sphere surface
//         All corners are isolated singularities at the surface's equator
//         x[N]=0, distributed with their original cubic symmetry.
//         The penalty for having this arrangement is that the north and south
//         hemispheres touch everywhere, with even Hilbert points defined as
//         north x[N]>0 and odd Hilbert points defined as south x[N]<0.
//         This will slow exploration because about half of the trial
//         points will be in the wrong hemisphere.  And there will be trouble
//         if an over-intelligent user tries to hide some information of his
//         own in this bit of lowest arithmetical significance.
//         The alternative, of collecting all the singularities together
//         down at the south pole, is worse.
//
// History:   JS        18 Sep 2003, 20 Oct 2003
//-----------------------------------------------------------------------------
//
void BayeShape(
double*    Coord,     //   O  coordinates              [N] or [N+1]
double*    Cube,      // I    hypercube position       [N]
int        N,         // I    dimension
int        Shape)     // I    choice of shape
{
static const double logG1[20] = { // log( (N/2)! )
            0.0000000000000000, -0.1207822376352452,  0.0000000000000000,
            0.2846828704729193,  0.6931471805599453,  1.2009736023470752,
            1.7917594692280550,  2.4537365708424432,  3.1780538303479458,
            3.9578139676187183,  4.7874917427820458,  5.6625620598571409,
            6.5792512120101012,  7.5343642367587300,  8.5251613610654147,
            9.5492672573009951, 10.6046029027452509, 11.6893334207972703,
           12.8018274800814691, 13.9406252194037652};
static const double  Z = (unsigned)(-1) + 1.0;          // 2^32
static const double  logpibytwo   =  0.45158270528945486473;

    unsigned  hemisphere;
    double    p, q, r, rr, s, t, inward, inwardN, outward, scale;
    int       i;

    switch( Shape )
    {
    case 0:
        ShapeIndex(N, Cube, Coord);
        break;

    case 1:
    case 2:
    case 3:
    case 4:
        for( i = 0; i < N; i++ )
            Coord[i] = -log(Cube[i]);
        if( Shape == 1 )
            break;

        r = 0.0;
        for( i = 0; i < N; i++ )
            r += Coord[i];                     // Old radius
        s = t = q = exp(-r/N);
        p = r * q;
        for( i = 1; i < N; i++ )
        {
            t *= p / i;
            s = s * q + t;
        }                                      // Outward cumulant
        s = (s * s < DBL_EPSILON) ? 1.0 - s / N
           : pow(1.0 - s, 1.0 / N);            // New radius
        s /= r;
        for( i = 0; i < N; i++ )               // Rescale
            Coord[i] *= s;
        if( Shape == 2 )
            break;

        t = 1.0;
        for( i = 0; i < N; i++ )
            t -= Coord[i];
        Coord[N] = t;
        if( Shape == 3 )
            break;

        for( i = 1; i <= N; i++ )
            Coord[i] += Coord[i-1];
        break;

    case 5:
    case 6:
    case 7:
    case 8:
        for( i = 0; i < N; i++ )
            Coord[i] = InvNorm(Cube[i]);
        if( Shape == 5 )
            break;

        rr = 0.0;
        for( i = 0; i < N; i++ )
        {
            Coord[i] = t = InvNorm(Cube[i]);
            rr += t * t;
        }
        r = sqrt(rr);
        scale = (N < 20) ? logG1[N] : logGamma(N/2. + 1.);
        scale = sqrt(2.0) / exp(scale / N);
        inwardN = r * scale;
        inward = pow(inwardN, N);
        outward = 1.0 - inward;                    // polar limit
        if( inward * inward > DBL_EPSILON )        // non-polar evaluation
        {
            if( N & 1 )
            {
                q = exp(-rr / (N + 1.0));
                p = q * rr;
                t = q / rr;
                s = 0.0;
                for( i = 1; i < N; i += 2 )
                {
                    t *= p / i;
                    s = s * q + t;
                }
                s = exp(logerf(r, 1) - 0.5 * (rr + logpibytwo)) + r * s;
            }
            else
            {
                t = s = q = exp(- rr / N);
                p = q * rr;
                for( i = 2; i < N; i += 2 )
                {
                    t *= p / i;
                    s = s * q + t;
                }
            }
            outward = s;
            inward = 1.0 - outward;
            inwardN = pow(inward, 1.0 / N);
        }
        if( Shape == 6 )
        {
            scale = (outward * outward < DBL_EPSILON)
                   ? (1.0 - outward / N) / r : inwardN / r;
            for( i = 0; i < N; i++ )
                Coord[i] *= scale;
            break;
        }

        ShapeSolve(N, outward, inwardN, &scale, &Coord[N]);
        scale /= r;
        for( i = 0; i < N; i++ )
            Coord[i] *= scale;
        if( Shape == 7 )
            break;

        hemisphere = 0;
        for( i = 0; i < N; i++ )
            hemisphere ^= (unsigned)(Cube[i] * Z);          // Hilbert parity
        if( hemisphere & 1 )
            Coord[N] = -Coord[N];
        break;

    default:
        for( i = 0; i < N; i++ )
            Coord[i] = Cube[i];
        break;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  ShapeSolve
//
// Purpose:   Latitude phi for which
//                   ShapeCumul(phi) / ShapeCumul(PI/2) = outward
//            Equator is    outward = 0, phi = 0
//            North pole is outward = 1, phi = PI/2
//
// Method:    FROWN(a) = ShapeCumul(a) is concave,
//                       FROWN''(a) < 0,
//                       so if a < aim, tangent meets aim at better lower bound
//            SMILE(b) = - (ShapeCumul(PI/2) - ShapeCumul(b))^(1/N) is convex,
//                       SMILE''(b) > 0,
//                       so if b > aim, tangent meets aim at better upper bound
//            Chop between lower and upper bounds
//
// History:   JS        18 Sep 2003
//-----------------------------------------------------------------------------
static
void ShapeSolve(
int      N,         // I    dimension of surface of sphere
double   outward,   // I    poleward cumulant
double   inwardN,   // I    (equator-ward cumulant)^(1/N)
double*  cosp,      //   O  cos(latitude)
double*  sinp)      //   O  sin(latitude)
{
static const double  pibytwo = 1.57079632679489661923;
    double  frac = 1.0 / N;
    double  a;      // left bound
    double  ya;     // FROWN(a)
    double  ma;     // FROWN'(a) > 0
    double  b;      // right bound
    double  yb;     // SMILE(b)
    double  mb;     // SMILE'(b) > 0
    double  c;      // central estimate
    double  yc;     // ShapeCumul(c)
    double  mc;     // ShapeCumul'(c)
    double  top;    // ShapeCumul(PI/2)
    double  aim;    // required value for FROWN
    double  bim;    // required value for SMILE
    int     iter;   // 5 chops suffice for full accuracy

    if( N == 1 )
    {
        a = outward * pibytwo;
        *cosp = cos(a);
        *sinp = sin(a);
        return;
    }
    if( N == 2 )
    {
        *cosp = sqrt(1.0 - outward * outward);
        *sinp = outward;
        return;
    }

    top = ShapeCumul(N, pibytwo, &mc);
    a = outward * top;
    b = a * a * (N - 2) / 6.0;
    if( b * b <= DBL_EPSILON )
    {
        *sinp = a * (1.0 + b);
        *cosp = sqrt(1.0 - *sinp * *sinp);
        return;
    }
    a = inwardN * pow(N * top, frac);
    b = a * a / (2.0 * N + 4.0);
    if( b * b <= DBL_EPSILON )
    {
        *cosp = a * (1.0 - b);
        *sinp = sqrt(1.0 - *cosp * *cosp);
        return;
    }

    aim = outward * top;
    bim = - pow(top - aim, frac);
    a = 0.0;      ya = 0.0;    ma = 1.0;
    b = pibytwo;  yb = 0.0;    mb = pow(frac, frac);
    for( iter = 0; iter < 5; iter++ )
    {
        c = (a + b) / 2.0;
        yc = ShapeCumul(N, c, &mc);
        if( yc < aim )
        {
            a = c;
            ya = yc;    ma = mc;
            a += (aim - ya) / ma;
        }
        else
        {
            b = c;
            yb = - pow(top - yc, frac);    mb = - mc * frac * yb / (top - yc);
            b += (bim - yb) / mb;
        }
    }
// protected interpolation
    ya = ShapeCumul(N, a, &ma) - aim;
    yb = ShapeCumul(N, b, &mb) - aim;
    if( ya >= yb )
        c = (a + b) / 2.0;                   // unusual
    else
    {
        c = (a * yb - b * ya) / (yb - ya);   // linear
        mc = (mb - ma) / (yb - ya);
        c += 0.5 * mc * (b - c) * (c - a);   // quadratic correction
    }
    if( c < a )  c = a;                      // final protection
    if( c > b )  c = b;                      // final protection
    *cosp = cos(c);
    *sinp = sin(c);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  ShapeCumul
//
// Purpose:     phi       N-1
//            INTEGRAL cos   (theta) d(theta)
//               0
//
// History:   JS        18 Sep 2003
//-----------------------------------------------------------------------------
static
double ShapeCumul(  //   O integral
int     N,          // I   dimension of surface of hemisphere
double  phi,        // I   latitude
double* deriv)      //   O cos(phi)^(N-1)
{
    double  r, s, t, u, sinp, cosp, cos2;
    int     i;

    sinp = sin(phi);
    cosp = cos(phi);
    cos2 = cosp * cosp;
    s = t = u = 1.0;
    for( i = (N & 1) + 2; i < N; i += 2 )
    {
        r = 1.0 / i;
        t *= 1.0 - r;
        u *= cos2;
        s += t * u;
    }
    t *= N - 1.0;
    s *= sinp;
    u *= cosp;
    if( N & 1 )
    {
        s = (phi + cosp * s) / 2.0;
        u *= cosp;
        t /= 2.0;
    }
    *deriv = u;
    return s / t;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  ShapeIndex
//
// Purpose:   Index vector y into increasing order
//            y[p[0]] <= y[p[1]] <= .... <= y[p[n-1]]
//
// History:   JS    27 Apr 1999, 27 Oct 2003
//-----------------------------------------------------------------------------
static
void  ShapeIndex(
int       N,    // I    dimension
double*   y,    // I    vector being indexed by value
double*   p)    //   O  index (integer values)
{
    int    j, k, l, m, q;
    float  r;

    if( N <= 1 )
    {
        p[0] = 0.0;
        return;
    }
    for( j = 0; j < N; j++ )
        p[j] = (double)j;
    l = N / 2;
    k = N - 1;
    while( 1 )
    {
        if( l > 0 )
        {
            q = (int)p[--l];
            r = y[q];
        }
        else
        {
            q = (int)p[k];
            r = y[q];
            p[k--] = p[0];
            if( k == 0 )
            {
                p[0] = (double)q;
                return;
            }
        }
        m = l;
        j = l + l + 1;
        while( j <= k )
        {
            if( j < k && y[(int)p[j]] < y[(int)p[j+1]] )
                j++;
            if( r < y[(int)p[j]] )
            {
                p[m] = p[j];
                m = j;
                j = j + j + 1;
            }
            else
                j = k + 1;
        }
        p[m] = (double)q;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  BayesAlloc
//
// Purpose:   Allocate memory for BayeSys (and MassInf)
//
// History:   JS        16 Oct 2002, 1 Feb 2003, 13 Sep 2003
//-----------------------------------------------------------------------------
static
int BayesAlloc(      //   O  0, or -ve error
CommonStr* Common,   // I(O) general information
ObjectStr* Objects,  //  (O) new ENSEMBLE of objects            [ENSEMBLE]
Node*      Links,    //  (O) new trees                          [ENSEMBLE]
Node*      Ocean)    //  (O) new ocean                                 [1]
{
    int        ENSEMBLE  = Common->ENSEMBLE;
    int        Ndim      = Common->Ndim;
    int        Valency   = Common->Valency;
    int        Nsize     = Ndim + Valency + 1;
    ObjectStr* Object;
    int        j, k;
    int        CALLvalue = 0;

// Allow subsequent CALLOC memory allocation to fail gracefully
    Common->offset  = NULL;
    Common->permute = NULL;

// Allocate
  // Common
    CALLOC(Common->offset,  Ndim, unsigned)
    CALLOC(Common->permute, Ndim, int)
  // Objects
    for( k = 0; k < ENSEMBLE; k++ )
    {
        Object = &Objects[k];
        Object->xLabel  = NULL;
        Object->yLabel  = NULL;
        Object->xOrigin = NULL;
        Object->yOrigin = NULL;
        Object->xTry    = NULL;
        Object->yTry    = NULL;
        Object->work    = NULL;
        Object->Cubes   = NULL;
        CALLOC(Object->xLabel,  Ndim, unsigned)
        CALLOC(Object->yLabel,  Ndim, unsigned)
        CALLOC(Object->xOrigin, Ndim, unsigned)
        CALLOC(Object->yOrigin, Ndim, unsigned)
        CALLOC(Object->xTry,    Ndim, unsigned)
        CALLOC(Object->yTry,    Ndim, unsigned)
        CALLOC(Object->work,    Ndim, unsigned)
        Object->Nstore = Common->MinAtoms + 2;
        CALLOC(Object->Cubes, Object->Nstore, double*)
        for( j = 0; j < Object->Nstore; j++ )
            CALLOC(Object->Cubes[j], Nsize, double)
        CALL( SetLink(Ndim, &Links[k]) )
    }
    if( Valency )
        CALL( FluxAlloc(Common, Objects) )
  // Ocean
    CALL( SetLink(Ndim, Ocean) )
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  BayesInit
//
// Purpose:   Initialise BayeSys (and MassInf)
//
// History:   JS        16 Oct 2002, 1 Feb 2003
//-----------------------------------------------------------------------------
static
int BayesInit(       //   O  0, or -ve error
CommonStr* Common,   // I O  general information
ObjectStr* Objects)  //   O  new ENSEMBLE of objects
{
    unsigned*  Rand      = Common->Rand;
    int        ENSEMBLE  = Common->ENSEMBLE;
    int        Ndim      = Common->Ndim;
    int        Valency   = Common->Valency;
    int        Nsize     = Ndim + Valency + 1;
    ObjectStr* Object;
    int        i, j, k;
    int        CALLvalue = 0;

// Initialise BayeSys
    Common->Nbits = 0;
    for( Rand[0] = 1; Rand[0]; Rand[0] <<= 1 )
        Common->Nbits++;
    CALL( RanInit(Rand, Common->Iseed) )
    Common->Iseed       = CALLvalue;
    Common->cool        = 0.0;
    Common->Evidence    = 0.0;
    Common->Information = 0.0;
    Common->Nsystem     = 0;
    Common->Success     = 0.0;
    Common->CPU         = 0.0;
    for( k = 0; k < ENSEMBLE; k++ )
    {
        Object = &Objects[k];
        j = Ranint(Rand);
        if( j < 0 )
            j = ~j;
        CALL( RanInit(Object->Rand, j) )
        Objects[k].Natoms = 0;
        for( j = 0; j < Object->Nstore; j++ )
            for( i = 0; i < Nsize; i++ )
                Object->Cubes[j][i] = 0.0;
        Object->reset = 1;
    }
    Topology(Common);
    if( Valency )
        CALL( FluxInit(Common, Objects) )
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PriorInit
//
// Purpose:   Initialise prior objects and extra loglikelihood values
//
// History:   JS        16 Oct 2002, 1 Feb 2003
//-----------------------------------------------------------------------------
static
int PriorInit(       //   O  0, or -ve error
CommonStr* Common,   // I O  general information
ObjectStr* Objects,  //   O  ENSEMBLE of objects
Node*      Links,    //  (O) empty on exit
int        NEXTRA,   // I    # extra loglikelihood values
double*    Lextra)   //   O  extra loglikelihood values   [NEXTRA]
{
    int       ENSEMBLE   = Common->ENSEMBLE;
    OperStr*  Operations = NULL; // list of operations   [ENSEMBLE]
    OperStr*  Oper;              // & individual operation
    int       i, j, k;
    int       CALLvalue  = 0;

// Pre-calibrate objects
    if( Common->Valency )
        CALL( FluxCalib0(Common, Objects, Links) )
    CALLOC(Operations, ENSEMBLE, OperStr)
// Find O(10) preliminary prior objects, preferably all with different L
    Operations->engine = 0;       // SetPrior
    Operations->iType  = -1;    Operations->i = 0;
    Operations->jType  = 0;
    Operations->kType  = 0;
    for( k = 0; k < NEXTRA; k++ )
    {
        for( i = 0; i < NEXTRA; i++ )
        {                     // try to ensure all L different but don't insist
            CALL( DoOperations(Operations, Common, Objects, Links, 1) )
            for( j = 0; j < k; j++ )
                if( Objects->Lhood == Lextra[j] )
                    break;
            if( j == k )
                break;
        }
        Lextra[k] = Objects->Lhood;
    }
// Set ENSEMBLE of prior objects
    for( i = 0; i < ENSEMBLE; i++ )
    {
        Oper = &Operations[i];
        Oper->engine =  0;               // SetPrior
        Oper->iType  = -1;    Oper->i = i;
        Oper->jType  =  0;
        Oper->kType  =  0;
    }
    CALL( DoOperations(Operations, Common, Objects, Links, ENSEMBLE) )
Exit:
    FREE(Operations)
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Control
//
// Purpose:   Rate-limited allowable cooling.
//            Aim for about Rate copy operations per object
//            (see Anneal for documentation on this).
//
//            Guard against accidental near-coalescence of likelihood values
//            (which would allow arbitrarily large cooling), by including
//            artificial extra loglikelihood values, perhaps from old ensembles
//
// History:   John Skilling   28 Jan 2002, 25 Mar 2002, 18 Jul 2002
//                            17 Sep 2002  Rate pertains to most likely object
//                            17 Dec 2002  Lextra kept different
//                            20 Aug 2003  Separate Sort2 call
//                             3 Feb 2004  Use <# copies> not #copies(Lmax)
//-----------------------------------------------------------------------------
static
int Control(        //   O  0, or -ve error
CommonStr* Common,  // I    General information
ObjectStr* Objects, // I    Sample objects, with Lhood               [ENSEMBLE]
int        Nextra,  // I    # extra objects
double*    Lextra,  // I O  Extra objects for stability                [Nextra]
double*    dcool)   //   O  cooling increment allowed by Rate
{
    double    Rate     = Common->Rate;
    int       ENSEMBLE = Common->ENSEMBLE;
    unsigned* Rand     = Common->Rand;
    double    Lmid;           // central likelihood value
    double    R;              // Rate, adjusted for finite allowable # copies
    int       N;              // total # likelihoods
    int       i, j;           // counters
    int       chop;           // binary chop counter
    double    a, copya;       // binary chop x and y (left)
    double    b, copyb;       // binary chop x and y (mid)
    double    c, copyc;       // binary chop x and y (right)
    double*   wa;             // binary chop data vector (left)
    double*   wb;             // binary chop data vector (mid)
    double*   wc;             // binary chop data vector (right)
    double*   work1  = NULL;  // workspace
    double*   work2  = NULL;  // workspace
    double*   work3  = NULL;  // workspace
    double*   Lvalue = NULL;  // workspace for sorted likelihoods
    double*   swap;           // exchange pointer
    double    Lmax;           // largest logL
    int       CALLvalue = 0;

// Default if all L equal or other anomaly
    *dcool = 0.0;

// Read all likelihoods
    N = ENSEMBLE + Nextra;
    CALLOC(Lvalue, N, double)
    CALLOC(work1, N, double)
    CALLOC(work2, N, double)
    CALLOC(work3, N, double)
    for( j = 0; j < ENSEMBLE; j++ )
        Lvalue[j] = Objects[j].Lhood;
    for( ; j < N; j++ )
        Lvalue[j] = Lextra[j-ENSEMBLE];
    if( N <= 1 )
        return 0.0;

// Offset loglikelihoods to MAX=0 for safety
    Lmax = Lvalue[0];
    for( i = 1; i < N; i++ )
        if( Lmax < Lvalue[i] )
            Lmax = Lvalue[i];
    for( i = 0; i < N; i++ )
        Lvalue[i] -= Lmax;
// Max # possible copies = max # recipients (per sample)
    R = N;
    for( i = 0; i < N; i++ )
        if( Lvalue[i] == 0.0 )
            R -= 1.0;
    R /= N;
    if( R > 0.0 )
    {
// Ensure # copies < max possible (the 3. is for rough backward compatibility)
        R = 1.0 / (3. / Rate + 1.0 / R);
// Guess initial value for dcool from normally distributed Lvalues
        a = b = 0.0;
        for( i = 0; i < N; i++ )
            a += Lvalue[i];
        a /= N;
        for( i = 0; i < N; i++ )
        {
            c = Lvalue[i] - a;
            b += c * c;
        }
        b = sqrt(b / N);
        a = b = 2.5066 * R / b;
// Bracket dcool by interval [a,b] a factor of 2 wide
        wa = work1;   wb = work2;
        for( i = 0; i < N; i++ )
            wa[i] = wb[i] = exp(a * Lvalue[i]);
        copya = copyb = Copies(wa, N);
        if( copya < R )
        {
            do
            {
                a = b;    copya = copyb;
                swap = wa;    wa = wb;    wb = swap;
                b = 2.0 * a;
                for( i = 0; i < N; i++ )
                    wb[i] = wa[i] * wa[i];
                copyb = Copies(wb, N);
            }
            while( copyb < R );
        }
        else
        {
            do
            {
                b = a;    copyb = copya;
                swap = wb;    wb = wa;    wa = swap;
                a = b / 2.0;
                for( i = 0; i < N; i++ )
                    wa[i] = sqrt(wb[i]);
                copya = Copies(wa, N);
            }
            while( copya > R );
        }
// Binary chop to accuracy 1 in 1000
        wc = work3;
        for( chop = 0; chop < 10; chop++ )
        {
            c = (a + b) / 2;
            for( i = 0; i < N; i++ )
                wc[i] = sqrt(wa[i] * wb[i]);
            copyc = Copies(wc, N);
            if( copyc < R )
            {
                a = c;    copya = copyc;
                swap = wa;    wa = wc;    wc = swap;
            }
            else
            {
                b = c;    copyb = copyc;
                swap = wb;    wb = wc;    wc = swap;
            }
        }
// Final linear interpolation
        *dcool = (b * (R - copya) + a * (copyb - R)) / (copyb - copya);
    }

    if( Nextra )
    {
// Keep Lextra up-to-date enough, trying to keep all different
        for( j = 0; j < Nextra; j++ )
        {
            Lmid = Objects[Rangrid(Rand, ENSEMBLE)].Lhood;
            for( i = 0; i < Nextra; i++ )
                if( Lextra[i] == Lmid )
                    break;
            if( i == Nextra )
            {
                Lextra[Rangrid(Rand, Nextra)] = Lmid;
                break;
            }
        }
    }

// Phenomenological correction for finite ENSEMBLE, reducing effective dcool
// (Undo correction at top of Anneal)
    *dcool *= (N - 1.0) / N;
Exit:
    FREE(work3);
    FREE(work2);
    FREE(work1);
    FREE(Lvalue);
    return  CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Copies
//
// Purpose:   <# copy operations per object> for weights w
//
//              = SUM | w - <w> |  /  2 SUM w
//
// History:   JS     3 Feb 2004
//-----------------------------------------------------------------------------
static
double Copies(
double*  w,      // I   weights
int      N)      // I   dimension
{
    double sumw, wbar, copy;
    int    i;

    sumw = 0.0;
    for( i = 0; i < N; i++ )
        sumw += w[i];
    wbar = sumw / N;
    copy = 0.0;
    for( i = 0; i < N; i++ )
        copy += fabs(w[i] - wbar);
    return  copy / (2.0 * sumw);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Anneal
//
// Purpose:   Anneal ensemble by dcool, by copying from source objects "src"
//            to destination objects "dest", as uniformly as possible
//            consistently with excess/deficit real multiplicities N:-
//
//               N(src)  = exp(dcool * (Lhood[src] - Lbar)) - 1.0
//                                                   for Lhood[src] > Lbar only
//               N(dest) = 1.0 - exp(dcool * (Lhood[dest] - Lbar))
//                                                  for Lhood[dest] < Lbar only
//            where Lbar is adjusted for equality in
//               # copy operations = SUM[src] N(src) = SUM[dest] N(dest).
//
//            A source will be copied either n or n+1 times, where
//               n = (int)N(src) = 0,1,2,...
//            A destination will be overwritten, with probability N(dest).
//
// History:   John Skilling   28 Jan 2002, 25 Mar 2002, 18 Jul 2002
//                            17 Sep 2002  Return # copies
//                             3 Feb 2004  Internal vector allocations
//-----------------------------------------------------------------------------
static
int Anneal(         //   O  #copies diagnostic, or -ve error
CommonStr* Common,  // I    CopyObject info
ObjectStr* Objects, // I    Sample objects                           [ENSEMBLE]
int        Nextra,  // I    # extra objects, for phenomenological de-correction
double     dcool)   // I    Cooling increment
{
    int       ENSEMBLE   = Common->ENSEMBLE;
    unsigned* Rand       = Common->Rand;
    OperStr*  Operations = NULL;             // list of operations   [ENSEMBLE]
    OperStr*  Oper;           // & individual operation
    double    Lmid;           // Central likelihood value
    double    r;              // Excess accumulator, also used for indexing
    double    s;              // Deficit accumulator
    double    p;              // COPY accumulator
    double    pmax;           // Accumulator limit
    double    w;              // Excess/deficit = weight - 1
    int       src;            // Identifier for duplication
    int       dest;           // Identifier for destruction
    int       i;              // ENSEMBLE counter
    int       copy;           // # copies
    double*   Lvalue = NULL;  // workspace for sorted likelihoods
    int*      Lindex = NULL;  // workspace for sorted likelihoods
    int       CALLvalue = 0;

    if( ENSEMBLE <= 1 )
        return CALLvalue;     // (nothing to do)

// Undo phenomenological correction for finite ENSEMBLE
// (Correction defined at bottom of Control)
    i = ENSEMBLE + Nextra;
    dcool /= (i - 1.0) / i;

// Read all likelihoods
    CALLOC(Lindex, ENSEMBLE, int)
    CALLOC(Lvalue, ENSEMBLE, double)
    for( i = 0; i < ENSEMBLE; i++ )
    {
        Lindex[i] = i;
        Lvalue[i] = Objects[i].Lhood;
    }
// Sort (index,likelihood) to increasing likelihood
    Sort2(Lindex, Lvalue, ENSEMBLE);
// Offset loglikelihoods to MAX=0 for safety and convenience
    for( i = 0; i < ENSEMBLE; i++ )
        Lvalue[i] -= Lvalue[ENSEMBLE - 1];

// Central likelihood value
    r = 0.0;
    for( i = 0; i < ENSEMBLE; i++ )
        r += exp(dcool * Lvalue[i]);
    Lmid = log(r / ENSEMBLE) / dcool;
// Find total excess and (equal) total deficit weights
    r = s = 0.0;
    for( i = 0; i < ENSEMBLE; i++ )
    {
        w = Lvalue[i] - Lmid;
        if( w > 0.0 )
            r += exp(dcool * w) - 1.0;
        if( w < 0.0 )
            s -= exp(dcool * w) - 1.0;
    }
    pmax = (r < s) ? r : s;   // (use the smaller for safety)
// Copy whenever cumulant excess/deficit weights step another integer
// Use ordered likelihoods as a refinement to cool more accurately
    CALLOC(Operations, ENSEMBLE, OperStr)
    copy = 0;
    r = s = 0.0;
    dest = src = -1;
    for( p = Randouble(Rand); p < pmax; p += 1.0 )
    {
        while( r < p )
        {
            w = Lvalue[++src] - Lmid;
            if( w > 0.0 )
                r += exp(dcool * w) - 1.0;
        }
        while( s < p )
        {
            w = Lvalue[++dest] - Lmid;
            if( w < 0.0 )
                s -= exp(dcool * w) - 1.0;
        }
        Oper = &Operations[copy];
        Oper->engine = 1;                         // CopyObject
        Oper->iType = -1;    Oper->i = Lindex[dest];
        Oper->jType =  1;    Oper->j = Lindex[src];
        Oper->kType =  0;
        copy++;
    }
    CALL( DoOperations(Operations, Common, Objects, NULL, copy) )
    CALLvalue = copy;
Exit:
    FREE(Operations)
    FREE(Lvalue)
    FREE(Lindex)
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Sort2
//
// Purpose:   Sort on value into increasing order
//
// History:   JS    20 Aug 2003     From Qsort of 29 Aug 1991
//-----------------------------------------------------------------------------
static
void Sort2(
int*     Lindex,    // I O  subsidiary to Lvalue
double*  Lvalue,    // I O  sort on this
int      N)         // I    dimension
{
    double    r;              // accumulator
    int       j, k, l, m, q;  // indexing counters
    l = N / 2;
    k = N - 1;
    while( 1 )
    {
        if( l > 0 )
        {
            q = Lindex[--l];    r = Lvalue[l];
        }
        else
        {
            q = Lindex[k];    r = Lvalue[k];
            Lindex[k] = Lindex[0];    Lvalue[k--] = Lvalue[0];
            if( k == 0 )
            {
                Lindex[0] = q;    Lvalue[0] = r;
                break;
            }
        }
        m = l;
        j = l + l + 1;
        while( j <= k )
        {
            if( j < k && Lvalue[j] < Lvalue[j+1] )
                j++;
            if( r < Lvalue[j] )
            {
                Lindex[m] = Lindex[j];    Lvalue[m] = Lvalue[j];
                m = j;
                j = j + j + 1;
            }
            else
                j = k + 1;
        }
        Lindex[m] = q;    Lvalue[m] = r;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FillOcean
//
// Purpose:   Width of atom = occupiable fraction of hypercube (as logarithm).
//
// Method:    Fill Ocean of previous atoms ('o' in diagram) along Hilbert line.
//            Estimate semi-width (diagram shows 2 atoms) for each current
//            atom ('|') in the Ocean, allowing a safety reduction below the
//            obvious estimate 0.5 * (# atoms in Ocean) / (# atoms in Object).
//
//            .o..oo...oo.oo.....oo.ooo....ooooooo..ooooo...ooo.o...o...o....o.
//                   <--|-->      <--|-->    <-|->   <-|->    <-----|----->
//
//            Construct symmetric window ('<-- -->' in diagram), to
//            whichever left or right displacement across Ocean is closer.
//            Then cancel the safety reduction to get "typical" width of atom.
//
//            Dense regions of Ocean have small widths,
//            sparse regions with isolated atoms have large widths.
//
// Notes: (1) Uses one standard Hilbert line so boundaries might become visible
//        (2) This procedure is for display, it's not solid Bayesian analysis.
//        (3) Not parallelisable.
//
// History:   JS        25 Jan 2003, 20 Oct 2003
//-----------------------------------------------------------------------------
static
int FillOcean(       //   O  0, or -ve error
CommonStr* Common,   // I O  general information
ObjectStr* Objects,  // I    ENSEMBLE of objects                [ENSEMBLE]
Node*      Ocean,    // I O  new ocean                                 [1]
int        NOCEAN)   // I    max # atoms in Ocean
{
static const double   Z        = (unsigned)(-1) + 1.0;   // 2^32
static const double   SAFE     = 0.3000;                 // safety
    int        ENSEMBLE = Common->ENSEMBLE;  // I    # objects of ensemble
    int        Nbits    = Common->Nbits;     // I    # bits per word
    int        Ndim     = Common->Ndim;      // I    # dimensions
    int        Valency  = Common->Valency;   // I    # fluxes
    unsigned*  Rand     = Common->Rand;      // I O  random generator
    int        Nsize    = Ndim + Valency + 1;
    ObjectStr* Object;          // object of ensemble
    int        Natoms;          // # atoms in object
    double*    Cube;            // location in hypercube
    unsigned*  Axes;            // location in hypercube, size 2^32
    Atom       atom[1];         // Hilbert length (standard orientation)
    Atom*      pCentre;         // position in Ocean
    int        semiwidth;       // atom count to windowframe either side
    Atom*      pLeft;           // left windowframe
    Atom*      pRight;          // right windowframe
    unsigned*  uleft;           // distance to left windowframe
    unsigned*  uright;          // distance to right windowframe
    unsigned*  uwide;           // distance to closer windowframe
    double     logvol;          // log(fraction of hypercube)
    int        i, j, k, n;      // counters
    int        CALLvalue = 0;

// For each atom in ensemble...
    for( k = 0; k < ENSEMBLE; k++ )
    {
        Object = &Objects[k];
        Natoms = Object->Natoms;
        semiwidth = (int)(0.5 * SAFE * (1.0 + NumAtoms(Ocean)) / Natoms) / 2;
        Axes   = Object->yLabel;
        atom->Label = Object->xLabel;
        uleft  = Object->xTry;
        uright = Object->yTry;
        for( j = 0; j < Natoms; j++ )
        {
// Erode old atoms
            if( NumAtoms(Ocean) && Rangrid(Rand, 2) )
            {
                i = Rangrid(Rand, NumAtoms(Ocean));
                pCentre = FindAtom(i, Ocean);
                CALL( DelAtom(pCentre, Ocean) )
            }
// Don't overfill Ocean
            while( NumAtoms(Ocean) >= NOCEAN - 1 )
            {
                i = Rangrid(Rand, NumAtoms(Ocean));
                pCentre = FindAtom(i, Ocean);
                CALL( DelAtom(pCentre, Ocean) )
            }
// New atom
            Cube = Object->Cubes[j];
            for( i = 0; i < Ndim; i++ )
                Axes[i] = (unsigned)(Z * Cube[i]);
            AxestoLine(atom->Label, Axes, Nbits, Ndim);
            CALL( InsAtom(atom, Ocean) )
            pCentre = FindHere(atom->Label, Ocean);  // necessarily present
// Window around it
            i = OrderNum(pCentre);
            n = i - semiwidth;
            while( n < 0 )
                n += NumAtoms(Ocean);
            pLeft = FindOrder(n, Ocean);
            n = i + semiwidth;
            while( n >= NumAtoms(Ocean) )
                n -= NumAtoms(Ocean);
            pRight = FindOrder(n, Ocean);
// Closer distance
            SubLabel(uleft, pCentre->Label, pLeft->Label, Ndim);
            SubLabel(uright, pRight->Label, pCentre->Label, Ndim);
            uwide = (CmpLabel(uright, uleft, Ndim) > 0) ? uleft : uright;
// Translate back to hypercube volume
            logvol = 0.0;
            for( i = Ndim-1; i >= 0; i-- )
            {
                if( uwide[i] )
                    logvol = log((double)uwide[i] + exp(logvol));
                logvol -= log(Z);                      // fraction of hypercube
            }
            logvol -= log(SAFE);                       // cancel safety factor
            logvol -= log(1.0 + (1.0 - SAFE) * exp(logvol)); // ensure vol < 1
// Store
            Cube[Nsize-1] = logvol;
        }
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  MCMCengines
//
// Purpose:   Apply MCMC engines for standard timestep
//
// History:   JS         2 Jan 2002 - 14 Jan 2003
//-----------------------------------------------------------------------------
static
int MCMCengines(     //   O  0, or -ve error
CommonStr* Common,   // I(O) general information
ObjectStr* Objects,  // I O  ENSEMBLE of objects
Node*      Links)    //(I O) empty on entry and exit
{
static double oldcool = 0.0;  // previous cool, initialise on first call only
static double Ncall[] = {1,1,1,1,1,1,1,1,1,1};  // # calls  weighted by cool
static double CPU[]   = {1,1,1,1,1,1,1,1,1,1};  // CPU time weighted by cool

    int       ENSEMBLE   = Common->ENSEMBLE;
    unsigned* Rand       = Common->Rand;
    double    cool       = Common->cool;
    OperStr*  Operations = NULL;      // list of operations             [total]
    OperStr*  Oper;                   // & individual operation
    OperStr   swap;                   // permutation exchange
    int       total;                  // # operations
    int       i, j, k;                // object numbers
    int       count;                  // operation counter
    int       engine;                 // engine id
    int       n3,n4,n5,n6,n7,n8,n9;   // numbers of engine calls per object
    double    t;                      // mean CPU time in LifeStory
    int       CALLvalue  = 0;
// Initialisation
    if( oldcool == 0.0 && cool > 0.0 )
    {
        t = 4.0000 * cool;
        Ncall[3]=Ncall[4]=Ncall[5]=Ncall[6]=Ncall[7]=Ncall[8]=Ncall[9] = t;
        CPU[3] = CPU[4] = CPU[5] = CPU[6] = CPU[7] = CPU[8] = CPU[9] = t;
    }

// Regenerate Links with revised Cube-to-Label mapping
    Topology(Common);

// Count MCMC operations
    n3 = n4 = n5 = n6 = n7 = n8 = n9 = 0;
    if( Common->Method & 2 && Common->MaxAtoms != 1 )
    {
        n4 = 1;    t = CPU[4] / Ncall[4];
    }
    else
    {
        n3 = 1;    t = CPU[3] / Ncall[3];
    }
    if( Common->Method & 4 && ENSEMBLE > 1 )
        n5 = (int)(t / (CPU[5] / Ncall[5]));
    if( Common->Method & 8 && ENSEMBLE > 1 )
        n6 = (int)(t / (CPU[6] / Ncall[6]));
    if( Common->Method & 16 )
        n7 = (int)(t / (CPU[7] / Ncall[7]));
    if( Common->Method & 32 )
        n8 = (int)(t / (CPU[8] / Ncall[8]));
    if( Common->Method & 64 && ENSEMBLE > 1 )
        n9 = (int)(t / (CPU[9] / Ncall[9]));
    total = (n3 + n4 + n5 + n6 + n7 + n8 + n9) * ENSEMBLE;
    CALLOC(Operations, total, OperStr)

// Set and perform MCMC entry operations
    for( count = 0; count < ENSEMBLE; count++ )
    {
        Oper = &Operations[count];
        Oper->engine = 2;                          // MCMCsetup
        Oper->iType = -1;    Oper->i = count;
        Oper->jType =  0;
        Oper->kType =  0;
    }
    CALL( DoOperations(Operations, Common, Objects, Links, ENSEMBLE) )

// Set MCMC operations
    Oper = Operations;
    for( count = 0; count < n3; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 3;                      // LifeStory1
            Oper->iType = -1;    Oper->i = i;
            Oper->jType =  0;
            Oper->kType =  0;
            Oper++;
        }
    for( count = 0; count < n4; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 4;                      // LifeStory2
            Oper->iType = -1;    Oper->i = i;
            Oper->jType =  0;
            Oper->kType =  0;
            Oper++;
        }
    for( count = 0; count < n5; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 5;                      // Chameleon1
            do j = Rangrid(Rand, ENSEMBLE);
            while( j == i );                       // objects must be different
            Oper->iType = -1;    Oper->i = i;
            Oper->jType = -1;    Oper->j = j;
            Oper->kType =  0;
            Oper++;
        }
    for( count = 0; count < n6; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 6;                      // Chameleon2
            do j = Rangrid(Rand, ENSEMBLE);
            while( j == i );                       // objects must be different
            Oper->iType = -1;    Oper->i = i;
            Oper->jType = -1;    Oper->j = j;
            Oper->kType =  0;
            Oper++;
        }
    for( count = 0; count < n7; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 7;                      // Leapfrog1
            j = Rangrid(Rand, ENSEMBLE);           // objects can be same
            Oper->iType = -1;                  Oper->i = i;
            Oper->jType = (j == i) ? 0 : 1;    Oper->j = j;
            Oper->kType = 0;
            Oper++;
        }
    for( count = 0; count < n8; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 8;                      // Leapfrog2
            j = Rangrid(Rand, ENSEMBLE);
            if( j == i )
                k = i;
            else       // either both or neither neighbour from original object
                do  k = Rangrid(Rand, ENSEMBLE);
                while( k == i );
            Oper->iType = -1;                  Oper->i = i;
            Oper->jType = (j == i) ? 0 : 1;    Oper->j = j;
            Oper->kType = (k == j) ? 0 : 1;    Oper->k = k;
            Oper++;
        }
    for( count = 0; count < n9; count++ )
        for( i = 0; i < ENSEMBLE; i++ )
        {
            Oper->engine = 9;                      // GuidedWalk
            do  j = Rangrid(Rand, ENSEMBLE);
            while( j == i );
            do  k = Rangrid(Rand, ENSEMBLE);
            while( k == i );
            Oper->iType = -1;                  Oper->i = i;
            Oper->jType =  1;                  Oper->j = j;
            Oper->kType = (k == j) ? 0 : 1;    Oper->k = k;
            Oper++;
        }

// Perform MCMC operations
// Individual engines have detailed balance, but only LifeStory explores
//  fully, so LifeStory alone ought to equilibrate (eventually).
// Other engines are slaved to LifeStory's CPU time: the slaving ratio
//  is locked after burn-in, so evolution has exact detailed balance.
// Permute operations for visible detailed balance.
    for( count = 0; count < total; count++ )
    {
        i = Rangrid(Rand, total - count);
        if( i )
            SWAP(Operations[count], Operations[count + i])
    }
// Perform in permuted order
    CALL( DoOperations(Operations, Common, Objects, Links, total) )

// Statistics for tuning
    for( count = 0; count < total; count++ )
    {
        Oper = &Operations[count];
        Common->Success += Oper->Success;
        Common->CPU     += Oper->CPU;
        if( cool > oldcool )
        {
            engine = (int)Oper->engine;
            CPU  [engine] = CPU  [engine] + cool * Oper->CPU;
            Ncall[engine] = Ncall[engine] + cool;
        }
    }

// Set and perform MCMC exit operations
    for( count = 0; count < ENSEMBLE; count++ )
    {
        Oper = &Operations[count];
        Oper->engine = 10;                         // MCMCexit
        Oper->iType = -1;    Oper->i = count;
        Oper->jType =  0;
        Oper->kType =  0;
    }
    CALL( DoOperations(Operations, Common, Objects, Links, ENSEMBLE) )
    oldcool = cool;
Exit:
    FREE(Operations)
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  BayesFree
//
// Purpose:   Free allocated memory for BayeSys (and MassInf)
//
// History:   JS        16 Oct 2002, 16 Nov 2002, 1 Feb 2003
//-----------------------------------------------------------------------------
static
void BayesFree (
CommonStr* Common,   // I(O) general information
ObjectStr* Objects,  //  (O) old objects
Node*      Links,    //  (O) old trees
Node*      Ocean)    //  (O) old ocean
{
    int        ENSEMBLE  = Common->ENSEMBLE;
    int        Valency   = Common->Valency;
    ObjectStr* Object;        // local object
    int        j;             // atom counter
    int        k;             // object counter

    FreeLink(Ocean);
    if( Valency )
        FluxFree(Common, Objects);
    for( k = 0; k < ENSEMBLE; k++ )
    {
        Object = &Objects[k];
        if( Object->Cubes )
            for( j = 0; j < Object->Nstore; j++ )
                FREE(Object->Cubes[j])
        FREE(Object->Cubes)
        FREE(Object->work)
        FREE(Object->yTry)
        FREE(Object->xTry)
        FREE(Object->yOrigin)
        FREE(Object->xOrigin)
        FREE(Object->yLabel)
        FREE(Object->xLabel)
        FreeLink(&Links[k]);
    }
    FREE(Common->permute)
    FREE(Common->offset)
}

//=============================================================================
//
//                               Engines
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  DoOperations
//
// Purpose:   Perform list of operations with parallelisable engine calls
//
// History:   JS        25 Oct 2002
//-----------------------------------------------------------------------------
#if ! PARALLEL
static
int DoOperations(        //   O  0, or -ve error
OperStr*   Operations,   // I O  all operations
CommonStr* Common,       // I O  general information
ObjectStr* Objects,      // I O  ENSEMBLE of objects
Node*      Links,        // I O  linked lists of labels
int        nOper)        // I    # operations
{
    int       op;                  // operation counter
    int       CALLvalue = 0;

    for( op = 0; op < nOper; op++ )
    {
        Operations[op].slave = 0;
        CALL( Do1operation(Operations, Common, Objects, Links, op) )
    }
Exit:
    return CALLvalue;
}
#endif

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Do1operation
//
// Purpose:   Distribution procedure for parallelisable engine calls
//
// History:   JS        25 Oct 2002, 1 Feb 2003
//-----------------------------------------------------------------------------
static
int Do1operation(        //   O  0, or -ve error
OperStr*   Operations,   // I O  all operations
CommonStr* Common,       // I O  general information
ObjectStr* Objects,      // I O  ENSEMBLE of objects
Node*      Links,        // I O  linked lists of labels
int        op)           // I    individual operation
{
    OperStr*   Oper    = &Operations[op];
    int        CALLvalue = 0;

    switch( Oper->engine )
    {
        case 0:    CALL( SetPrior  (Oper, Common, Objects, Links) )    break;
        case 1:    CALL( CopyObject(Oper, Common, Objects       ) )    break;
        case 2:    CALL( MCMCsetup (Oper, Common, Objects, Links) )    break;
        case 3:    CALL( LifeStory1(Oper, Common, Objects, Links) )    break;
        case 4:    CALL( LifeStory2(Oper, Common, Objects, Links) )    break;
        case 5:    CALL( Chameleon1(Oper, Common, Objects, Links) )    break;
        case 6:    CALL( Chameleon2(Oper, Common, Objects, Links) )    break;
        case 7:    CALL( Leapfrog1 (Oper, Common, Objects, Links) )    break;
        case 8:    CALL( Leapfrog2 (Oper, Common, Objects, Links) )    break;
        case 9:    CALL( GuidedWalk(Oper, Common, Objects, Links) )    break;
        case 10:   CALL( MCMCexit  (Oper, Common, Objects, Links) )    break;
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SetPrior                                       Parallel engine #0
//
// Purpose:   Sample the prior to set an Object of ENSEMBLE.
//            The prior on the number of atoms is as follows.
//
//  (1) Alpha = 0.
//      The prior is UNIFORM (and MAX must be finite).
//
//                Pr(N) = 1 / (MAX+1-MIN)
//
//                <N> = (MAX+MIN) / 2,     var(N) = (MAX-MIN) (MAX-MIN+2) / 12
//
//  (2) Alpha > 0.
//      The prior is BINOMIAL.
//
//                   (MAX-MIN)!      N-MIN      MAX-N              Alpha
//      Pr(N)  =  ----------------  q      (1-q)      ,   q =  -------------
//                (N-MIN)!(MAX-N)!                             Alpha+MAX-MIN
//
//      <N> = (1-q) MIN + q MAX  ,      var(N) = (MAX-MIN) q (1-q)
//
//
//      Special case  MAX = infinity  becomes POISSON in (N-MIN).
//                               -alpha      N-MIN
//                    Pr(N)  =  e       alpha      / (N-MIN)!   ,   N >= MIN
//
//                    <N> = MIN + alpha  ,      var(N) = alpha
//
//
//  (3) Alpha < 0.
//      The prior is truncated GEOMETRIC.
//
//                  N-MIN            MAX-MIN+1
//      Pr(N)  =   a     (1-a) / (1-a         )  ,   a = |Alpha| / (|Alpha|+1)
//
//      The mean <N> and variance var(N) have explicit but uninformative forms.
//
//
//      Special case  MAX = infinity  becomes GEOMETRIC with ratio  a .
//                                     N-MIN
//                    Pr(N)  =  (1-a) a       ,  N >= MIN
//
//                    <N> = MIN + alpha  ,      var(N) = alpha (alpha + 1)
//
//
// History:   JS         24 Jan 2002, 30 Sep 2002, 16 Oct 2002, 14 Jan 2003
//                        1 Feb 2003, 8 Feb 2002
//-----------------------------------------------------------------------------
static
int SetPrior(          //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I    general information
ObjectStr* Objects,    //   O  sample objects
Node*      Links)      //(I O) empty on entry and exit
{
    Node*      Link     = &Links  [Oper->i];
    ObjectStr* Object   = &Objects[Oper->i];
    unsigned*  Label    = Object->xLabel;
    unsigned*  Rand     = Object->Rand;
    int        MinAtoms = Common->MinAtoms;
    int        MaxAtoms = Common->MaxAtoms;
    int        Ndim     = Common->Ndim;
    double     Alpha    = Common->Alpha;
    double     Ltry;
    int        j;
    int        Natoms;
    Atom*      pAtom;
    int        CALLvalue = 0;

// # atoms
    Natoms = MinAtoms;
    if( MaxAtoms > MinAtoms )
    {
        j = MaxAtoms - MinAtoms;
        if( Alpha == 0.0 )
            Natoms += Rangrid(Rand, 1 + j);
        else if( Alpha > 0.0 )
            Natoms += Ranbinom(Rand, j, Alpha / (Alpha + j));
        else  // Alpha < 0.0
            Natoms += Rangeom(Rand, -Alpha, 1 + j);
    }
    else if( MaxAtoms == 0 )
    do
    {
        Natoms = MinAtoms;
        if( Alpha > 0.0 )
            Natoms += Ranpoiss(Rand, Alpha);
        else // Alpha < 0.0 (Alpha=0 excluded)
            Natoms += (int)(log(Randouble(Rand)) / log(Alpha / (Alpha - 1.0)));
    } while( Natoms < MinAtoms );         // technical overflow protection
    Object->Natoms = Natoms;
    CALL( SafeCube(Common, Object) )
// Sample the Cube and avoid spatial coincidences by filling Link
    do
    {
        Object->Natoms = 0;               // no Cubes exist
        CALL( Empty(Common, Object) )
        for( j = 0; j < Natoms; j++ )
        {
            do
            {
                do
                {
                    RanLabel(Label, Ndim, Rand);
                    LabeltoCube(Object, Label, Common);
                } while( FindHere(Label, Link) );
                CALL( Try1(&Ltry, Common, Object) )
            } while( CALLvalue == 0 );    // repeat until atom accepted by user
            CALL( Insert1(Common, Object, Link) )
        }
        for( pAtom = BeginLink(Link); pAtom; pAtom = BeginLink(Link) )
            CALL( DelAtom(pAtom, Link) )  // empty Link for later use
    } while( j < Natoms );
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  CopyObject                                     Parallel engine #1
//
// Purpose:   Copy source ensemble object to destination, "*dest = *src".
//            Ideally, destination object should be rebuilt afterwards to
//            recover any nuisance parameters.
//
// History:   JS          25 Apr 2001, 2 Jan 2002, 21 Mar 2002, 30 Sep 2002,
//                         1 Feb 2003, 8 Feb 2003, 20 Aug 2003
//-----------------------------------------------------------------------------
static
int CopyObject(        //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) General information
ObjectStr* Objects)    // I O  Ensemble of objects
{
    int        Ndim      = Common->Ndim;
    int        Valency   = Common->Valency;
    int        Nsize     = Ndim + Valency + 1;
    int        dest      = Oper->i;
    int        src       = Oper->j;
    ObjectStr* yObject   = &Objects[dest];
    ObjectStr* xObject   = &Objects[src];
    int        Natoms    = xObject->Natoms;
    int        i, j;
    int        CALLvalue = 0;

    if( dest != src )
    {
        yObject->Lhood = xObject->Lhood;
        yObject->Natoms = Natoms;
        CALL( SafeCube(Common, yObject) )
        for( j = 0; j < Natoms; j++ )
            for( i = 0; i < Nsize; i++ )
                yObject->Cubes[j][i] = xObject->Cubes[j][i];
        if( Valency )
            yObject->FluxUnit = xObject->FluxUnit;
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  MCMCsetup                                      Parallel engine #2
//
// Purpose:   Re-generate Links in new topology
//
// History:   JS        16 Oct 2002, 8 Feb 2003, 20 Aug 2003
//-----------------------------------------------------------------------------
static
int MCMCsetup(         //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I    general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      //   O  new trees
{
    int        MinAtoms = Common->MinAtoms;
    int        MaxAtoms = Common->MaxAtoms;
    ObjectStr* Object   = &Objects[Oper->i];
    Node*      Link     = &Links  [Oper->i];
    int        Natoms   = Object->Natoms;
    int        j;
    int        CALLvalue = 0;

    Object->Natoms = 0;
    CALL( Empty(Common, Object) )
    Object->reset = 0;                            // keep MassInf fluxes
    for( j = 0; j < Natoms; j++ )
        CALL( Insert1(Common, Object, Link) )
    if( Object->Natoms < MinAtoms )
        return E_BAYESYS_SYSERR;
    if( MaxAtoms && Object->Natoms > MaxAtoms )
        return E_BAYESYS_SYSERR;
    Object->reset = 1;                            // recover default state
    if( Common->Valency )
        FluxCalib(Common, Object);                // Re-calibrate object
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  LifeStory1                                     Parallel engine #3
//
// Purpose:   MCMC engine, cutdown alternate to 2-atom version LifeStory2.
//
//            Birth/Move/Death for one atom in an object.  Explores all space.
//            Success defined as birth or death if this is allowed,
//            else as movement.
//
// History:   JS         2 Jan 2002 - 8 Feb 2003
//-----------------------------------------------------------------------------
static
int LifeStory1(        //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      // I O  linked lists of labels
{
// Sampling interval: each atom is investigated about TIMESTEP times
static const double   TIMESTEP = 2.0000;
    int        Ndim     = Common->Ndim;
    double     cool     = Common->cool;
    int        Nbits    = Common->Nbits;
    int        MinAtoms = Common->MinAtoms;
    int        MaxAtoms = Common->MaxAtoms;
    ObjectStr* Object   = &Objects[Oper->i];
    Node*      Link     = &Links  [Oper->i];
    unsigned*  xLabel   = Object->xLabel;
    unsigned*  Origin   = Object->xOrigin;
    unsigned*  xTry     = Object->xTry;
    unsigned*  Rand     = Object->Rand;
    int        n;
    int        del;
    int        jbits;
    double     accept;
    double     Ltry;
    double     L0;
    double     L1;
    double     b;
    double     d;
    double     Clock;
    Atom*      pLeft;
    Atom*      pRight;
    unsigned*  uLeft;
    unsigned*  uRight;
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;

    b = Birth(Common, Object->Natoms);
    d = Death(Common, Object->Natoms);
    for( Clock  = -log(Randouble(Rand)) / (b + d) ;
         Clock  < TIMESTEP ;
         Clock += -log(Randouble(Rand)) / (b + d) )
    {
        n = Object->Natoms;
        if( (b + d) * Randouble(Rand) > b )
        {
// Death: Natoms-1 or Natoms
            del = Rangrid(Rand, Object->Natoms);
            CopyLabel(xLabel, FindAtom(del,Link)->Label, Ndim);
            CALL( Delete1(del, &L1, NULL, Common, Object, Link) )     CPU += 1;
        }
        else
        {
// Birth: Natoms or Natoms+1
            CALL( SafeCube(Common, Object) )
            RanLabel(xLabel, Ndim, Rand);
            if( FindHere(xLabel, Link) )
                continue;
            LabeltoCube(Object, xLabel, Common);
            CALL( Try1(&L1, Common, Object) )                         CPU += 1;
            if( CALLvalue == 0 )          // abort if location rejected by user
                continue;
        }
        L0 = Object->Lhood;
// Neighbours of selected position
        if( !(pLeft  = FindLeft (xLabel, Link)) )   pLeft  = EndLink(Link);
        if( !(pRight = FindRight(xLabel, Link)) )   pRight = BeginLink(Link);
        uLeft  = pLeft  ? pLeft->Label  : NULL;
        uRight = pRight ? pRight->Label : NULL;
// Prepare to move
        Ltry = (MaxAtoms == MinAtoms) ? cool*L1 : PLUS(cool*L0, cool*L1);
        accept = Ltry + log(Randouble(Rand));
        RanLabel(Origin, Ndim, Rand);
        for( jbits = Nbits * Ndim; jbits >= 0; jbits -= 1 )
        {
// Slice sampling for optional atom
            NewLabel(xTry, xLabel, Origin, Ndim, jbits, Nbits, Rand);
            if( Object->Natoms && Outside(xTry, uLeft, uRight, Ndim) )
                continue;
            LabeltoCube(Object, xTry, Common);
            CALL( Try1(&L1, Common, Object) )                         CPU += 1;
            if( CALLvalue == 0 )          // abort if location rejected by user
                continue;
            Ltry = (MaxAtoms == MinAtoms) ? cool*L1 : PLUS(cool*L0, cool*L1);
            if( Ltry >= accept )
                break;
        }
// Birth/Death/Move
        if( log(Randouble(Rand)) < cool * L1 - Ltry )
        {
            CALL( Insert1(Common, Object, Link) )                     CPU += 1;
        }

// Update event rates
        b = Birth(Common, Object->Natoms);
        d = Death(Common, Object->Natoms);
        if( (MaxAtoms != MinAtoms && Object->Natoms != n)
         || (MaxAtoms == MinAtoms && jbits > 0) )                   ++ Success;
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  LifeStory2                                     Parallel engine #4
//
// Purpose:   LifeStory MCMC diffusion engine
//
//            Birth/Move/Death for two adjacent atoms.  Explores all space.
//            Success defined as birth or death if this is allowed,
//            else as movement.
//
// History:   JS         2 Jan 2002 - 8 Feb 2003
//-----------------------------------------------------------------------------
static
int LifeStory2(        //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      // I O  linked lists of labels
{
// Sampling interval: each atom is investigated about TIMESTEP times
static const double   TIMESTEP = 2.0000;
    int        Ndim     = Common->Ndim;
    double     cool     = Common->cool;
    int        Nbits    = Common->Nbits;
    int        MinAtoms = Common->MinAtoms;
    int        MaxAtoms = Common->MaxAtoms;
    ObjectStr* Object   = &Objects[Oper->i];
    Node*      Link     = &Links  [Oper->i];
    unsigned*  xLabel   = Object->xLabel;
    unsigned*  yLabel   = Object->yLabel;
    unsigned*  xOrigin  = Object->xOrigin;
    unsigned*  yOrigin  = Object->yOrigin;
    unsigned*  xTry     = Object->xTry;
    unsigned*  yTry     = Object->yTry;
    unsigned*  Rand     = Object->Rand;
    int        del;
    int        n;
    int        jbits;
    double     accept;
    double     Ltry;
    double     L1;
    double     L2;
    double     b;
    double     d;
    double     Clock;
    Atom*      pAtom;
    Atom*      pLeft;
    Atom*      pRight;
    unsigned*  uLeft;
    unsigned*  uRight;
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;

    b = Birth(Common, Object->Natoms);
    d = Death(Common, Object->Natoms);
    for( Clock  = -log(Randouble(Rand)) / (b + d) ;
         Clock  < TIMESTEP ;
         Clock += -log(Randouble(Rand)) / (b + d) )
    {
        n = Object->Natoms;
        if( (b + d) * Randouble(Rand) > b )
        {
// Death: Natoms-1 or Natoms
            L2 = Object->Lhood;
            del = Rangrid(Rand, Object->Natoms);
            CopyLabel(xLabel, FindAtom(del,Link)->Label, Ndim);
            CALL( Delete1(del, NULL, NULL, Common, Object, Link) )    CPU += 1;
        }
        else
        {
// Birth: Natoms or Natoms+1
            RanLabel(xLabel, Ndim, Rand);
            if( FindHere(xLabel, Link) )
                continue;
            LabeltoCube(Object, xLabel, Common);
            CALL( Try1(&L2, Common, Object) )                         CPU += 1;
            if( CALLvalue == 0 ) //must abort here if location rejected by user
                continue;
            CALL( SafeCube(Common, Object) )
        }
// Kill stimulation neighbour
        if( Rangrid(Rand, 2) )
        {
            if( !(pAtom = FindLeft(xLabel,Link)) )  pAtom = EndLink(Link);
        }
        else
        {
            if( !(pAtom = FindRight(xLabel,Link)) )  pAtom = BeginLink(Link);
        }
        del = Storage(pAtom);
        CALL( Delete1(del, &L1, &L2, Common, Object, Link) )          CPU += 1;
        Ltry = (MaxAtoms == MinAtoms) ? cool*L2 : PLUS(cool*L1, cool*L2);
// Neighbours of selected pair
        CubetoLabel(yLabel, Object, Common);
        if( !(pLeft  = FindLeft (yLabel, Link)) )   pLeft  = EndLink(Link);
        if( !(pRight = FindRight(yLabel, Link)) )   pRight = BeginLink(Link);
        uLeft  = pLeft  ? pLeft->Label  : NULL;
        uRight = pRight ? pRight->Label : NULL;
// Prepare to move
        accept = Ltry + log(Randouble(Rand));
        RanLabel(xOrigin, Ndim, Rand);
        RanLabel(yOrigin, Ndim, Rand);
        for( jbits = Nbits * Ndim; jbits >= 0; jbits -= 1 )
        {
// Slice sampling for replacement atom and optional atom
            NewLabel(xTry, xLabel, xOrigin, Ndim, jbits, Nbits, Rand);
            if( Object->Natoms && Outside(xTry, uLeft, uRight, Ndim) )
                continue;
            NewLabel(yTry, yLabel, yOrigin, Ndim, jbits, Nbits, Rand);
            if( CmpLabel(xTry, yTry, Ndim) == 0 )
                continue;
            if( Object->Natoms && Outside(yTry, uLeft, uRight, Ndim) )
                continue;
            Object->Natoms++;
            LabeltoCube(Object, xTry, Common);
            Object->Natoms--;
            LabeltoCube(Object, yTry, Common);
            CALL( Try2(&L1, &L2, Common, Object) )                    CPU += 2;
            if( CALLvalue == 0 )          // abort if location rejected by user
                continue;
            Ltry = (MaxAtoms == MinAtoms) ? cool*L2 : PLUS(cool*L1, cool*L2);
            if( Ltry >= accept )
                break;
        }
// Birth/Death/Move
        if( log(Randouble(Rand)) < cool * L2 - Ltry )
        {
            CALL( Insert2(Common, Object, Link) )                     CPU += 2;
        }
        else
        {
            CALL( Insert1(Common, Object, Link) )                     CPU += 1;
        }
// Update event rates
        b = Birth(Common, Object->Natoms);
        d = Death(Common, Object->Natoms);
        if( (MaxAtoms != MinAtoms && Object->Natoms != n)
         || (MaxAtoms == MinAtoms && jbits > 0) )                   ++ Success;
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Chameleon1                                     Parallel engine #5
//
// Purpose:   MCMC engine.
//
//            Tries to move atoms from one object to another,
//            WITHOUT geometrical exploration.  Success defined as acceptance.
//
//                 ========X===============X==================X=========
//                                         |
//                                         |
//                                        \|/
//                 ====X==============X====:=====X========X=============
//
// Note:      ENSEMBLE should be at least 2, obviously.
//
// History:   JS         31 Jan 2002 - 8 Feb 2003
//-----------------------------------------------------------------------------
static
int Chameleon1(        //   O  0 = op, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      // I O  linked lists of labels
{
// Sampling interval: each atom is investigated about TIMESTEP times
static const double   TIMESTEP = 1.0000;
    double     cool     = Common->cool;
    ObjectStr* iObject  = &Objects[Oper->i];
    ObjectStr* jObject  = &Objects[Oper->j];
    ObjectStr* kObject;
    Node*      iLink    = &Links  [Oper->i];
    Node*      jLink    = &Links  [Oper->j];
    Node*      kLink;
    unsigned*  Rand     = iObject->Rand;
    int        del;
    double     accept;
    double     L1;
    double     Ltry;
    double*    swap;
    double     b;
    double     d;
    double     constant;
    double     Clock;
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;

    b = Birth(Common, iObject->Natoms) * Death(Common, jObject->Natoms);
    d = Death(Common, iObject->Natoms) * Birth(Common, jObject->Natoms);
    if( b + d <= 0.0 )
        goto Exit;
    constant = (iObject->Natoms + jObject->Natoms) / 2.0;
    for( Clock  = -log(Randouble(Rand)) * constant / (b + d);
         Clock <= TIMESTEP;
         Clock += -log(Randouble(Rand)) * constant / (b + d) )
    {
        if( Randouble(Rand) > b / (b + d) )
        {
            kObject = iObject;   iObject = jObject;   jObject = kObject;
            kLink   = iLink;     iLink   = jLink;     jLink   = kLink;
        }
// Try atom j ---> i
        CALL( SafeCube(Common, iObject) )
// Ensure deletion not already present in destination
        del = Rangrid(Rand, jObject->Natoms);
        if( FindHere(FindAtom(del, jLink)->Label, iLink) )
            continue;
// Delete from source
        CALL( Delete1(del, &L1, NULL, Common, jObject, jLink) )         ++ CPU;
        accept = cool * (L1 + iObject->Lhood) + log(Randouble(Rand));
// Try movement
        SWAP(iObject->Cubes[iObject->Natoms],
             jObject->Cubes[jObject->Natoms])
        CALL( Try1(&Ltry, Common, iObject) )                            ++ CPU;
        Ltry = cool * (Ltry + jObject->Lhood);
// Reject?
        if( CALLvalue == 0 || Ltry < accept )
        {
            SWAP(iObject->Cubes[iObject->Natoms],
                 jObject->Cubes[jObject->Natoms])
            CALL( Insert1(Common, jObject, jLink) )                     ++ CPU;
        }
// Accept?
        else
        {                                                           ++ Success;
            CALL( Insert1(Common, iObject, iLink) )                     ++ CPU;
        }
        b = Birth(Common, iObject->Natoms) * Death(Common, jObject->Natoms);
        d = Death(Common, iObject->Natoms) * Birth(Common, jObject->Natoms);
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Chameleon2                                     Parallel engine #6
//
// Purpose:   MCMC engine.
//
// Method:    Tries to exchange positions of nearby atoms in different objects,
//            but does NOT explore all space.  Success defined as acceptance.
//
//                 =====X=====:======X============X=======
//                           /|\                Select
//                            |                   |
//                            |                   |
//                         Exchange              \|/
//                 ===========X===========X=======:=======
//
//            For reasonable success rate, the exchanging atom is taken to be
//            the left or right neighbour of the original atom, and hence
//            "close".
//
// Note:      ENSEMBLE should be at least 2, obviously.
//
// History:   JS         31 Jan 2002 - 8 Feb 2003
//-----------------------------------------------------------------------------
static
int Chameleon2(        //   O  0 = op, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node *     Links)      // I O  linked lists of labels
{
// Sampling interval: each pair is investigated about TIMESTEP times
static const double   TIMESTEP = 0.5000;
    double     cool     = Common->cool;
    ObjectStr* iObject  = &Objects[Oper->i];
    ObjectStr* jObject  = &Objects[Oper->j];
    ObjectStr* kObject;
    Node*      iLink    = &Links  [Oper->i];
    Node*      jLink    = &Links  [Oper->j];
    Node*      kLink;
    unsigned*  Rand     = iObject->Rand;
    int        iNatoms  = iObject->Natoms;
    int        jNatoms  = jObject->Natoms;
    int        iatom;    // first object
    Atom*      iAtom;    // first object
    int        jatom;    // other object
    Atom*      jAtom;    // other object
    double     accept;
    double     Ltry;
    double     L2;
    double*    swap;
    int        count;
    int        r;
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;

    for( count = 0; count <= (iNatoms + jNatoms) * TIMESTEP; count++ )
    {
// Select random atom
        r = Rangrid(Rand, iNatoms + jNatoms);
        iatom = r;
        if( r >= iNatoms )
        {
            iatom = r - iNatoms;
            kObject = iObject;   iObject = jObject;   jObject = kObject;
            kLink   = iLink;     iLink   = jLink;     jLink   = kLink;
            r       = iNatoms;   iNatoms = jNatoms;   jNatoms = r;
        }
        iAtom = FindAtom(iatom, iLink);
        if( FindHere(iAtom->Label, jLink) )  // must not be in other object
            continue;
        if( Rangrid(Rand, 2) )
        {
            if( !(jAtom = FindLeft(iAtom->Label, jLink)) )
                jAtom = EndLink(jLink);
        }
        else
        {
            if( !(jAtom = FindRight(iAtom->Label, jLink)) )
                jAtom = BeginLink(jLink);
        }
        if( FindHere(jAtom->Label, iLink) )  // must not be in other object
            continue;
        jatom = Storage(jAtom);
// Exchange coordinates
        CALL( Delete1(iatom, &Ltry, NULL, Common, iObject, iLink) )     ++ CPU;
        CALL( Delete1(jatom, &L2,   NULL, Common, jObject, jLink) )     ++ CPU;
        accept = cool * (Ltry + L2) + log(Randouble(Rand));
        SWAP(iObject->Cubes[iObject->Natoms],
             jObject->Cubes[jObject->Natoms])
        CALL( Try1(&Ltry, Common, iObject) )                            ++ CPU;
        if( CALLvalue > 0 )
        {                              // care with user rejections
            CALL( Try1(&L2, Common, jObject) )                          ++ CPU;
            Ltry += L2;
        }
// Reject?
        if( CALLvalue == 0 || cool * Ltry < accept )
           SWAP(iObject->Cubes[iObject->Natoms],
                jObject->Cubes[jObject->Natoms])
        else                                                        ++ Success;
// Rebuild
        CALL( Insert1(Common, iObject, iLink) )                         ++ CPU;
        CALL( Insert1(Common, jObject, jLink) )                         ++ CPU;
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Leapfrog1                                      Parallel engine #7
//
// Purpose:   Geometrical MCMC engine.
//
// Method:    Tries to move points by inverting w.r.t. a catalyst atom,
//            but does NOT explore all space.  Success defined as acceptance.
//
//              old position        catalyst         new position
//                 x1 >-----------------0-----------------> x2
//
//                            old + new = 2 * catalyst
//
//            For reasonable success rate, the catalyst is taken to be either
//            the left or right neighbour of the original atom, and hence
//            "close".  For detailed balance, the new position must also have
//            its right or left neighbour (respectively) being the centre,
//            otherwise the trial is aborted.  Because a Hilbert curve may find
//            similarly close atoms in any quadrant, the efficiency is roughly
//                   1  in  MIN(2^Ndim, Natoms)  tries .
//
// Note:      ENSEMBLE should be at least 3, otherwise
//            Leapfrog1 will not of itself expand the bubble of atoms.
//
// History:   JS         2 Jan 2002 - 8 Feb 2003
//-----------------------------------------------------------------------------
static
int Leapfrog1(         //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      // I O  linked lists of labels
{
// Sampling interval: each atom is investigated about TIMESTEP times
static const double   TIMESTEP = 1.0000;
    int        Ndim     = Common->Ndim;
    double     cool     = Common->cool;
    ObjectStr* iObject  = &Objects[Oper->i];
    ObjectStr* jObject  = &Objects[Oper->j];
    Node*      iLink    = &Links  [Oper->i];
    Node*      jLink    = &Links  [Oper->j];
    unsigned*  xTry     = iObject->xTry;
    unsigned*  Rand     = iObject->Rand;
    int        Natoms   = iObject->Natoms;
    int        del;      // selection
    Atom*      pAtom;    // selection
    double*    cube;     // selection
    int        left;     // catalyst
    Atom*      pCentre;  // catalyst
    double*    centre;   // catalyst
    double*    trial;
    double     Ltry;
    double     accept;
    double     t;
    double*    swap;
    int        r;
    int        count;
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;

    for( count = 0; count <= Natoms * TIMESTEP; count++ )
    {
// Select original atom
        del = Rangrid(Rand, Natoms);
        cube = iObject->Cubes[del];
        pAtom = FindAtom(del, iLink);
// Select centre
        left = Rangrid(Rand, 2);
        if( left )    // left-definite
        {
            pCentre = FindRight(pAtom->Label, jLink);
            if( ! pCentre )    pCentre = BeginLink(jLink);
            pCentre = PrevAtom(pCentre);
            if( ! pCentre )    pCentre = EndLink(jLink);
        }
        else          // right-definite
        {
            pCentre = FindLeft(pAtom->Label, jLink);
            if( ! pCentre )    pCentre = EndLink(jLink);
            pCentre = NextAtom(pCentre);
            if( ! pCentre )    pCentre = BeginLink(jLink);
        }
        centre = jObject->Cubes[Storage(pCentre)];

// Invert original cube through neighbour
        trial = iObject->Cubes[Natoms];
        for( r = 0; r < Ndim; r++ )
        {
            t = cube[r] - centre[r];
            if( t < 0.0 )
                t += 1.0;
            t = centre[r] - t;
            if( t < 0.0 )
                t += 1.0;
            trial[r] = t;
        }

// Ensure not already present
        CubetoLabel(xTry, iObject, Common);
        if( FindHere(xTry, iLink) )
            continue;
// Ensure in range for detailed balance
        if( left )    // right-definite
        {
            pAtom = FindLeft(xTry, jLink);
            if( ! pAtom )    pAtom = EndLink(jLink);
            pAtom = NextAtom(pAtom);
            if( ! pAtom )    pAtom = BeginLink(jLink);
        }
        else          // left-definite
        {
            pAtom = FindRight(xTry, jLink);
            if( ! pAtom )    pAtom = BeginLink(jLink);
            pAtom = PrevAtom(pAtom);
            if( ! pAtom )    pAtom = EndLink(jLink);
        }
        if( pAtom != pCentre )
            continue;

// Delete from source
        CALL( Delete1(del, &Ltry, NULL, Common, iObject, iLink) )       ++ CPU;
        accept = cool * Ltry + log(Randouble(Rand));
        SWAP(iObject->Cubes[Natoms], iObject->Cubes[Natoms-1])
        CALL( Try1(&Ltry, Common, iObject) )                            ++ CPU;
// Reject?
        if( CALLvalue == 0 || cool * Ltry < accept )
            SWAP(iObject->Cubes[Natoms], iObject->Cubes[Natoms-1])
        else                                                        ++ Success;

// Rebuild
        CALL( Insert1(Common, iObject, iLink) )                         ++ CPU;
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Leapfrog2                                      Parallel engine #8
//
// Purpose:   Geometrical MCMC engine.
//
// Method:    Tries to move points by reflecting between catalyst atoms,
//            but does NOT explore all space.  Success defined as acceptance.
//
//                 Lcatalyst            old position
//                                       /
//                                      /
//                                     /
//                          new position           Rcatalyst
//
//                   old + new = Lcatalyst + Rcatalyst
//
//            For reasonable success rate, the catalysts are taken to be
//            the left and right neighbours of the original atom, and hence
//            "close".  For detailed balance, the new position must also have
//            its left and right neighbours (respectively) being the catalysts,
//            otherwise the trial is aborted.  Because a Hilbert curve may find
//            similarly close atoms in any quadrant, the efficiency is roughly
//                   1  in  MIN(2^Ndim, Natoms)  tries .
//
// Note:      ENSEMBLE should be at least 4, otherwise
//            Leapfrog2 will not of itself expand the bubble of atoms.
//
// History:   JS         2 Jan 2002 - 8 Feb 2003
//-----------------------------------------------------------------------------
static
int Leapfrog2(         //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      // I O  linked lists of labels
{
// Sampling interval: each atom is investigated about TIMESTEP times
static const double   TIMESTEP = 1.0000;
    int        Ndim     = Common->Ndim;
    double     cool     = Common->cool;
    ObjectStr* iObject  = &Objects[Oper->i];
    ObjectStr* jObject  = &Objects[Oper->j];
    ObjectStr* kObject  = &Objects[Oper->k];
    Node*      iLink    = &Links  [Oper->i];
    Node*      jLink    = &Links  [Oper->j];
    Node*      kLink    = &Links  [Oper->k];
    unsigned*  xTry     = iObject->xTry;
    unsigned*  Rand     = iObject->Rand;
    int        Natoms   = iObject->Natoms;
    int        del;      // selection
    Atom*      pAtom;    // selection
    double*    cube;     // selection
    double*    left;     // left  catalyst
    Atom*      pLeft;    // left  catalyst
    double*    right;    // right catalyst
    Atom*      pRight;   // right catalyst
    double*    trial;
    double     Ltry;
    double     accept;
    double     t;
    double*    swap;
    int        r;
    int        count;
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;

    for( count = 0; count <= Natoms * TIMESTEP; count++ )
    {
// Select original atom
        del = Rangrid(Rand, Natoms);
        cube = iObject->Cubes[del];
        pAtom = FindAtom(del, iLink);
// Left-definite and Right-definite neighbours
        pLeft = FindRight(pAtom->Label, jLink);
        if( ! pLeft )    pLeft = BeginLink(jLink);
        pLeft = PrevAtom(pLeft);
        if( ! pLeft )    pLeft = EndLink(jLink);
        left = jObject->Cubes[Storage(pLeft)];
        pRight = FindLeft(pAtom->Label, kLink);
        if( ! pRight )    pRight = EndLink(kLink);
        pRight = NextAtom(pRight);
        if( ! pRight )    pRight = BeginLink(kLink);
        right = kObject->Cubes[Storage(pRight)];

// Invert original cube through centre
        trial = iObject->Cubes[Natoms];
        for( r = 0; r < Ndim; r++ )
        {
            t = cube[r] - left[r];
            if( t < 0.0 )
                t += 1.0;
            t = right[r] - t;
            if( t < 0.0 )
                t += 1.0;
            trial[r] = t;
        }

// Ensure not already present
        CubetoLabel(xTry, iObject, Common);
        if( FindHere(xTry, iLink) )
            continue;
// Ensure same neighbours for detailed balance
        pAtom = FindRight(xTry, jLink);
        if( ! pAtom )    pAtom = BeginLink(jLink);
        pAtom = PrevAtom(pAtom);
        if( ! pAtom )    pAtom = EndLink(jLink);
        if( pAtom != pLeft )
            continue;
        pAtom = FindLeft(xTry, kLink);
        if( ! pAtom )    pAtom = EndLink(kLink);
        pAtom = NextAtom(pAtom);
        if( ! pAtom )    pAtom = BeginLink(kLink);
        if( pAtom != pRight )
            continue;

// Delete from source
        CALL( Delete1(del, &Ltry, NULL, Common, iObject, iLink) )       ++ CPU;
        accept = cool * Ltry + log(Randouble(Rand));
        SWAP(iObject->Cubes[Natoms], iObject->Cubes[Natoms-1])
        CALL( Try1(&Ltry, Common, iObject) )                            ++ CPU;
// Reject?
        if( CALLvalue == 0 || cool * Ltry < accept )
            SWAP(iObject->Cubes[Natoms], iObject->Cubes[Natoms-1])
        else                                                        ++ Success;
// Rebuild
        CALL( Insert1(Common, iObject, iLink) )                         ++ CPU;
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GuidedWalk                                     Parallel engine #9
//
// Purpose:   Geometrical MCMC engine.
//
// Method:    Tries to offset atoms along a line parallel to the displacement
//            vector between two neighbours from other objects.
//                         ________________________
//                        |                        |
//                        |                       .|
//                        |                    ... |
//                        |                  ..    |
//                        |               ...      |
//                        |             o.         |
//                        |          o.o           |
//                        |       ooo              |
//                        |     oX                 |   X = original point
//                        |  o..       R           |   R = catalyst
//                        |..                      |   L = catalyst
//                        |     L                  |   Line parallel to R-L
//                        |________________________|
//            Points "o" of line have same Hilbert neighbours L and R as have X
//            so these destinations are in detailed balance with X.
//
//            It is likely (but not guaranteed) that a starting set of atom
//            positions is capable of exploring the full space, but the number
//            of atoms is constant, so not explored.
//
//            Catalysts L and R are taken to be the left and right neighbours
//            of the original atom X, drawn from other object(s).
//            The line spans unit change in its most rapidly varying coordinate
//            and has 2^32 constituent points. The phase of the original point
//            X along the line is random.  Trial positions are chosen by slice
//            sampling along the line, so there will nearly always be movement.
//
//            Hopefully, the cloud of neighbours of X from different objects
//            follows the local shape of where X may lie, in which case the
//            engine can fill out "unmeasured" directions exponentially fast.
//
// Note:      ENSEMBLE should be at least 2 and Ndim should also be at least 2.
//
// History:   JS         14 Jan 2003, 8 Feb 2002, 19 Jul 2003
//-----------------------------------------------------------------------------
static
int GuidedWalk(        //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I(O) general information
ObjectStr* Objects,    // I O  ENSEMBLE of objects
Node*      Links)      // I O  linked lists of labels
{
// Sampling interval: each atom is investigated about TIMESTEP times
static const double   TIMESTEP = 1.0000;
static const double   Z        = (unsigned)(-1) + 1.0;         // 2^32
    int        Ndim     = Common->Ndim;
    double     cool     = Common->cool;
    unsigned*  offset   = Common->offset;
    ObjectStr* iObject  = &Objects[Oper->i];
    ObjectStr* jObject  = &Objects[Oper->j];
    ObjectStr* kObject  = &Objects[Oper->k];
    Node*      iLink    = &Links  [Oper->i];
    Node*      jLink    = &Links  [Oper->j];
    Node*      kLink    = &Links  [Oper->k];
    unsigned*  xTry     = iObject->xTry;
    unsigned*  xLabel   = iObject->xLabel;
    unsigned*  yLabel   = iObject->yLabel;
    unsigned*  Rand     = iObject->Rand;
    int        Natoms   = iObject->Natoms;
    int        del;      // selection
    Atom*      pAtom;    // selection
    double*    cube;     // old location
    double*    arrow;    // displacement arrow
    double*    trial;    // new location
    double*    left;     // left  catalyst
    double*    right;    // right catalyst
    Atom*      pLeft;    // left  catalyst
    Atom*      pRight;   // right catalyst
    double     Ltry;     // likelihood
    double     accept;   // acceptance level
    double     t;        // temporary
    double*    swap;     // exchange coordinates
    int        r;        // coordinate index
    int        m;        // dominant index of arrow
    unsigned   arc;      // length along displacement arrow
    unsigned   old;      // old position on displacement arrow
    unsigned   origin;   // arbitrary for slice sampling
    unsigned   mask;     // bits for slice sampling
    unsigned*  jLeft;    // left  boundary for object j
    unsigned*  jRight;   // right boundary for object j
    unsigned*  kLeft;    // left  boundary for object k
    unsigned*  kRight;   // right boundary for object k
    int        count;    // # atoms investigated
    int        CPU       = 1;
    int        Success   = 0;
    int        CALLvalue = 0;
#undef  STAIRCASE
#define STAIRCASE(i,s)  ((arrow[i] >= 0.0) ?  (unsigned)( s * arrow[i] ) \
                              : (unsigned)0 - (unsigned)((0-s) * arrow[i]))

    arrow = jObject->Cubes[NumAtoms(jLink)];
    for( count = 0; count <= Natoms * TIMESTEP; count++ )
    {
// Select original atom
        del = Rangrid(Rand, Natoms);
        pAtom = FindAtom(del, iLink);
        CopyLabel(xLabel, pAtom->Label, Ndim);

// Neighbours
        if( !(pLeft  = FindLeft (xLabel, jLink)) )   pLeft  = EndLink(jLink);
        if( !(pRight = FindRight(xLabel, kLink)) )   pRight = BeginLink(kLink);
        left  = jObject->Cubes[Storage(pLeft)];
        right = kObject->Cubes[Storage(pRight)];
        if( !(pAtom = NextAtom(pLeft)) )     pAtom = BeginLink(jLink);
        jLeft = pLeft->Label;    jRight = pAtom->Label;
        if( !(pAtom = PrevAtom(pRight)) )    pAtom = EndLink(kLink);
        kLeft = pAtom->Label;    kRight = pRight->Label;

// Displacement arrow
        t = 0.0;
        m = 0;
        for( r = 0; r < Ndim; r++ )
        {
            arrow[r] = right[r] - left[r];
            if( t < fabs(arrow[r]) )
            {
                t = fabs(arrow[r]);
                m = r;
            }
        }
        if( t == 0.0 )
            continue;
// Biggest reasonable slice-sampling mask
        arc = (unsigned)(Z * t);
        mask = 0;
        do  mask = mask + mask + 1;
        while( mask < arc );
// Normalised arrow
        t = arrow[m];
        for( r = 0; r < Ndim; r++ )
            arrow[r] /= t;

// Delete from source
        CALL( Delete1(del, &Ltry, NULL, Common, iObject, iLink) )       ++ CPU;
        accept = cool * Ltry + log(Randouble(Rand));
        SWAP(iObject->Cubes[Natoms-1], iObject->Cubes[Natoms])
        cube  = iObject->Cubes[Natoms];
        trial = iObject->Cubes[Natoms-1];

// Project original cube by displacement vector
        for( r = 0; r < Ndim; r++ )
            xLabel[r] = (unsigned)(Z * cube[r]);
        xLabel[m] -= offset[m];
        old = xLabel[m];
        for( r = 0; r < Ndim; r++ )
            xLabel[r] -= STAIRCASE(r, old);
        origin = Ranint(Rand);
        do
        {
            mask >>= 1;
            arc = ((old - origin) ^ (Ranint(Rand) & mask)) + origin;
            for( r = 0; r < Ndim; r++ )
                yLabel[r] = xLabel[r] + STAIRCASE(r, arc);
            yLabel[m] += offset[m];
            for( r = 0; r < Ndim; r++ )
                trial[r] = (0.5 + (double)yLabel[r]) / Z;
            CubetoLabel(xTry, iObject, Common);
// Ensure same neighbours for detailed balance
            if( Outside(xTry, jLeft, jRight, Ndim)
             || Outside(xTry, kLeft, kRight, Ndim) || FindHere(xTry, iLink) )
                continue;
            CALL( Try1(&Ltry, Common, iObject) )                        ++ CPU;
// Reject?
            if( CALLvalue > 0 && cool * Ltry >= accept )
                break;
        } while( mask );
// Rebuild: loop can only bottom out with trial = original cube
        CALL( Insert1(Common, iObject, iLink) )           ++Success;    ++ CPU;
    }
Exit:
    Oper->CPU     = CPU;
    Oper->Success = Success;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  MCMCexit                                      Parallel engine #10
//
// Purpose:   Empty Links
//
// History:   JS        16 Oct 2002, 1 Feb 2003
//-----------------------------------------------------------------------------
static
int MCMCexit(          //   O  0, or -ve error
OperStr*   Oper,       // I O  operation
CommonStr* Common,     // I    general information
ObjectStr* Objects,    //(I)   ENSEMBLE of objects
Node*      Links)      // I O  empty
{
    Node*  Link = &Links[Oper->i];
    Atom*  pAtom;
    int    CALLvalue = 0;

#if DEBUG                                      // Check Cube/Link consistency
    ObjectStr* Object = &Objects[Oper->i];
        CALLvalue = Object->Natoms;
        if( NumAtoms(Link) != Object->Natoms )
            return E_BAYESYS_DEBUG;
        for( Object->Natoms--; Object->Natoms >= 0; Object->Natoms-- )
        {
            CubetoLabel(Object->xTry, Object, Common);
            if( ! FindHere(Object->xTry, Link) )
                return E_BAYESYS_DEBUG;
        }
        Object->Natoms = CALLvalue;
#endif
    for( pAtom = BeginLink(Link); pAtom; pAtom = BeginLink(Link) )
        CALL( DelAtom(pAtom, Link) )
Exit:
    return (Common && Objects) ? CALLvalue : CALLvalue; //avoid compiler whinge
}

//=============================================================================
//
//            Event library for birth and death rates
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Birth
//
// Purpose:   d(birth) / dtime
//
// History:   JS         30 Sep 2002
//-----------------------------------------------------------------------------
static
double Birth(          //   O  birth rate
CommonStr* Common,     // I    general information
int        Natoms)     // I    # atoms
{
    int    MinAtoms = Common->MinAtoms;
    int    MaxAtoms = Common->MaxAtoms;
    double Alpha    = Common->Alpha;
    double r;      // birth rate

    if( MaxAtoms == MinAtoms )
        r = 0.0;
    else if( Alpha == 0.0 )            // MaxAtoms=0 excluded as improper
        r = (Natoms < MaxAtoms) ? Natoms + 1 - MinAtoms : 0;
    else if( Alpha > 0.0 )
    {
        r = Alpha;
        if( MaxAtoms )
            r *= (double)(MaxAtoms - Natoms) / (MaxAtoms - MinAtoms);
    }
    else  // Alpha < 0.0
    {
        r = (Natoms == MaxAtoms)
            ? 0.0
            : (Natoms - MinAtoms + 1.0) * Alpha / (Alpha - 1.0);
    }
    return r;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Death
//
// Purpose:   d(death) / dtime
//
// History:   JS         30 Sep 2002
//-----------------------------------------------------------------------------
static
double Death(          //   O  death rate
CommonStr* Common,     // I    general information
int        Natoms)     // I    # atoms
{
    int  MinAtoms = Common->MinAtoms;
    int  MaxAtoms = Common->MaxAtoms;

    if( MaxAtoms == MinAtoms )
        return 1.0;
    else
        return Natoms - MinAtoms;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SafeCube
//
// Purpose:   Ensure sufficient memory to add another atom into Object->Cubes
//
// History:   JS         2 Jan 2002
//-----------------------------------------------------------------------------
static
int SafeCube(         //   O  0, or -ve allocation error
CommonStr* Common,    // I    info and workspace
ObjectStr* Object)    // I O  sample object
{
    int     Ndim      = Common->Ndim;
    int     Valency   = Common->Valency;
    int     Nsize     = Ndim + Valency + 1;
    int     i, j, n;
    int     CALLvalue = 0;

    if( Ndim > 0 )
    if( Object->Natoms >= Object->Nstore - 1 )
    {
        n = Object->Natoms + Object->Natoms / 2 + 2;
        REALLOC(Object->Cubes, n, double*)
        for( j = Object->Nstore; j < n; j++ )
        {
            Object->Cubes[j] = NULL;   // (allows CALLOC to fail gracefully)
            CALLOC(Object->Cubes[j], Nsize, double)
            for( i = 0; i < Nsize; i++ )
                Object->Cubes[j][i] = 0.0;
        }
       Object->Nstore = n;
    }
Exit:
    return CALLvalue;
}

//=============================================================================
//
//                          Atom control
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Empty
//
// Purpose:   Set empty object.
//
// History:   JS         8 Feb 2002, 11 Oct 2003
//-----------------------------------------------------------------------------
static
int Empty(            //   O  0, or -ve error
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I O  sample object
{
    double Lhood;
    int    CALLvalue = 0;

    CALL( UserEmpty(&Object->Lhood, Common, Object) )
    if( Common->Valency )
    {
        FluxEmpty(&Lhood, Common, Object);
        Object->Lhood += Lhood;
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Try1
//
// Purpose:   Try inserting one new atom from Object->Cubes[Object->Natoms],
//            with NO update.
//
// History:   JS         8 Feb 2002
//-----------------------------------------------------------------------------
static
int Try1(             //   O  +ve = OK, 0 = DO NOT USE, -ve = error
double*    Ltry,      //   O  loglikelihood after adding trial atom
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I    sample object
{
    double Utry      = 0.0;
    double Ftry      = 0.0;
    int    CALLvalue = 0;

    CALL( UserTry1(&Utry, Common, Object) )
    if( Common->Valency )
        CALL( FluxTry1(&Ftry, Common, Object) )
    *Ltry = Object->Lhood + Utry + Ftry;
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Try2
//
// Purpose:   Try inserting two new atoms from Object->Cubes[Object->Natoms]
//            and Object->Cubes[Object->Natoms + 1],  with NO update.
//
// History:   JS         8 Feb 2002
//-----------------------------------------------------------------------------
static
int Try2(             //   O  +ve = OK, 0 = DO NOT USE, -ve = error
double*    Ltry1,     //   O  loglikelihood after adding first trial atom
double*    Ltry2,     //   O  loglikelihood after adding both trial atoms
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I    sample object
{
    double Utry1     = 0.0;
    double Ftry1     = 0.0;
    double Utry2     = 0.0;
    double Ftry2     = 0.0;
    int    CALLvalue = 0;

    CALL( UserTry2(&Utry1, &Utry2, Common, Object) )
    if( Common->Valency )
        CALL( FluxTry2(&Ftry1, &Ftry2, Common, Object) )
    *Ltry1 = Object->Lhood + Utry1 + Ftry1;
    *Ltry2 = Object->Lhood + Utry2 + Ftry2;
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Insert1
//
// Purpose:   Insert one new atom from Object->Cubes[Object->Natoms],
//            with full update.
//
// History:   JS         2 Jan 2002, 8 Feb 2002
//-----------------------------------------------------------------------------
static
int Insert1(          //   O  0, or -ve error
CommonStr* Common,    // I(O) general information
ObjectStr* Object,    // I O  sample object
Node*      Link)      // I O  linked list of labels
{
    Atom    atom[1];
    double  dLhood;
    int     CALLvalue = 0;

    atom->Label = Object->xTry;
    CubetoLabel(atom->Label, Object, Common);
    CALL( InsAtom(atom, Link) )
    if( ! CALLvalue )
        return E_BAYESYS_SYSERR;
    Object->Natoms ++;
    CALL( UserInsert1(&dLhood, Common, Object) )
    Object->Lhood += dLhood;
    if( Common->Valency )
    {
        CALL( FluxInsert1(&dLhood, Common, Object) )
        Object->Lhood += dLhood;
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Insert2
//
// Purpose:   Insert two new atoms from Object->Cubes[Object->Natoms] and
//            Object->Cubes[Object->Natoms + 1], with full update.
//
// History:   JS         2 Jan 2002, 8 Feb 2002
//-----------------------------------------------------------------------------
static
int Insert2(          //   O  0, or -ve error
CommonStr* Common,    // I(O) general information
ObjectStr* Object,    // I O  sample object
Node*      Link)      // I O  linked list of labels
{
    Atom    atom[1];
    double  dLhood;
    int     CALLvalue = 0;

    atom->Label = Object->xTry;
    CubetoLabel(atom->Label, Object, Common);
    CALL( InsAtom(atom, Link) )
    if( ! CALLvalue )
        return E_BAYESYS_SYSERR;
    Object->Natoms ++;

    CubetoLabel(atom->Label, Object, Common);
    CALL( InsAtom(atom, Link) )
    if( ! CALLvalue )
        return E_BAYESYS_SYSERR;
    Object->Natoms ++;

    CALL( UserInsert2(&dLhood, Common, Object) )
    Object->Lhood += dLhood;
    if( Common->Valency )
    {
        CALL( FluxInsert2(&dLhood, Common, Object) )
        Object->Lhood += dLhood;
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Delete1
//
// Purpose:   Delete one old atom from Object->Cubes[del], with full update,
//            and optional return of loglikelihoods after 1 or 2 re-insertions.
//
// History:   JS         2 Jan 2002, 8 Feb 2003
//-----------------------------------------------------------------------------
static
int Delete1(          //   O  0, or -ve error
int        del,       // I    serial # for deletion
double*    L1,        //  (O) loglikelihood after deletion then 1 insertion
double*    L2,        //  (O) loglikelihood after deletion then 2 insertions
CommonStr* Common,    // I    general information
ObjectStr* Object,    // I O  sample object
Node*      Link)      // I O  linked list of labels
{
    Atom*   pAtom;
    double* swap;
    double  dLhood;
    int     CALLvalue = 0;

    if( L1 )
        *L1 = Object->Lhood;
    pAtom = FindAtom(del, Link);
    CALL( DelAtom(pAtom, Link) )
    Object->Natoms --;
    SWAP(Object->Cubes[del], Object->Cubes[Object->Natoms])
    if( Common->Valency )
    {
        CALL( FluxDelete1(&dLhood, Common, Object) )
        Object->Lhood += dLhood;
    }
    CALL( UserDelete1(&dLhood, Common, Object) )
    Object->Lhood += dLhood;
    if( L1 && Common->Valency )
    {
        if( L2 )
            CALL( Try2(L1, L2, Common, Object) )
        else
            CALL( Try1(L1, Common, Object) )
    }
Exit:
    return CALLvalue;
}

//=============================================================================
//
//                  Label library for multi-bit labels
//      Label is multi-bit integer stored as vector of unsigned, [0] high.
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Topology
//
// Purpose:   Set "hypercube axes to extended-integer label" mapping,
//            preserving a degree of locality.
//
//   Mapping scales each coordinate is to an unsigned integer plus random
//   offset, reduced to range [0,2^32) by wraparound continuity.
//   Then either
//            Simple raster in random order, in which successive dimensions
//            take precedence.
//   or
//            Space-filling Hilbert curve in random orientation
//            (usually the preferred choice).
//
// History:   JS          2 Jan 2002, 16 Oct 2002
//-----------------------------------------------------------------------------
static
void Topology(
CommonStr* Common)    // I O  workspace
{
    int       Ndim    = Common->Ndim;       // I    # dimensions
    int*      permute = Common->permute;    //   O  permutation
    unsigned* offset  = Common->offset;     //   O  random
    unsigned* Rand    = Common->Rand;       // I O  random generator state

    if( Ndim <= 0 )
        return;
    RanLabel(offset, Ndim, Rand);
    Ranperm(Rand, Ndim, permute);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  CubetoLabel
//
// Purpose:   Translate hypercube axes to extended-integer label.
//            Inverse of LabeltoCube.
//
// History:   JS          25 April 2001, 2 Jan 2002, 16 Oct 2002
//-----------------------------------------------------------------------------
static
void CubetoLabel(
unsigned*  Label,     //   O  extended-integer label               [Ndim]
ObjectStr* Object,    // I    object containing input cube
CommonStr* Common)    // I    general information
{
static const double Z = (unsigned)(-1) + 1.0;          // 2^32
    int       method  = Common->Method;                // I    mapping switch
    int       Ndim    = Common->Ndim;                  // I    # dimensions
    int       Nbits   = Common->Nbits;                 // I    # bits per word
    int*      permute = Common->permute;               // I    permutation
    unsigned* offset  = Common->offset;                // I    random
    unsigned* work    = Object->work;                  //  (O) workspace
    double*   Cube    = Object->Cubes[Object->Natoms]; // I    hypercube[Ndim]
    int       i;

    if( Ndim <= 0 )
        return;
    for( i = 0; i < Ndim; i++ )                 // Read from "double" in (0,1)
        work[i] = (unsigned)(Z * fmod(fabs(Cube[permute[i]]), 1.)) + offset[i];
    if( method & 1 )
        AxestoLine(Label, work, Nbits, Ndim);   // Hilbert curve
    else
        for( i = 0; i < Ndim; i++ )             // Simple raster
            Label[i] = work[i];
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  LabeltoCube
//
// Purpose:   Translate extended-integer label to hypercube axes.
//            Inverse of CubetoLabel.
//
// History:   JS          25 April 2001, 2 Jan 2002, 6 Jan 2002, 16 Oct 2002
//-----------------------------------------------------------------------------
static
void LabeltoCube(
ObjectStr* Object,    //   O  object containing output cube
unsigned*  Label,     // I    extended-integer label               [Ndim]
CommonStr* Common)    // I(O) workspace
{
static const double Z = (unsigned)(-1) + 1.0;          // 2^32
    int       method  = Common->Method;                // I    mapping switch
    int       Ndim    = Common->Ndim;                  // I    # dimensions
    int       Nbits   = Common->Nbits;                 // I    # bits per word
    int*      permute = Common->permute;               // I    permutation
    unsigned* offset  = Common->offset;                // I    random
    unsigned* work    = Object->work;                  //  (O) workspace
    double*   Cube    = Object->Cubes[Object->Natoms]; //   O  hypercube[Ndim]
    unsigned  temp;
    int       i;

    if( Ndim <= 0 )
        return;
    if( method & 1 )
        LinetoAxes(work, Label, Nbits, Ndim);   // Hilbert curve
    else
        for( i = 0; i < Ndim; i++ )             // Simple raster
            work[i] = Label[i];
    for( i = 0; i < Ndim; i++ )
    {
        temp = work[i] - offset[i];             // force unsigned arith here
        Cube[permute[i]] = (temp + 0.5) / Z;    // Write to "double" in (0,1)
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  RanLabel
//
// Purpose:   Set     r = random
//
// History:   JS          25 Apr 2001
//-----------------------------------------------------------------------------
static
void RanLabel(
unsigned*  r,      //   O  value                  [Ndim]
int        Ndim,   // I    dimension
unsigned*  Rand)   // I O  random generator state
{
    int j;
    for( j = 0; j < Ndim; j++ )
        r[j] = Ranint(Rand);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  CopyLabel
//
// Purpose:   Set     w = u
//
// History:   JS          28 Jan 2002
//-----------------------------------------------------------------------------
static
void CopyLabel(
unsigned*  w,      //   O  w = u                   [Ndim]
unsigned*  u,      // I    can be overwritten      [Ndim]
int        Ndim)   // I    dimension
{
    int  i;
    for( i = 0; i < Ndim; i++ )
    {
        w[i] = u[i];
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  NewLabel
//
// Purpose:   Set     New = ((Old - Origin) ^ Random) + Origin
//                   nbits   nbits  nbits     rbits     nbits
//
// History:   JS          2 Jan, 12 Jan 2002
//-----------------------------------------------------------------------------
static
void NewLabel(
unsigned*  New,    //   O  new label              [Ndim]
unsigned*  Old,    // I    old label              [Ndim]
unsigned*  Origin, // I    random origin          [Ndim]
int        Ndim,   // I    dimension
int        rbits,  // I    # low-order bits to be randomised
int        Nbits,  // I    # bits per word          (nbits = Ndim * Nbits)
unsigned*  Rand)   // I O  random generator state
{
    int  shift;
    int  i;

    if( rbits < 1 )
    {
        for( i = 0; i < Ndim; i++ )
            New[i] = Old[i];
    }
    else
    {
        SubLabel(New, Old, Origin, Ndim);
        shift = Nbits - 1 - (rbits - 1) % Nbits;
        i = Ndim - 1 - (rbits - 1) / Nbits;
        New[i] = New[i] ^ ((unsigned)Ranint(Rand) >> shift);
        for( i++; i < Ndim; i++ )
            New[i] = New[i] ^ (unsigned)Ranint(Rand);
        AddLabel(New, New, Origin, Ndim);
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Outside
//
// Purpose:   Is point strictly inside interval, or is it outside or on edge?
//            Interval [Left,Right] may wraparound.
//
// History:   JS          2 Jan 2002
//-----------------------------------------------------------------------------
static
int Outside(       //   O  0 = strictly inside, 1 = outside or on edge
unsigned*  Where,  // I    where is this point?     [Ndim]
unsigned*  Left,   // I    left edge                [Ndim]
unsigned*  Right,  // I    right edge               [Ndim]
int        Ndim)   // I    dimension
{
    if( CmpLabel(Right, Left, Ndim) > 0 ) //  [0...Left........Right....2^32-1]
    {                                     //           .inside.
        if( (CmpLabel(Where, Left,  Ndim) <= 0)
         || (CmpLabel(Where, Right, Ndim) >= 0) )
            return 1;
    }
    else                                  //  [0........Right....Left...2^32-1]
    {                                     //    .inside.             .inside.
        if( (CmpLabel(Where, Left,  Ndim) <= 0)
         && (CmpLabel(Where, Right, Ndim) >= 0) )
            return 1;
    }
    return 0;
}

//=============================================================================
//
//                             MassInf procedures
//
//               MAIN
//           ______|_______
//          |              |
//          |   BayeSys3   |
//          |______________|
//                 |      |
//                 |      |
//                 |      |_______________________________________________
//                 |     |                      MassInf                   |
//                 |     | FluxEmpty               |                      |
//                 |     | FluxTry1   ,FluxTry2    | FluxAlloc, FluxInit  |
//                 |     | FluxInsert1,FluxInsert2 | FluxCalib0,FluxCalib |
//                 |     | FluxDelete1             | FluxFree             |
//                 |     |_________________________|______________________|
//                 |                   |
//                 |                   |
//     ____________|___________________|_______
//    |                       USER             |
//    | UserMonitor             |              |
//    | UserEmpty               |              |
//    | UserTry1   ,UserTry2    |   UserFoot   |
//    | UserInsert1,UserInsert2 |              |
//    | UserDelete1             |              |
//    |_________________________|______________|
//=============================================================================
//
//                      MassInf ancillary library
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxAlloc
//
// Purpose:   Allocate memory for MassInf
//
// History:   JS       10 Feb 2003, 15 Aug 2003
//-----------------------------------------------------------------------------
static
int FluxAlloc(           //   O  0, or -ve error
CommonStr* Common,       // I    general info
ObjectStr* Objects)      // I O  allocate memory
{
    int        ENSEMBLE  = Common->ENSEMBLE;
    int        MassInf   = Common->MassInf;
    int        Valency   = Common->Valency;
    int        Ndata     = Common->Ndata;
    int        k;
    int        CALLvalue = 0;

    Common->Counts = NULL;
    if( MassInf >= 100 )
        CALLOC(Common->Counts, Ndata, int)
    for( k = 0; k < ENSEMBLE; k++ )
    {
        Objects[k].Mock    = NULL;
        Objects[k].g1      = NULL;
        Objects[k].g2      = NULL;
        Objects[k].A11     = NULL;
        Objects[k].A12     = NULL;
        Objects[k].A22     = NULL;
        Objects[k].Xindex  = NULL;
        Objects[k].nbits   = NULL;
        Objects[k].nbitx   = NULL;
        Objects[k].ibits   = NULL;
        Objects[k].ibitx   = NULL;
        Objects[k].zbits   = NULL;
        Objects[k].zbitx   = NULL;
        Objects[k].Foot    = NULL;
        Objects[k].flags   = NULL;
        CALLOC(Objects[k].Mock,   Ndata,   double)
        CALLOC(Objects[k].g1,     Valency, double)
        CALLOC(Objects[k].g2,     Valency, double)
        CALLOC(Objects[k].A11,    Valency, double)
        CALLOC(Objects[k].A12,    Valency, double)
        CALLOC(Objects[k].A22,    Valency, double)
        CALLOC(Objects[k].Xindex, Valency, int)
        CALLOC(Objects[k].nbits,  Valency, int)
        CALLOC(Objects[k].nbitx,  Valency, int)
        if( k == 0 || PARALLEL )
        {
            CALLOC(Objects[k].ibits,  Ndata, int)
            CALLOC(Objects[k].ibitx,  Ndata, int)
            CALLOC(Objects[k].zbits,  Ndata, double)
            CALLOC(Objects[k].zbitx,  Ndata, double)
            CALLOC(Objects[k].Foot ,  Ndata, double)
            CALLOC(Objects[k].flags,  Ndata, int)
        }
        else
        {
            Objects[k].ibits = Objects[0].ibits;
            Objects[k].ibitx = Objects[0].ibitx;
            Objects[k].zbits = Objects[0].zbits;
            Objects[k].zbitx = Objects[0].zbitx;
            Objects[k].Foot  = Objects[0].Foot;
            Objects[k].flags = Objects[0].flags;
        }
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxInit
//
// Purpose:   Initialise system for MassInf
//
// History:   JS       10 Feb 2003
//-----------------------------------------------------------------------------
static
int FluxInit(            //   O  0, or -ve error
CommonStr* Common,       // I    general info
ObjectStr* Objects)      // I O  initialise object info
{
    int        ENSEMBLE  = Common->ENSEMBLE;
    int        Ndata     = Common->Ndata;
    int        j, k;

    for( k = 0; k < ENSEMBLE; k++ )
        for( j = 0; j < Ndata; j++ )
        {
            Objects[k].flags[j] = -1;
            Objects[k].Foot[j] = 0.0;
        }
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxCalib0
//
// Purpose:   Set internal hyperparameters
//            FluxUnit0 in Common and FluxUnit in Objects.
//
// History:    JS        1 Feb 2003, 11 Jul 2003, 12 Aug 2003, 20 Aug 2003
//                      12 Sep 2003
//-----------------------------------------------------------------------------
static
int FluxCalib0(          //   O  0, or -ve error
CommonStr* Common,       // I O  set hyperparameters
ObjectStr* Objects,      //   O  set hyperparameters
Node*      Links)        //  (O) empty on exit
{
    int        ENSEMBLE = Common->ENSEMBLE;     // I    # objects
    int        MassInf  = Common->MassInf;      // I    type of data
    int        Ndata    = Common->Ndata;        // I    # data
    double*    Data     = Common->Data;         // I    data            [Ndata]
    double*    Acc      = Common->Acc;          // I    accuracies      [Ndata]
    double*    Mock     = Objects->Mock;        //  (O) mock data       [Ndata]
    int*       Counts   = Common->Counts;       //  (O) (Poisson data)  [Ndata]
    int        N        = 12;     // adequate # samples
    double     Datanorm = 0.0;    // ||data||
    double     Mocknorm = 0.0;    // ||mockdata||
    OperStr    Operations[1];     // set prior
    int        i, k;              // counter
    int        CALLvalue = 0;

// Poisson Counts = 0 at cool=0
    if( MassInf >= 100 )
        for( k = 0; k < Ndata; k++ )
            Counts[k] = 0;
// Common parameter(s) from preliminary peek at data, if necessary
    if( Common->FluxUnit0 == 0.0 )
    {
        Objects->FluxUnit  = 1.0;      // Preliminary flux unit
        Operations->engine = 0;        // SetPrior
        Operations->iType  = -1;    Operations->i = 0;
        Operations->jType  = 0;
        Operations->kType  = 0;
        if( MassInf >= 100 )           // Poisson
        {
            for( i = 0; i < Ndata; i++ )
                Datanorm += Data[i] + Acc[i];
            for( k = 0; k < N; k++ )
            {
                CALL( DoOperations(Operations, Common, Objects, Links, 1) )
                for( i = 0; i < Ndata; i++ )
                    Mocknorm += Mock[i] + Acc[i];
            }
            Mocknorm /= N;
        }
        else                           // Gaussian
        {
            for( i = 0; i < Ndata; i++ )
                Datanorm += Data[i] * Acc[i] * Acc[i] * Data[i];
            Datanorm = sqrt(Datanorm);
            for( k = 0; k < N; k++ )
            {
                CALL( DoOperations(Operations, Common, Objects, Links, 1) )
                for( i = 0; i < Ndata; i++ )
                    Mocknorm += Mock[i] * Acc[i] * Acc[i] * Mock[i];
            }
            Mocknorm = sqrt(Mocknorm / N);
        }
        Common->FluxUnit0 = -( (Datanorm > 0.0 && Mocknorm > 0.0)
                               ? Datanorm/Mocknorm : 1.0 );
    }
// Object parameters
    if( Common->FluxUnit0 > 0.0 )
        for( k = 0; k < ENSEMBLE; k++ )
            Objects[k].FluxUnit = Common->FluxUnit0;
    else
        for( k = 0; k < ENSEMBLE; k++ )
            Objects[k].FluxUnit = Ran1posX(Objects[k].Rand,
                                           1.0 / (-Common->FluxUnit0), 0.0);
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxCalib
//
// Purpose:   Reset after cooling.
//
//            Common Counts = (int)(cool * Poisson data)
//
//            Object FluxUnit q  =  mean prior flux per atom
//            sampled according to
//                Pr(q | F*,Data)  ~  Prior(q) Likelihood(Data | qF*)
//            where
//                Prior(q)  ~  q exp(-q/q0),    q0 = prior hyperparameter
//                F* = Mock/q = current mock signal as generated from unit q=1
//
// History:   JS  2 Jan 2002 - 20 Aug 2003
//-----------------------------------------------------------------------------
static
void FluxCalib(
CommonStr* Common,       // I    read hyperparameters
ObjectStr* Object)       // I O  reset hyperparameters
{
    int        MassInf  = Common->MassInf;      // I    type of data
    double     cool     = Common->cool;         // I    Annealing parameter
    double     q0       = Common->FluxUnit0;    // I    Prior hyperparameter
    int        Ndim     = Common->Ndim;         // I    # coordinates
    int        Valency  = Common->Valency;      // I    # fluxes
    int        Ndata    = Common->Ndata;        // I    # data
    double*    Data     = Common->Data;         // I    data            [Ndata]
    double*    Acc      = Common->Acc;          // I    accuracies      [Ndata]
    int*       Counts   = Common->Counts;       //  (O) (Poisson data)  [Ndata]
    double*    Mock     = Object->Mock;         // I O  mock data       [Ndata]
    int        Natoms   = Object->Natoms;       // I    # atoms
    double     q        = Object->FluxUnit;     // I    flux unit: F*=Mock/q
    double*    Flux;                       // fluxes                 [Valency]
    double     a;                          // autoscaled accuracy
    double     FF;                         // Mock.Acc.Acc.Mock
    double     FD;                         // Mock.Acc.Acc.Data
    double     g;                          // gradient   F*.Acc.Acc.Data
    double     A;                          // curvature  F*.Acc.Acc.F*
    double     s;                          // rescale by (new q)/(old q)
    int        i, j;                       // counters

// Poisson data ?
    if( MassInf >= 100 )
    {
        a = 0.5 + cool * (Data[0] + Acc[0]);
        Counts[0] = (int)a;
        for( i = 1; i < Ndata; i++ )
        {
            a += cool * (Data[i] + Acc[i]) - Counts[i-1];
            Counts[i] = (int)a;
        }
    }
// New flux unit
    FF = FD = 0.0;
    for( i = 0; i < Ndata; i++ )
    {
        if( MassInf >= 100 )         // Poisson (approx is OK here)
        {
            FF += Mock[i] * Mock[i] / (Data[i] + Acc[i]);
            FD += Mock[i] * Data[i] / (Data[i] + Acc[i]);
        }
        else                         // Gaussian
        {
            FF += Mock[i] * Acc[i] * Acc[i] * Mock[i];
            FD += Mock[i] * Acc[i] * Acc[i] * Data[i];
        }
    }
    g = cool * FD / q;
    A = cool * FF / (q * q);
    if( q0 <= 0.0 )
        Object->FluxUnit = Ran1posX(Object->Rand, 1.0 / fabs(q0) - g, A);
// Rescale observables
    s = Object->FluxUnit / q;
    for( i = 0; i < Ndata; i++ )
        Mock[i] *= s;
    for( i = 0; i < Natoms; i++ )
    {
        Flux = Object->Cubes[i] + Ndim;
        for( j = 0; j < Valency; j++ )
            Flux[j] *= s;
    }
    Object->Lhood = (MassInf >= 100 )
                   ? PoissLhood(Ndata, Mock, Data, Acc)
                   : GaussLhood(Ndata, Mock, Data, Acc);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxFree
//
// Purpose:   Free memory for MassInf
//
// History:   JS       10 Feb 2003
//-----------------------------------------------------------------------------
static
void FluxFree(
CommonStr* Common,       // I    general info
ObjectStr* Objects)      // I(O) free memory
{
    int        ENSEMBLE  = Common->ENSEMBLE;
    int        MassInf   = Common->MassInf;
    int        k;             // object counter

    for( k = 0; k < ENSEMBLE; k++ )
    {
        if( k == 0 || PARALLEL )
        {
            FREE(Objects[k].flags)
            FREE(Objects[k].Foot)
            FREE(Objects[k].zbitx)
            FREE(Objects[k].zbits)
            FREE(Objects[k].ibitx)
            FREE(Objects[k].ibits)
        }
        FREE(Objects[k].nbitx)
        FREE(Objects[k].nbits)
        FREE(Objects[k].Xindex)
        FREE(Objects[k].A22)
        FREE(Objects[k].A12)
        FREE(Objects[k].A11)
        FREE(Objects[k].g2)
        FREE(Objects[k].g1)
        FREE(Objects[k].Mock)
    }
    if( MassInf >= 100 )
        FREE(Common->Counts)
}

//=============================================================================
//
//                      MassInf outer procedures
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxEmpty
//
// Purpose:   Set Lhood = logLikelihood(fluxes) from Data,  for empty sample
//            with 0 atoms in Object, and initialise its other information.
//
//            Lhood := log L(empty)
//
//            The number 0 has already been placed in Object->Natoms.
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 18 Aug 2003
//-----------------------------------------------------------------------------
static
void FluxEmpty(
double*    Lhood,     //   O  loglikelihood
CommonStr* Common,    // I O  general information
ObjectStr* Object)    // I(O) sample object
{
    int      MassInf = Common->MassInf;    // MassInf switches
    int      Ndata   = Common->Ndata;      // I    # data
    double*  Data    = Common->Data;       // I    Data            [Ndata]
    double*  Acc     = Common->Acc;        // I    Accuracies      [Ndata]
    double*  Mock    = Object->Mock;       //   O  Mock data       [Ndata]
    int      k;                            //      data counter

    for( k = 0; k < Ndata; k++ )
        Mock[k] = 0.0;
    *Lhood = (MassInf >= 100)
            ? PoissLhood(Ndata, Mock, Data, Acc)
            : GaussLhood(Ndata, Mock, Data, Acc);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxTry1
//
// Purpose:   d(logLikelihood(fluxes)) after supposedly adding one new atom
//            to Object, without any side effects on it.
//                                                                      cool
//              dLtry := (1/cool) DELTA log INTEGRAL dPrior(z) L(...,x,z)
//
//            The new atom x has already been placed after the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms ].
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 12 Aug 2003, 18 Aug 2003
//-----------------------------------------------------------------------------
static
int FluxTry1(         //   O  +ve = OK, 0 = DO NOT USE, -ve = error abort
double*    dLtry,     //   O  d(logLikelihood) from inserting one atom
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I    sample object (DO NOT UPDATE Object info)
{
    int       MassInf = Common->MassInf;    // MassInf switches
    int       Valency = Common->Valency;    // # fluxes per atom
    double    ProbON  = Common->ProbON;     // Pr(individual flux != 0)
    double    cool    = Common->cool;       // Annealing parameter
    double*   Data    = Common->Data;       // Data                     [Ndata]
    double*   Acc     = Common->Acc;        // Accuracies               [Ndata]
    int*      Counts  = Common->Counts;     // Annealed Poisson counts  [Ndata]
    double*   Mock    = Object->Mock;       // Mock data                [Ndata]
    double*   g1      = Object->g1;         // grad chisquared        [Valency]
    double*   A11     = Object->A11;        // grad grad chisquared   [Valency]
    double    q       = Object->FluxUnit;   // Unit of flux
    int*      nbits   = Object->nbits;      // # fragments            [Valency]
    int*      ibits   = Object->ibits;      // Fragment positions     [<=Ndata]
    double*   zbits   = Object->zbits;      // Fragment quantities    [<=Ndata]
    int*      flags   = Object->flags;      // flags for overlap        [Ndata]
    double*   Cube;                        // Position of new atom      [Ndim]
    int       code;
    int       CALLvalue = 0;

    *dLtry = 0.0;                           // default output
// Footprint
    Cube = Object->Cubes[Object->Natoms];
    CALL( AtomBits(Cube, Common, ibits, zbits, nbits, flags) )
    if( CALLvalue == 0 )
        goto Exit;
    code = CALLvalue;

// Gaussian
    if( MassInf < 100 )
    {
// Construct quadratic chisquared
        SetGrad1(Mock, Data, Acc, Valency, nbits, ibits, zbits, g1,A11);
// Marginal likelihood
        *dLtry = GaussTry1(MassInf, cool, q, ProbON, Valency, g1, A11);
    }

// Poisson
    else
        CALL( PoissTry1(MassInf, cool, q, ProbON, Mock, Data, Acc, Counts,
                        Valency, nbits, ibits, zbits, dLtry) )
    return code;
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxTry2
//
// Purpose:   d(logLikelihood(fluxes)) after supposedly adding one, then
//            two new atoms to Object, without any side effects on it.
//
//                                                                cool
//    dLtry1 := (1/cool) DELTA log INTEGRAL dPrior(z1) L(...,x1,z1)
//                                                                         cool
//    dLtry2 := (1/cool) DELTA log INTEGRAL dPrior(z1,z2) L(...,x1,z1,x2,z2)
//
//            The new atoms x1,x2 have already been put after the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms ] and
//                               at Object->Cubes[ Object->Natoms + 1 ].
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 12 Aug 2003, 18 Aug 2003
//-----------------------------------------------------------------------------
static
int FluxTry2(         //   O  +ve = OK, 0 = DO NOT USE, -ve = error abort
double*    dLtry1,    //   O  d(logLikelihood) from inserting first atom
double*    dLtry2,    //   O  d(logLikelihood) from inserting both atoms
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I    sample object (DO NOT UPDATE Object info)
{
    int       MassInf = Common->MassInf;    // MassInf switches
    int       Valency = Common->Valency;    // # fluxes per atom
    double    ProbON  = Common->ProbON;     // Pr(individual flux != 0)
    double    cool    = Common->cool;       // Annealing parameter
    double*   Data    = Common->Data;       // Data                     [Ndata]
    double*   Acc     = Common->Acc;        // Accuracies               [Ndata]
    int*      Counts  = Common->Counts;     // Annealed Poisson counts  [Ndata]
    double*   Mock    = Object->Mock;       // Mock data                [Ndata]
    double*   Foot    = Object->Foot;       // Footprint workspace (=0) [Ndata]
    int*      flags   = Object->flags;      // Valency identifiers (-1) [Ndata]
    double*   g1      = Object->g1;         // grad chisquared        [Valency]
    double*   g2      = Object->g2;         // grad chisquared        [Valency]
    double*   A11     = Object->A11;        // grad grad chisquared   [Valency]
    double*   A12     = Object->A12;        // grad grad chisquared   [Valency]
    double*   A22     = Object->A22;        // grad grad chisquared   [Valency]
    int*      Xindex  = Object->Xindex;     // index to crossterms    [Valency]
    double    q       = Object->FluxUnit;   // Unit of flux
    int*      nbits   = Object->nbits;      // # fragments            [Valency]
    int*      nbitx   = Object->nbitx;      // # fragments            [Valency]
    int*      ibits   = Object->ibits;      // Fragment positions     [<=Ndata]
    int*      ibitx   = Object->ibitx;      // Fragment positions     [<=Ndata]
    double*   zbits   = Object->zbits;      // Fragment quantities    [<=Ndata]
    double*   zbitx   = Object->zbitx;      // Fragment quantities    [<=Ndata]
    double*   Cube1;                       // Position of new atom      [Ndim]
    double*   Cube2;                       // Position of new atom      [Ndim]
    int       code;
    int       CALLvalue = 0;

    *dLtry1 = *dLtry2 = 0.0;                // default output
// Footprints
    Cube1 = Object->Cubes[Object->Natoms];
    Cube2 = Object->Cubes[Object->Natoms + 1];
    CALL( AtomBits(Cube1, Common, ibits, zbits, nbits, flags) )
    if( CALLvalue == 0 )
        goto Exit;
    CALL( AtomBits(Cube2, Common, ibitx, zbitx, nbitx, flags) )
    if( CALLvalue == 0 )
        goto Exit;
    code = CALLvalue;
    SetIndex(Valency, nbits, ibits, nbitx, ibitx, flags, Xindex);

// Gaussian
    if( MassInf < 100 )
    {
// Construct quadratic chisquared
        SetGrad2(Mock, Data, Acc, Foot, Valency,
                 nbits, ibits, zbits, nbitx, ibitx, zbitx,
                 g1, g2, A11, A12, A22, Xindex);
// Marginal likelihoods
        GaussTry2(MassInf, cool, q, ProbON, Valency,
                  g1, g2, A11, A12, A22, Xindex, dLtry1, dLtry2);
    }

// Poisson
    else
        CALL( PoissTry2(MassInf, cool, q, ProbON, Mock, Data,Acc,Counts, Foot,
                        Valency, nbits, ibits, zbits, nbitx, ibitx, zbitx,
                        Xindex, dLtry1, dLtry2) )
    return code;
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxInsert1
//
// Purpose:   Insert 1 new atom with fluxes into Object, keeping it up-to-date,
//            and set d(loglikelihood(fluxes)).
//                                              cool
//            Sample z from  Prior(z) L(...,x,z)
//
//            and set   dL := DELTA logL(...,x,z)  at the sampled z.
//
//            The new atom x has already been placed in the last location of
//            the updated Atoms list, at Object->Cubes[ Object->Natoms - 1 ].
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 12 Aug 2003, 18 Aug 2003
//-----------------------------------------------------------------------------
static
int FluxInsert1(        //   O  0, or -ve return code
double*    dL,          //   O  d(loglikelihood(fluxes))
CommonStr* Common,      // I O  general information
ObjectStr* Object)      // I O  updated sample object
{
    int       MassInf = Common->MassInf;    // MassInf switches
    int       Ndim    = Common->Ndim;       // dimension
    int       Valency = Common->Valency;    // # fluxes per atom
    double    ProbON  = Common->ProbON;     // Pr(individual flux != 0)
    double    cool    = Common->cool;       // Annealing parameter
    double*   Data    = Common->Data;       // Data                     [Ndata]
    double*   Acc     = Common->Acc;        // Accuracies               [Ndata]
    int*      Counts  = Common->Counts;     // Annealed Poisson counts  [Ndata]
    double*   Mock    = Object->Mock;       // Mock data                [Ndata]
    int*      flags   = Object->flags;      // Valency identifiers (-1) [Ndata]
    double*   g1      = Object->g1;         // grad chisquared        [Valency]
    double*   A11     = Object->A11;        // grad grad chisquared   [Valency]
    int       reset   = Object->reset;      // Reset inserted MassInf fluxes?
    double    q       = Object->FluxUnit;   // Unit of flux
    unsigned* Rand    = Object->Rand;       // Random generator
    int*      nbits   = Object->nbits;      // # fragments            [Valency]
    int*      ibits   = Object->ibits;      // Fragment positions     [<=Ndata]
    double*   zbits   = Object->zbits;      // Fragment quantities    [<=Ndata]
    double*   Cube;                        // Position of new atom      [Ndim]
    double*   Flux;                         // Flux of new atom       [Valency]
    int       CALLvalue = 0;

    *dL = 0.0;                              // default output
// Footprint
    Cube = Object->Cubes[Object->Natoms - 1];
    Flux = Cube + Ndim;
    CALL( AtomBits(Cube, Common, ibits, zbits, nbits, flags) )
    if( CALLvalue == 0 )
        return E_BAYESYS_SYSERR;
    CALL( Overlap1(nbits, ibits, flags, Valency) )     //  Abort?

// Gaussian
    if( MassInf < 100 )
    {
// Construct quadratic chisquared
       SetGrad1(Mock, Data, Acc, Valency, nbits, ibits, zbits, g1, A11);
// Optionally reset Flux
        if( reset )
            GaussInsert1(MassInf, Rand, cool, q, ProbON, Valency,
                        g1, A11, Flux);
// Modify Lhood
        *dL = GaussLhood1(Valency, Flux, g1, A11, +1);
    }

// Poisson
    else
    {
        if( reset )
            CALL( PoissInsert1(MassInf, Rand, cool, q, ProbON, Mock, Data, Acc,
                               Counts, Valency, nbits, ibits, zbits, Flux) )
        *dL = PoissLhood1(Mock, Data, Acc,
                               Valency, nbits, ibits, zbits, Flux, +1);
    }

// Modify mock data
    InsBits(Mock, Flux, nbits, ibits, zbits, Valency);
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxInsert2
//
// Purpose:   Insert 2 new atoms with fluxes into Object, keeping it up-to-date
//            and set d(loglikelihood(fluxes)).
//                                                              cool
//            Sample z1,z2 from  Prior(z1,z2) L(...,x1,z1,x2,z2)
//
//            and set  dL := DELTA logL(...,x1,z1,x2,z2)  at the sampled z1,z2.
//
//            The new atoms x1,x2 have already been put in the last locations
//            of the updated Atoms list at Object->Cubes[ Object->Natoms - 2 ]
//                                  and at Object->Cubes[ Object->Natoms - 1 ]
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 12 Aug 2003, 18 Aug 2003
//-----------------------------------------------------------------------------
static
int FluxInsert2(        //   O  0, or -ve return code
double*    dL,          //   O  d(loglikelihood(fluxes))
CommonStr* Common,      // I O  general information
ObjectStr* Object)      // I O  updated sample object
{
    int       MassInf = Common->MassInf;    // MassInf switches
    int       Ndim    = Common->Ndim;       // dimension
    int       Valency = Common->Valency;    // # fluxes per atom
    double    ProbON  = Common->ProbON;     // Pr(individual flux != 0)
    double    cool    = Common->cool;       // Annealing parameter
    double*   Data    = Common->Data;       // Data                     [Ndata]
    double*   Acc     = Common->Acc;        // Accuracies               [Ndata]
    int*      Counts  = Common->Counts;     // Annealed Poisson counts  [Ndata]
    double*   Mock    = Object->Mock;       // Mock data                [Ndata]
    double*   Foot    = Object->Foot;       // Footprint workspace (=0) [Ndata]
    int*      flags   = Object->flags;      // Valency identifiers (=-1)[Ndata]
    double*   g1      = Object->g1;         // grad chisquared        [Valency]
    double*   g2      = Object->g2;         // grad chisquared        [Valency]
    double*   A11     = Object->A11;        // grad grad chisquared   [Valency]
    double*   A12     = Object->A12;        // grad grad chisquared   [Valency]
    double*   A22     = Object->A22;        // grad grad chisquared   [Valency]
    int*      Xindex  = Object->Xindex;     // index to crossterms    [Valency]
    double    q       = Object->FluxUnit;   // Unit of flux
    unsigned* Rand    = Object->Rand;       // Random generator state
    int*      nbits   = Object->nbits;      // # fragments            [Valency]
    int*      nbitx   = Object->nbitx;      // # fragments            [Valency]
    int*      ibits   = Object->ibits;      // Fragment positions     [<=Ndata]
    int*      ibitx   = Object->ibitx;      // Fragment positions     [<=Ndata]
    double*   zbits   = Object->zbits;      // Fragment quantities    [<=Ndata]
    double*   zbitx   = Object->zbitx;      // Fragment quantities    [<=Ndata]
    double*   Cube1;                       // Position of new atom      [Ndim]
    double*   Cube2;                       // Position of new atom      [Ndim]
    double*   Flux1;                        // Flux of new atom       [Valency]
    double*   Flux2;                        // Flux of new atom       [Valency]
   int       CALLvalue = 0;

    *dL = 0.0;                              // default output
// Footprints
    Cube1 = Object->Cubes[Object->Natoms - 1];
    Cube2 = Object->Cubes[Object->Natoms - 2];
    Flux1 = Cube1 + Ndim;
    Flux2 = Cube2 + Ndim;
    CALL( AtomBits(Cube2, Common, ibitx, zbitx, nbitx, flags) )
    if( CALLvalue == 0 )
        return E_BAYESYS_SYSERR;
    CALL( AtomBits(Cube1, Common, ibits, zbits, nbits, flags) )
    if( CALLvalue == 0 )
        return E_BAYESYS_SYSERR;
    CALL( Overlap2(nbits, ibits, nbitx, ibitx, flags, Xindex, Valency) )

// Gaussian
    if( MassInf < 100 )
    {
// Construct quadratic chisquared
        SetGrad2(Mock, Data, Acc, Foot, Valency,
                 nbits, ibits, zbits, nbitx, ibitx, zbitx,
                 g1, g2, A11, A12, A22, Xindex);
// Reset Fluxes
        GaussInsert2(MassInf, Rand, cool, q, ProbON, Valency,
                     g1, g2, A11, A12, A22, Xindex, Flux1, Flux2);
// Modify Lhood
       *dL = GaussLhood2(Valency, g1, g2, A11, A12, A22, Xindex, Flux1, Flux2);
    }

// Poisson
    else
    {
        CALL( PoissInsert2(MassInf, Rand, cool, q, ProbON,
                           Mock, Data, Acc, Counts, Foot,
                           Valency, nbits, ibits, zbits, nbitx, ibitx, zbitx,
                           Xindex, Flux1, Flux2) )
        *dL = PoissLhood2(Mock, Data, Acc, Foot,
                           Valency, nbits, ibits, zbits, nbitx, ibitx, zbitx,
                           Xindex, Flux1, Flux2);
    }

// Modify mock data
    InsBits(Mock, Flux1, nbits, ibits, zbits, Valency);
    InsBits(Mock, Flux2, nbitx, ibitx, zbitx, Valency);
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FluxDelete1
//
// Purpose:   Delete 1 old atom with fluxes from Object, keeping it up-to-date,
//            and set d(loglikelihood(fluxes)).
//
//            dL := d logL(...)
//
//            The old atom has been placed after the last location of
//            the updated Atoms list, at Object->Cubes[ Object->Natoms ].
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 12 Aug 2003, 18 Aug 2003
//-----------------------------------------------------------------------------
static
int FluxDelete1(        //   O  0, or -ve return code
double*    dL,          //   O  d(loglikelihood(fluxes))
CommonStr* Common,      // I    general information
ObjectStr* Object)      // I O  updated sample object
{
    int       MassInf = Common->MassInf;    // MassInf switches
    int       Valency = Common->Valency;    // # fluxes per atom
    int       Ndim    = Common->Ndim;       // dimension
    double*   Data    = Common->Data;       // Data                     [Ndata]
    double*   Acc     = Common->Acc;        // Accuracies               [Ndata]
    double*   Mock    = Object->Mock;       // Mock data                [Ndata]
    double*   g1      = Object->g1;         // grad chisquared        [Valency]
    double*   A11     = Object->A11;        // grad grad chisquared   [Valency]
    int*      nbits   = Object->nbits;      // # fragments            [Valency]
    int*      ibits   = Object->ibits;      // Fragment positions     [<=Ndata]
    double*   zbits   = Object->zbits;      // Fragment quantities    [<=Ndata]
    int*      flags   = Object->flags;      // flags for overlap        [Ndata]
    double*   Cube;                        // Position of new atom      [Ndim]
    double*   Flux;                         // Flux of new atom       [Valency]
    int       CALLvalue = 0;

    *dL = 0.0;                              // default output
// Footprint
    Cube = Object->Cubes[Object->Natoms];
    Flux = Cube + Ndim;
    CALL( AtomBits(Cube, Common, ibits, zbits, nbits, flags) )

// Gaussian
    if( MassInf < 100 )
    {
// Construct quadratic chisquared
        SetGrad1(Mock, Data, Acc, Valency, nbits, ibits, zbits, g1, A11);
// Modify Lhood
        *dL = GaussLhood1(Valency, Flux, g1, A11, -1);
    }

// Poisson
    else
        *dL = PoissLhood1(Mock, Data, Acc, Valency, nbits, ibits, zbits,
                          Flux, -1);

// Modify mock data
    DelBits(Mock, Flux, nbits, ibits, zbits, Valency);
Exit:
    return CALLvalue;
}

//=============================================================================
//             MassInf application procedures for user Footprints
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  AtomBits
//
// Purpose:   Wrapper for user's UserFoot, offering some protection against
//            user writing unusable fragments:
//                nbits[ Valency ]       elements obey nbits[.] >= 0
//                                                and  SUM(nbits[.]) <= Ndata
//                ibits[ SUM(nbits[.]) ] elements obey 0 <= ibits[.] < Ndata
//                zbits[ SUM(nbits[.]) ]
//
// History:   JS          21 Mar 2002, 3 Oct 2002, 12 Aug 2003
//-----------------------------------------------------------------------------
static
int AtomBits(         //   O  +ve = OK
                      //       0  = DO NOT USE this position
                      //      -ve = ERROR
double*    Cube,      // I    Atom position
CommonStr* Common,    // I    Definition of footprint
int*       ibits,     //   O  Fragment coords on data, serially per valency
double*    zbits,     //   O  Fragment quantities,     serially per valency
int*       nbits,     //   O  # fragments, for each valency
int*       flags)     //(I O) squash any un-necessary overlap from user
{
    int  Valency = Common->Valency;
    int  Ndata   = Common->Ndata;
    int  i, j, k, n, v;
    int  CALLvalue = 0;

// Footprint from user
    CALL( UserFoot(Cube, Common, ibits, zbits, nbits) )

    i = j = k = n = 0;
    for( v = 0; v < Valency; v++ )
    {
// Check for overwriting
        if( nbits[v] < 0 )
            return E_MASSINF_NBITS;
        n += nbits[v];
        if( n > Ndata )
            return E_MASSINF_NBITS;
// Ignore any zeros and overlay any overlaps
        for( ; i < n; i++ )
            if( zbits[i] != 0.0 )
            {
                if( ibits[i] < 0 || ibits[i] >= Ndata )
                    return E_MASSINF_RANGE;
                if( flags[ibits[i]] < 0 )
                {
                    flags[ibits[i]] = j;
                    ibits[j] = ibits[i];
                    zbits[j] = zbits[i];
                    j++;
                }
                else
                    zbits[flags[ibits[i]]] += zbits[i];
            }
        nbits[v] = j - k;
        for( ; k < j; k++ )
            flags[ibits[k]] = -1;
     }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Overlap1
//
// Purpose:   Ensure (before inserting) that atom valencies do not overlap.
//
//            Dot products of individual footprints must be diagonal, else
//            error return.
//                    --------------------------
//                   | A11[0]  .     .     .    |
//                   |                          |
//                   |   .   A11[1]  .     .    |
//                   |                          |
//                   |   .     .   A11[2]  .    |
//                   |                          |
//                   |   .     .     .   A11[3] |
//                    --------------------------
//                           --> Valency
//
// History:   JS           2 Jan 2002, 21 Mar 2002
//-----------------------------------------------------------------------------
static
int Overlap1(      //   O  0, or -ve user error
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
int*     flags,    // I O  Valency identifiers (-1)  [Ndata]
int      Valency)  // I    # fluxes per atom
{
    int      i;        // bits counter
    int      k;        // data counter
    int      v;        // valency counter
    int      n;        // bits limiter

// Check Atom valencies do not overlap by setting valency flags
    i = n = 0;
    for( v = 0; v < Valency; v++ )
    {
        n += nbits[v];
        for( ; i < n; i++ )
        {
            k = ibits[i];
            if( flags[k] == -1 )
                flags[k] = v;
            else if( flags[k] != v )
                return E_MASSINF_OVERLAP;
        }
    }

// Reset flags = -1 (off)
    for( n--; n >= 0; n-- )
        flags[ibits[n]] = -1;
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Overlap2
//
// Purpose:   Ensure (before inserting) that valencies of two atoms do not
//            wrongly overlap.
//
//            Dot products of individual footprints must be diagonal, and
//            each valency of an atom can overlap with at most one valency
//            of another, else error return.
//
//                   Atom1                      Atom2
//         -----------------------------------------------------
//        | A11[0]  .     .     .    |   .     .   A12[0]  .    | Xindex[0] = 2
//        |                          |                          |
//        |   .   A11[1]  .     .    |   .     .     .     .    | Xindex[1] off
//  Atom1 |                          |                          |
//        |   .     .   A11[2]  .    |   .     .     .     .    | Xindex[2] off
//        |                          |                          |
//        |   .     .     .   A11[3] |   .   A12[3]  .     .    | Xindex[3] = 1
//        |--------------------------|--------------------------|
//        |   .     .     .     .    | A22[0]  .     .     .    |
//        |                          |                          |
//        |   .     .     .   A12[3] |   .   A22[1]  .     .    |
//  Atom2 |                          |                          |
//        | A12[0]  .     .     .    |   .     .   A22[2]  .    |
//        |                          |                          |
//        |   .     .     .     .    |   .     .     .   A22[3] |
//         -----------------------------------------------------
//                       --> Valency                 -->Valency
//
//            Cross-terms between the two atoms (at most 1 per row or column)
//            are indexed by Xindex, as in this example with Valency = 4.
//
// History:   JS           2 Jan 2002
//-----------------------------------------------------------------------------
static
int Overlap2(      //   O  0, or -ve user error
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
int*     nbitx,    // I    # fragments             [Valency]
int*     ibitx,    // I    Fragment positions      [<=Ndata]
int*     flags,    // I O  Valency identifiers (-1)  [Ndata]
int*     Xindex,   //   O  index to cross-terms    [Valency]
int      Valency)  // I    # fluxes per atom
{
    int      i;        // bits counter
    int      k;        // data counter
    int      v1;       // valency counter
    int      v2;       // valency counter
    int      n1;       // bits limiter
    int      n2;       // bits limiter
    int      check;    // crossterm checker

// Check Atom1 valencies do not overlap by setting valency flags
    i = n1 = 0;
    for( v1 = 0; v1 < Valency; v1++ )
    {
        n1 += nbits[v1];
        for( ; i < n1; i++ )
        {
            k = ibits[i];
            if( flags[k] == -1 )
                flags[k] = v1;
            else if( flags[k] != v1 )
                return E_MASSINF_OVERLAP;
        }
    }

// Check no more than one Atom1 overlap with each Atom2 valency
    i = n2 = 0;
    for( v2 = 0; v2 < Valency; v2++ )
    {
        check = -1;
        n2 += nbitx[v2];
        for( ; i < n2; i++ )
        {
            k = ibitx[i];
            if( flags[k] != -1 )
            {
                if( check == -1 )
                    check = flags[k];
                else if( check != flags[k] )
                    return E_MASSINF_OVERLAP;
            }
        }
    }

// Reset flags = -1 (off)
    for( n1--; n1 >= 0; n1-- )
        flags[ibits[n1]] = -1;

// Check Atom2 valencies do not overlap
    i = n2 = 0;
    for( v2 = 0; v2 < Valency; v2++ )
    {
        n2 += nbitx[v2];
        for( ; i < n2; i++ )
        {
            k = ibitx[i];
            if( flags[k] == -1 )
                flags[k] = v2;
            else if( flags[k] != v2 )
                return E_MASSINF_OVERLAP;
        }
    }

// Check no more than one Atom2 overlap with each Atom1 valency
// Set Xindex[.] = Atom2 valency overlapping with Atom1 valency (-1 if none)
    i = n1 = 0;
    for( v1 = 0; v1 < Valency; v1++ )
    {
        check = -1;
        n1 += nbits[v1];
        for( ; i < n1; i++ )
        {
            k = ibits[i];
            if( flags[k] != -1 )
            {
                if( check == -1 )
                    check = flags[k];
                else if( check != flags[k] )
                    return E_MASSINF_OVERLAP;
            }
        }
        Xindex[v1] = check;
    }

// Reset flags = -1 (off)
    for( n2--; n2 >= 0; n2-- )
        flags[ibitx[n2]] = -1;
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SetIndex
//
// Purpose:   Set Xindex cross-terms between valencies.
//
// History:   JS          18 Aug 2003 from SetGrad2
//-----------------------------------------------------------------------------
static
void SetIndex(
int      Valency, // I    # fluxes per atom
int*     nbits,   // I    # fragments             [Valency]
int*     ibits,   // I    Fragment positions      [<=Ndata]
int*     nbitx,   // I    # fragments             [Valency]
int*     ibitx,   // I    Fragment positions      [<=Ndata]
int*     flags,   // I O  Valency identifiers (=-1) [Ndata]
int*     Xindex)  //   O  index to cross-terms    [Valency]
{
    int      i;   // bits counter
    int      k;   // data counter
    int      v1;  // valency counter
    int      v2;  // valency counter
    int      m1;  // bits limiter
    int      n1;  // bits limiter
    int      n2;  // bits limiter

// Set valency flags
    i = n2 = 0;
    for( v2 = 0; v2 < Valency; v2++ )
    {
        n2 += nbitx[v2];
        for( ; i < n2; i++ )
            flags[ibitx[i]] = v2;
    }

// Set Xindex[.] = Atom2 valency overlapping with Atom1 valency (-1 if none)
    m1 = n1 = 0;
    for( v1 = 0; v1 < Valency; v1++ )
    {
        Xindex[v1] = -1;
        n1 += nbits[v1];
        for( i = m1; i < n1; i++ )
        {
            k = ibits[i];
            if( flags[k] >= 0 )
            {
                Xindex[v1] = flags[k];
                break;
            }
        }
        m1 += nbits[v1];
    }

// Reset flags = -1 (off)
    n2 = 0;
    for( v2 = 0; v2 < Valency; v2++ )
        n2 += nbitx[v2];
    for( i = 0; i < n2; i++ )
        flags[ibitx[i]] = -1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  InsBits
//
// Purpose:   Insert atom fluxes into mock data
//
// History:   JS           2 Jan 2002
//-----------------------------------------------------------------------------
static
void InsBits(
double*  Mock,     // I O  Mock data
double*  Flux,     // I    fluxes                 [Valency]
int*     nbit,     // I    # fragments            [Valency]
int*     ibit,     // I    fragment positions     [<=Ndata]
double*  zbit,     // I    fragment quantities    [<=Ndata]
int      Valency)  // I    # fluxes per atom
{
    int  i, v, n;

    i = n = 0;
    for( v = 0; v < Valency; v++ )
    {
        n += nbit[v];
        for( ; i < n; i++ )
            Mock[ibit[i]] += Flux[v] * zbit[i];
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  DelBits
//
// Purpose:   Delete atom fluxes from mock data
//
// History:   JS           2 Jan 2002
//-----------------------------------------------------------------------------
static
void DelBits(
double*  Mock,     // I O  Mock data
double*  Flux,     // I    fluxes                 [Valency]
int*     nbit,     // I    # fragments            [Valency]
int*     ibit,     // I    fragment positions     [<=Ndata]
double*  zbit,     // I    fragment quantities    [<=Ndata]
int      Valency)  // I    # fluxes per atom
{
    int  i, v, n;

    i = n = 0;
    for( v = 0; v < Valency; v++ )
    {
        n += nbit[v];
        for( ; i < n; i++ )
            Mock[ibit[i]] -= Flux[v] * zbit[i];
    }
}

//=============================================================================
//          MassInf application procedures for quadratic chisquared
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SetGrad1
//
// Purpose:   Set quadratic coeffs for increment of logLikelihood from 1 atom,
//            allowing for multiple Valency.
//
//               d(Lhood)  =  - g1.Flux  -  Flux.A11.Flux / 2
//
// History:   JS           2 Jan 2002
//                        14 Aug 2003 Rely on no overlay in received bits
//-----------------------------------------------------------------------------
static
void SetGrad1(
double*  Mock,    // I    Mock data                 [Ndata]
double*  Data,    // I    Data                      [Ndata]
double*  Acc,     // I    Accuracies                [Ndata]
int      Valency, // I    # fluxes per atom
int*     nbits,   // I    # fragments             [Valency]
int*     ibits,   // I    Fragment positions      [<=Ndata]
double*  zbits,   // I    Fragment quantities     [<=Ndata]
double*  g1,      //   O  grad chisquared         [Valency]
double*  A11)     //   O  grad grad chisquared    [Valency]
{
    int      i;   // bits counter
    int      k;   // data counter
    int      v;   // valency counter
    int      m;   // bits limiter
    int      n;   // bits limiter
    double   z;   // local footprint

    m = n = 0;
    for( v = 0; v < Valency; v++ )
    {
        g1[v] = A11[v] = 0.0;
        n += nbits[v];
        for( i = m; i < n; i++ )
        {
            k = ibits[i];
            z = zbits[i] * Acc[k];
            g1[v] += z * Acc[k] * (Mock[k] - Data[k]);
            A11[v] += z * z;
        }
        m += nbits[v];
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SetGrad2
//
// Purpose:   Set quadratic coeffs for increment of logLikelihood from 2 atoms,
//            allowing for multiple Valency.
//
//    d(Lhood) = - g1.Flux1 - g2.Flux2
//               - (Flux1.A11.Flux1 + 2 Flux1.A12.Flux2 + Flux2.A22.Flux2) / 2
//
// History:   JS           2 Jan 2002
//                        18 Aug 2003 Rely on no overlay in received bits
//-----------------------------------------------------------------------------
static
void SetGrad2(
double*  Mock,    // I    Mock data                 [Ndata]
double*  Data,    // I    Data                      [Ndata]
double*  Acc,     // I    Accuracies                [Ndata]
double*  Foot,    // I O  Footprint workspace (=0)  [Ndata]
int      Valency, // I    # fluxes per atom
int*     nbits,   // I    # fragments             [Valency]
int*     ibits,   // I    Fragment positions      [<=Ndata]
double*  zbits,   // I    Fragment quantities     [<=Ndata]
int*     nbitx,   // I    # fragments             [Valency]
int*     ibitx,   // I    Fragment positions      [<=Ndata]
double*  zbitx,   // I    Fragment quantities     [<=Ndata]
double*  g1,      //   O  grad chisquared         [Valency]
double*  g2,      //   O  grad chisquared         [Valency]
double*  A11,     //   O  grad grad chisquared    [Valency]
double*  A12,     //   O  grad grad chisquared    [Valency]
double*  A22,     //   O  grad grad chisquared    [Valency]
int*     Xindex)  // I    index to cross-terms    [Valency]
{
    int      i;   // bits counter
    int      k;   // data counter
    int      v1;  // valency counter
    int      v2;  // valency counter
    int      m1;  // bits limiter
    int      n1;  // bits limiter
    int      n2;  // bits limiter
    double   z;   // local flux

// For each Atom1 valency
    m1 = n1 = 0;
    for( v1 = 0; v1 < Valency; v1++ )
    {
        g1[v1] = g2[v1] = A11[v1] = A12[v1] = A22[v1] = 0.0;
// Footprint
        n1 += nbits[v1];
        for( i = m1; i < n1; i++ )
            Foot[ibits[i]] = zbits[i];
// Index to cross term (if any)
        if( Xindex[v1] >= 0 )
        {
            i = 0;
            for( v2 = 0; v2 < Xindex[v1]; v2++ )
                i += nbitx[v2];
            n2 = i + nbitx[v2];
            for( ; i < n2; i++ )
            {
                k = ibitx[i];
                A12[v1] += zbitx[i] * Acc[k] * Acc[k] * Foot[k];
            }
        }
// Diagonal term
        for( i = m1; i < n1; i++ )
        {
            k = ibits[i];
            z = Foot[k] * Acc[k] * Acc[k];
            g1[v1]  += z * (Mock[k] - Data[k]);
            A11[v1] += z * Foot[k];
            Foot[k] = 0.0;
        }
        m1 += nbits[v1];
    }

// For each Atom2 valency, diagonal term
    i = n2 = 0;
    for( v2 = 0; v2 < Valency; v2++ )
    {
        n2 += nbitx[v2];
        for( ; i < n2; i++ )
        {
            k = ibitx[i];
            z = zbitx[i] * Acc[k] * Acc[k];
            g2[v2]  += z * (Mock[k] - Data[k]);
            A22[v2] += z * zbitx[i];
        }
    }
}

//=============================================================================
//          MassInf application procedures for quadratic probabilities
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussLhood
//                                       -1  -chisquared/2
// Purpose:   Set Lhood = log(L),   L = Z   e
//
// History:   JS          18 Aug 2003
//-----------------------------------------------------------------------------
static
double GaussLhood( //   O  log(L)
int      Ndata,    // I    # data
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc)      // I    Accuracies                [Ndata]
{
static const double sqrt2pi = 2.50662827463100050240;
    double  C     = 0.0;
    double  Lnorm = 0.0;
    double  temp;
    int     k;

    for( k = 0; k < Ndata; k++ )
    {
        if( Acc[k] > 0.0 )
        {
            temp   = Acc[k];
            Lnorm += log(temp / sqrt2pi);
            temp  *= Mock[k] - Data[k];
            C     += temp * temp;
        }
    }
    return Lnorm - C / 2.0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussTry1
//
//              (1/cool) log INTEGRAL dPrior (Lhood / L0)^cool
//
//            for entire atom with Gaussian likelihood
//
// History:   JS   19 Jul 2003, 14 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
double GaussTry1(      //   O  DELTA logL
int      MassInf,      // I    MassInf prior
double   cool,         // I    Annealing coefficient
double   q,            // I    Flux unit
double   ProbON,       // I    Pr(individual flux != 0)
int      Valency,      // I    # fluxes per atom
double*  g1,           // I    grad chisquared         [Valency]
double*  A11)          // I    grad grad chisquared    [Valency]
{
    double  dL = 0.0;
    int     v;

    for( v = 0; v < Valency; v++ )
        dL += Gauss1Marginal(MassInf, cool, q, ProbON, g1[v], A11[v]);
    return dL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussTry2
//
//              (1/cool) log INTEGRAL dPrior (Lhood / L0)^cool
//
//            for two entire atoms with Gaussian likelihood
//
// History:   JS   19 Jul 2003, 14 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
void GaussTry2(
int      MassInf,      // I    MassInf prior
double   cool,         // I    Annealing coefficient
double   q,            // I    Flux unit
double   ProbON,       // I    Pr(individual flux != 0)
int      Valency,      // I    # fluxes per atom
double*  g1,           // I    grad chisquared         [Valency]
double*  g2,           // I    grad chisquared         [Valency]
double*  A11,          // I    grad grad chisquared    [Valency]
double*  A12,          // I    grad grad chisquared    [Valency]
double*  A22,          // I    grad grad chisquared    [Valency]
int*     Xindex,       // I    index to cross-terms    [Valency]
double*  dLtry1,       //   O  DELTA logL for first atom
double*  dLtry2)       //   O  DELTA logL for both atoms
{
    int     v1, v2;
    *dLtry1 = *dLtry2 = 0.0;

// Trial Lhood increment
    for( v1 = 0; v1 < Valency; v1++ )
        *dLtry1 += Gauss1Marginal(MassInf, cool, q, ProbON, g1[v1], A11[v1]);

// Trial Lhood increment
    for( v1 = 0; v1 < Valency; v1++ )
    {
        v2 = Xindex[v1];
        if( v2 >= 0 )                     // pick up Atom1/Atom2 cross term...
            *dLtry2 += Gauss2Marginal(MassInf, cool, q, ProbON,
                                 g1[v1], g2[v2], A11[v1], A12[v1], A22[v2]);
        else                              // ... or pick up lone Atom1
            *dLtry2 += Gauss1Marginal(MassInf, cool,q,ProbON, g1[v1],A11[v1]);
    }
    for( v2 = 0; v2 < Valency; v2++ )
    {
        for( v1 = 0; v1 < Valency; v1++ )
           if( Xindex[v1] == v2 )
                 break;
        if( v1 == Valency )               // remaining lone Atom2
            *dLtry2 += Gauss1Marginal(MassInf, cool,q,ProbON, g2[v2],A22[v2]);
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussInsert1
//
//            Sample fluxes from    Prior * Lhood^cool
//            for one entire atom with Gaussian likelihood
//
// History:   JS   19 Jul 2003, 14 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
void GaussInsert1(
int      MassInf,    // I    MassInf prior
Rand_t   Rand,       // I O  random generator state
double   cool,       // I    Annealing coefficient
double   q,          // I    Flux unit
double   ProbON,     // I    Pr(individual flux != 0)
int      Valency,    // I    # fluxes per atom
double*  g1,         // I    grad chisquared         [Valency]
double*  A11,        // I    grad grad chisquared    [Valency]
double*  Flux)       //   O  fluxes                  [Valency]
{
    int  v;
    for( v = 0; v < Valency; v++ )
        Flux[v] = Gauss1Sample(MassInf, Rand, cool, q, ProbON, g1[v], A11[v]);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussInsert2
//
//            Sample fluxes from    Prior * Lhood^cool
//            for two entire atoms with Gaussian likelihood
//
// History:   JS   19 Jul 2003, 14 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
void GaussInsert2(
int      MassInf,    // I    MassInf prior
Rand_t   Rand,       // I O  random generator state
double   cool,       // I    Annealing coefficient
double   q,          // I    Flux unit
double   ProbON,     // I    Pr(individual flux != 0)
int      Valency,    // I    # fluxes per atom
double*  g1,         // I    grad chisquared         [Valency]
double*  g2,         // I    grad chisquared         [Valency]
double*  A11,        // I    grad grad chisquared    [Valency]
double*  A12,        // I    grad grad chisquared    [Valency]
double*  A22,        // I    grad grad chisquared    [Valency]
int*     Xindex,     // I    index to cross-terms    [Valency]
double*  Flux1,      //   O  fluxes of atom1         [Valency]
double*  Flux2)      //   O  fluxes of atom2         [Valency]
{
    int     v1, v2;

    for( v1 = 0; v1 < Valency; v1++ )
    {
        v2 = Xindex[v1];
        if( v2 >= 0 )
            Gauss2Sample(MassInf, Rand, cool, q, ProbON,
                         g1[v1], g2[v2], A11[v1], A12[v1], A22[v2],
                         &Flux1[v1], &Flux2[v2]);
        else
            Flux1[v1] = Gauss1Sample(MassInf, Rand, cool, q, ProbON,
                                     g1[v1], A11[v1]);
    }
    for( v2 = 0; v2 < Valency; v2++ )
    {
        for( v1 = 0; v1 < Valency; v1++ )
            if( Xindex[v1] == v2 )
                break;
        if( v1 == Valency )
            Flux2[v2] = Gauss1Sample(MassInf, Rand, cool, q, ProbON,
                                     g2[v2], A22[v2]);
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussLhood1
//
//            DELTA( log Lhood )  for inserting or deleting
//            one entire atom with Gaussian likelihood
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
double GaussLhood1(
int      Valency,   // I    # fluxes per atom
double*  Flux ,     // I    fluxes                  [Valency]
double*  g1,        // I    grad chisquared         [Valency]
double*  A11,       // I    grad grad chisquared    [Valency]
int      flag)      // I    +1 = insert, -1 = delete
{
    double  dL = 0.0;
    int     v;
    if( flag > 0 )
        for( v = 0; v < Valency; v++ )
            dL -= Flux[v] * (Flux[v] * A11[v] / 2.0 + g1[v]);
    else
        for( v = 0; v < Valency; v++ )
            dL -= Flux[v] * (Flux[v] * A11[v] / 2.0 - g1[v]);
    return dL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  GaussLhood2
//
//            DELTA( log Lhood )
//            for inserting two entire atoms with Gaussian likelihood
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
double GaussLhood2(
int      Valency,   // I    # fluxes per atom
double*  g1,        // I    grad chisquared         [Valency]
double*  g2,        // I    grad chisquared         [Valency]
double*  A11,       // I    grad grad chisquared    [Valency]
double*  A12,       // I    grad grad chisquared    [Valency]
double*  A22,       // I    grad grad chisquared    [Valency]
int*     Xindex,    // I    index to cross-terms    [Valency]
double*  Flux1,     // I    fluxes of atom1         [Valency]
double*  Flux2)     // I    fluxes of atom2         [Valency]
{
    double  dL = 0.0;
    int     v1, v2;

    for( v1 = 0; v1 < Valency; v1++ )
    {
        v2 = Xindex[v1];
        if( v2 >= 0 )
            dL -= g1[v1] * Flux1[v1] + g2[v2] * Flux2[v2]
                 + A11[v1] * Flux1[v1] * Flux1[v1] / 2.0
                 + A12[v1] * Flux1[v1] * Flux2[v2]
                 + A22[v2] * Flux2[v2] * Flux2[v2] / 2.0;
        else
            dL -= g1[v1] * Flux1[v1] + A11[v1] * Flux1[v1] * Flux1[v1] / 2.0;
    }
    for( v2 = 0; v2 < Valency; v2++ )
    {
        for( v1 = 0; v1 < Valency; v1++ )
            if( Xindex[v1] == v2 )
                break;
        if( v1 == Valency )
            dL -= g2[v2] * Flux2[v2] + A22[v2] * Flux2[v2] * Flux2[v2] / 2.0;
    }
    return dL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Gauss1Marginal
//                                               -cool(g*z + A*z*z/2)
//               (1/cool) log INTEGRAL dz Pr(z) e
//
// History:   JS   19 Jul 2003, 18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
double Gauss1Marginal( //   O  value, calculated OK even if cool=0
int      MassInf,      // I    MassInf prior (0,1,2,3)
double   cool,         // I    Annealing coefficient
double   q,            // I    Flux unit
double   ProbON,       // I    Pr(individual flux != 0)
double   g,            // I    gradient
double   A)            // I    curvature
{
    double  p;
    double  ProbOFF;

    p = gauss1marginal(MassInf, cool, q, g, A);
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        if( cool * cool * p * p < DBL_EPSILON )
            p *= ProbON;
        else
            p = PLUS(log(ProbOFF), log(ProbON) + cool * p) / cool;
    }
    return  p;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Gauss2Marginal
//                                                       -cool(g.x + x.A.x/2)
//                (1/cool) log INTEGRAL dxdy Pr(x)Pr(y) e
//
//            where  g.x = g1*x + g2*y,  x.A.x = A11*x*x + 2*A12*x*y + A22*y*y
//
// History:   JS   19 Jul 2003, 18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
double Gauss2Marginal( //   O  value, calculated OK even if cool=0
int      MassInf,      // I    MassInf prior (0,1,2,3)
double   cool,         // I    Annealing coefficient
double   q,            // I    Flux unit
double   ProbON,       // I    Pr(individual flux != 0)
double   g1,           // I    gradient
double   g2,           // I    gradient
double   A11,          // I    curvature
double   A12,          // I    curvature
double   A22)          // I    curvature
{
    double  p1, p2, p12;
    double  ProbOFF;

    p12 = gauss2marginal(MassInf, cool, q, g1, g2, A11, A12, A22);
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        p1  = gauss1marginal(MassInf, cool, q, g1, A11);
        p2 = gauss1marginal(MassInf, cool, q, g2, A22);
        if( cool * cool * (p1*p1 + p2*p2 + p12*p12) < DBL_EPSILON )
        {
            p12 = ProbON * (ProbOFF * (p1 + p2) + ProbON * p12);
            p1 *= ProbON;
        }
        else
        {
            p1 = PLUS(log(ProbOFF), log(ProbON) + cool*p1);
            p2 = PLUS(log(ProbOFF) + cool*p2, log(ProbON) + cool*p12);
            p12 = PLUS(log(ProbOFF) + p1, log(ProbON) + p2) / cool;
        }
    }
    return p12;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Gauss1Sample
//                                                -cool(g*z + A*z*z/2)
//            Sample single flux z from    Pr(z) e
//
// History:   JS   19 Jul 2003, 18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
double Gauss1Sample( //   O  flux z
int      MassInf,    // I    MassInf prior (0,1,2,3)
Rand_t   Rand,       // I O  random generator state
double   cool,       // I    cooling parameter
double   q,          // I    Flux unit
double   ProbON,     // I    Pr(individual flux != 0)
double   g,          // I    gradient
double   A)          // I    curvature
{
    double  p, ProbOFF;
    int     k = 1;

    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        p = log(ProbON) + cool * gauss1marginal(MassInf, cool, q, g, A);
        if( log(Randouble(Rand)) > p - PLUS(log(ProbOFF), p) )
            k = 0;
    }
    return  k ? gauss1sample(MassInf, Rand, cool, q, g, A) : 0.0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Gauss2Sample
//                                                         -cool(g.x + x.A.x/2)
//            Sample single fluxes (x,y) from  Pr(x)Pr(y) e
//            where  g.x = g1*x + g2*y,  x.A.x = A11*x*x + 2*A12*x*y + A22*y*y
//
// History:   JS   19 Jul 2003, 18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
void Gauss2Sample(
int      MassInf,    // I    MassInf prior (0,1,2,3)
Rand_t   Rand,       // I O  random generator state
double   cool,       // I    cooling parameter
double   q,          // I    Flux unit
double   ProbON,     // I    Pr(individual flux != 0)
double   g1,         // I    gradient
double   g2,         // I    gradient
double   A11,        // I    curvature
double   A12,        // I    curvature
double   A22,        // I    curvature
double*  Flux1,      //   O  flux
double*  Flux2)      //   O  flux
{
    double  p, p1, p2, p12, a, b, c, r, ProbOFF;
    int     k = 3;

    *Flux1 = *Flux2 = 0.0;
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        p   = log(ProbOFF) + log(ProbOFF);
        p1  = log(ProbOFF) + log(ProbON)
             + cool * gauss1marginal(MassInf, cool, q, g1, A11);
        p2  = log(ProbOFF) + log(ProbON)
             + cool * gauss1marginal(MassInf, cool, q, g2, A22);
        p12 = log(ProbON) + log(ProbON)
             + cool * gauss2marginal(MassInf, cool, q, g1, g2, A11, A12, A22);
        a = PLUS(p, p1);
        b = PLUS(p2, p12);
        c = PLUS(a, b);
        p1  = exp(p1 - c);
        p2  = exp(p2 - c);
        p12 = exp(p12 - c);
        r = Randouble(Rand);
        if( r < p12 )
            k = 3;
        else if( r < p12 + p2 )
            k = 2;
        else if( r < p12 + p2 + p1 )
            k = 1;
        else
            k = 0;
    }
    if( k == 1 )
        *Flux1 = gauss1sample(MassInf, Rand, cool, q, g1, A11);
    if( k == 2 )
        *Flux2 = gauss1sample(MassInf, Rand, cool, q, g2, A22);
    if( k == 3 )
        gauss2sample(MassInf, Rand, cool, q, g1,g2, A11,A12,A22, Flux1, Flux2);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  gauss1marginal
//                                              -cool(g*z + A*z*z/2)
//              (1/cool) log INTEGRAL dz Pr(z) e
//
//            where z is non-zero flux
//
// History:   JS    2 Jan 2002
//                 19 Jul 2003   "switch" recoding
//-----------------------------------------------------------------------------
static
double gauss1marginal( //   O  value, calculated OK even if cool=0
int     MassInf,       // I    Prior (0=monkeys,1=pos,2=posneg,3=gauss)
double  cool,          // I    Annealing coefficient
double  q,             // I    Flux unit
double  g,             // I    gradient
double  A)             // I    curvature
{
    double  t, pos, neg;
    g *= q;
    A *= q * q;

    switch( MassInf % 10 )
    {
    case 0:
         return -(g + A / 2.0);

    case 1:
        t = cool * (fabs(g) + A);
        if( t * t > DBL_EPSILON )
            return logerf(1.0 + cool * g, cool * A) / cool;
        else
            return -(g + A);

    case 2:
        t = cool * (cool * g * g + A);
        if( t * t > DBL_EPSILON )
        {
            pos = logerf(1.0 + cool * g, cool * A);
            neg = logerf(1.0 - cool * g, cool * A);
            return (log(0.5) + PLUS(pos, neg)) / cool;
        }
        else
            return  cool * g * g - A;

    case 3:
        t = cool * A;
        if( t * t > DBL_EPSILON )
            t = log(1.0 + t) / cool;
        else
            t = A;
        return 0.5 * (cool * g * g / (1.0 + cool * A) - t);
    }
    return  0.0;  // calling error
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  gauss2marginal
//                                                     -cool(g.x + x.A.x/2)
//                (1/cool) log INTEGRAL dxdy P(x)P(y) e
//
//            where  g.x = g1*x + g2*y,  x.A.x = A11*x*x + 2*A12*x*y + A22*y*y
//            and  x,y  are non-zero fluxes
//
// History:   JS    2 Jan 2002
//                 19 Jul 2003   "switch" recoding
//-----------------------------------------------------------------------------
static
double gauss2marginal( //   O  value, calculated OK even if cool=0
int     MassInf,       // I    Prior (0=monkeys,1=pos,2=posneg,3=gauss)
double  cool,          // I    Annealing coefficient
double  q,             // I    Flux unit
double  g1,            // I    gradient
double  g2,            // I    gradient
double  A11,           // I    curvature
double  A12,           // I    curvature
double  A22)           // I    curvature
{
    double  t, pospos, posneg, negpos, negneg, pos, neg, det, gBg;
    g1 *= q;
    g2 *= q;
    A11 *= q * q;
    A12 *= q * q;
    A22 *= q * q;
    t = sqrt(A11 * A22);
    if( A12 > t )
        A12 = t;
    if( A12 < -t )
        A12 = -t;

    switch( MassInf % 10 )
    {
    case 0:
        return -(g1 + g2 + 0.5 * A11 + A12 + 0.5 * A22);

    case 1:
        t = cool * (fabs(g1) + fabs(g2) + A11 + A22);
        if( t * t > DBL_EPSILON )
            return logerf2(1.0 + cool * g1, 1.0 + cool * g2,
                           cool * A11, cool * A12, cool * A22) / cool;
        else
            return -(g1 + g2 + A11 + A12 + A22);

    case 2:
        t = cool * (cool * (g1 * g1 + g2 * g2) + A11 + A22);
        if( t * t > DBL_EPSILON )
        {
            pospos = logerf2(1.0 + cool * g1, 1.0 + cool * g2,
                             cool * A11, cool * A12, cool * A22);
            posneg = logerf2(1.0 + cool * g1, 1.0 - cool * g2,
                             cool * A11, cool * A12, cool * A22);
            negpos = logerf2(1.0 - cool * g1, 1.0 + cool * g2,
                             cool * A11, cool * A12, cool * A22);
            negneg = logerf2(1.0 - cool * g1, 1.0 - cool * g2,
                             cool * A11, cool * A12, cool * A22);
            pos = PLUS(pospos, posneg);
            neg = PLUS(negpos, negneg);
            return (log(0.25) + PLUS(pos, neg)) / cool;
        }
        else
            return cool * (g1 * g1 + g2 * g2) - (A11 + A12 + A22);

    case 3:
        t = cool * (A11 + A22);
        det = A11 * A22 - A12 * A12;
        if( det < 0.0 )
            det = 0.0;
        det = 1.0 + t + cool * cool * det;
        gBg = (g1 * g1 * (1.0 + cool * A22) - 2.0 * g1 * g2 * cool * A12
              + g2 * g2 * (1.0 + cool * A11)) / det;
        if( t * t > DBL_EPSILON )
            t = log(det) / cool;
        else
            t = A11 + A22;
        return 0.5 * (cool * gBg - t);
    }
    return  0.0;  // calling error
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  gauss1sample
//                                                 -cool(g*z + A*z*z/2)
//            Sample non-zero flux z from    P(z) e
//
// History:   JS    2 Jan 2002
//                 17 Dec 2002   MassInf=3 debugged
//                 19 Jul 2003   "switch" recoding
//                 13 Aug 2003   q into output
//-----------------------------------------------------------------------------
static
double gauss1sample( //   O  flux z
int      MassInf,    // I    Prior (0=monkeys,1=pos,2=posneg,3=gauss)
Rand_t   Rand,       // I O  random generator state
double   cool,       // I    cooling parameter
double   q,          // I    Flux unit
double   g,          // I    gradient
double   A)          // I    curvature
{
    g *= cool * q;
    A *= cool * q * q;
    switch( MassInf % 10 )
    {
    case 0:
        return  q;

    case 1:
        return  q * Ran1pos(Rand, g + 1.0, A);

    case 2:
        return  q * Ran1posneg(Rand, g, 1.0, A);

    case 3:
        return  q * Rangauss(Rand) / sqrt(1.0 + A) - g / (1.0 + A);
    }
    return  0.0;  // calling error
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  gauss2sample
//                                                         -cool(g.x + x.A.x/2)
//            Sample non-zero fluxes (x,y) from  P(x)P(y) e
//
//            where  g.x = g1*x + g2*y,  x.A.x = A11*x*x + 2*A12*x*y + A22*y*y
//
// History:   JS    2 Jan 2002
//                 19 Jul 2003   "switch" recoding
//                 13 Aug 2003   q into output
//-----------------------------------------------------------------------------
static
void gauss2sample(
int      MassInf,    // I    Prior (0=monkeys,1=pos,2=posneg,3=gauss)
Rand_t   Rand,       // I O  random generator state
double   cool,       // I    cooling parameter
double   q,          // I    Flux unit
double   g1,         // I    gradient
double   g2,         // I    gradient
double   A11,        // I    curvature
double   A12,        // I    curvature
double   A22,        // I    curvature
double*  x,          //   O  flux
double*  y)          //   O  flux
{
    double  t;

    g1 *= cool * q;
    g2 *= cool * q;
    A11 *= cool * q * q;
    A12 *= cool * q * q;
    A22 *= cool * q * q;
    t = sqrt(A11 * A22);
    if( A12 > t )
        A12 = t;
    if( A12 < -t )
        A12 = -t;
    switch( MassInf % 10 )
    {
    case 0:
        *x = *y = 1.0;
        break;

    case 1:
        Ran2pos(Rand, g1 + 1.0, g2 + 1.0, A11, A12, A22, x, y);
        break;

    case 2:
        Ran2posneg(Rand, g1, g2, 1.0, 1.0, A11, A12, A22, x, y);
        break;

    case 3:
        Ran2gauss(Rand, g1, g2, 1.0 + A11, A12, 1.0 + A22, x, y);
        break;
    default:
        *x = *y = 0.0;   // calling error
    }
    *x *= q;
    *y *= q;
}

//=============================================================================
//          MassInf application procedures for Poisson probabilities
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissLhood
//
// Purpose:   Set Lhood = log(L) where
//
//                             -(F+B)     D+B
//                L = PRODUCT e      (F+B)    / (D+B)!
//
//            D = data, F = mock data, B = background > 0 wherever D+B > 0.
//            Maximum likelihood is F = D, as wanted.
//            For small F, L remains finite, so empty object is non-singular.
//            For small B, L becomes exp(-F) F^D / D!, standard Poisson.
//
// History:   JS          18 Aug 2003, 25 Oct 2003
//-----------------------------------------------------------------------------
static
double PoissLhood( //   O  log(L)
int      Ndata,    // I    # data
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc)      // I    Background                [Ndata]
{
    double  C = 0.0;
    double  D, F;
    int     k;

    for( k = 0; k < Ndata; k++ )
    {
        D = Data[k] + Acc[k];
        F = Mock[k] + Acc[k];
        if( F > 0.0 )               // else D known to be 0 also
            C += D * log(F) - F - logGamma(1.0 + D);
    }
    return C;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissTry1
//
//              (1/cool) log INTEGRAL dPrior (Lhood / L0)^cool
//
//            for entire atom with Poisson likelihood
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int PoissTry1(     //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
int      Valency,  // I    # fluxes per atom
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
double*  zbits,    // I    Fragment quantities     [<=Ndata]
double*  dLtry)    //   O  (1/cool) INT dprior (L/L0)^cool
{
    double  dL = 0.0;
    double  t;
    int     k, v;
    int     CALLvalue = 0;

    for( k = v = 0; v < Valency; k += nbits[v++] )
    {
        CALL( Poiss1Marginal(MassInf, cool, q, ProbON, Mock, Data, Acc, Counts,
                             nbits[v], &ibits[k], &zbits[k], &t) )
        dL += t;
    }
    *dLtry = dL;
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissTry2
//
//              (1/cool) log INTEGRAL dPrior (Lhood / L0)^cool
//
//            for two entire atoms with Poisson likelihood
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int PoissTry2(     //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      Valency,  // I    # fluxes per atom
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
double*  zbits,    // I    Fragment quantities     [<=Ndata]
int*     nbitx,    // I    # fragments             [Valency]
int*     ibitx,    // I    Fragment positions      [<=Ndata]
double*  zbitx,    // I    Fragment quantities     [<=Ndata]
int*     Xindex,   // I    index to cross terms    [Valency]
double*  dLtry1,   //   O  (1/cool) INT dprior[1] (L/L0)^cool
double*  dLtry2)   //   O  (1/cool) INT dprior[2] (L/L0)^cool
{
    double  dL1, dL2, t1, t2;
    int     i, k1, k2, v1, v2;
    int     CALLvalue = 0;

    dL1 = dL2 = 0.0;
// Trial Lhood increments
    for( k1 = v1 = 0; v1 < Valency; k1 += nbits[v1++] )
    {
        v2 = Xindex[v1];
        if( v2 >= 0 )                     // pick up Atom1/Atom2 cross term...
        {
            for( k2 = i = 0; i < v2; k2 += nbitx[i++] ) ;
            CALL( Poiss2Marginal(MassInf, cool, q, ProbON,
                                 Mock, Data, Acc, Counts, Foot,
                                 nbits[v1], &ibits[k1], &zbits[k1],
                                 nbitx[v2], &ibitx[k2], &zbitx[k2],
                                 &t1, &t2) )
            dL2 += t2;
        }
        else                              // ... or pick up lone Atom1
        {
            Poiss1Marginal(MassInf, cool, q, ProbON, Mock, Data, Acc, Counts,
                     nbits[v1], &ibits[k1], &zbits[k1], &t1);
            dL2 += t1;
        }
        dL1 += t1;
    }
    for( k2 = v2 = 0; v2 < Valency; k2 += nbitx[v2++] )
    {
        for( v1 = 0; v1 < Valency; v1++ )
            if( Xindex[v1] == v2 )
                break;
        if( v1 == Valency )               // remaining lone Atom2
        {
            Poiss1Marginal(MassInf, cool, q, ProbON, Mock, Data, Acc, Counts,
                     nbitx[v2], &ibitx[k2], &zbitx[k2], &t2);
            dL2 += t2;
        }
    }
    *dLtry1 = dL1;
    *dLtry2 = dL2;
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissInsert1
//
//            Sample fluxes from    Prior * Lhood^cool
//            for one entire atom with Poisson likelihood
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int PoissInsert1(  //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
Rand_t   Rand,     // I O  Random generator state
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
int      Valency,  // I    # fluxes per atom
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
double*  zbits,    // I    Fragment quantities     [<=Ndata]
double*  Flux)     //   O  Sample fluxes           [Valency]
{
    int  k, v;
    int  CALLvalue = 0;

    for( k = v = 0; v < Valency; k += nbits[v++] )
        CALL( Poiss1Sample(MassInf, Rand, cool, q, ProbON, Mock, Data, Acc,
                           Counts, nbits[v], &ibits[k], &zbits[k], &Flux[v]) )
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissInsert2
//
//            Sample fluxes from    Prior * Lhood^cool
//            for two entire atoms with Poisson likelihood
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int PoissInsert2(  //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
Rand_t   Rand,     // I O  Random generator state
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      Valency,  // I    # fluxes per atom
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
double*  zbits,    // I    Fragment quantities     [<=Ndata]
int*     nbitx,    // I    # fragments             [Valency]
int*     ibitx,    // I    Fragment positions      [<=Ndata]
double*  zbitx,    // I    Fragment quantities     [<=Ndata]
int*     Xindex,   // I    Index to cross-terms    [Valency]
double*  Flux1,    //   O  Fluxes of 1st atom      [Valency]
double*  Flux2)    //   O  Fluxes of 2nd atom      [Valency]
{
    int  i, k1, k2, v1, v2;
    int  CALLvalue = 0;

    for( k1 = v1 = 0; v1 < Valency; k1 += nbits[v1++] )
    {
        v2 = Xindex[v1];
        if( v2 >= 0 )
        {
            for( k2 = i = 0; i < v2; k2 += nbitx[i++] ) ;
            CALL( Poiss2Sample(MassInf, Rand, cool, q, ProbON,
                               Mock, Data, Acc, Counts, Foot,
                               nbits[v1], &ibits[k1], &zbits[k1],
                               nbitx[v2], &ibitx[k2], &zbitx[k2],
                               &Flux1[v1], &Flux2[v2]) )
        }
        else
            CALL( Poiss1Sample(MassInf, Rand, cool, q, ProbON,
                               Mock, Data, Acc, Counts,
                               nbits[v1], &ibits[k1], &zbits[k1], &Flux1[v1]) )
    }
    for( k2 = v2 = 0; v2 < Valency; k2 += nbitx[v2++] )
    {
        for( v1 = 0; v1 < Valency; v1++ )
            if( Xindex[v1] == v2 )
                break;
        if( v1 == Valency )
            CALL( Poiss1Sample(MassInf, Rand, cool, q, ProbON,
                               Mock, Data, Acc, Counts,
                               nbitx[v2], &ibitx[k2], &zbitx[k2], &Flux2[v2]) )
    }
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissLhood1
//
//            DELTA( log Lhood )  for inserting or deleting
//            one entire atom with Poisson likelihood
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
double PoissLhood1(//   O  DELTA(logL)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int      Valency,  // I    # fluxes per atom
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
double*  zbits,    // I    Fragment quantities     [<=Ndata]
double*  Flux,     // I    Fluxes of atom          [Valency]
int      flag)     // I    +1 = insert, -1 = delete
{
    double  dL = 0.0;
    int     k, v;

    if( flag < 0 )
        for( v = 0; v < Valency; v++ )
            Flux[v] = -Flux[v];
    for( k = v = 0; v < Valency; k += nbits[v++] )
        dL += poiss1lhood(Flux[v], Mock, Data, Acc,
                          nbits[v], &ibits[k], &zbits[k]);
    if( flag < 0 )
        for( v = 0; v < Valency; v++ )
            Flux[v] = -Flux[v];
    return  dL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PoissLhood2
//
//            DELTA( log Lhood )
//            for inserting two entire atoms with Poisson likelihood
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
double PoissLhood2(//   O  DELTA(logL)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      Valency,  // I    # fluxes per atom
int*     nbits,    // I    # fragments             [Valency]
int*     ibits,    // I    Fragment positions      [<=Ndata]
double*  zbits,    // I    Fragment quantities     [<=Ndata]
int*     nbitx,    // I    # fragments             [Valency]
int*     ibitx,    // I    Fragment positions      [<=Ndata]
double*  zbitx,    // I    Fragment quantities     [<=Ndata]
int*     Xindex,   // I    index to cross-terms    [Valency]
double*  Flux1,    // I    Fluxes of 1st atom      [Valency]
double*  Flux2)    // I    Fluxes of 2nd atom      [Valency]
{
    double  dL = 0.0;
    double  t1, t2;
    int     i, k1, k2, v1, v2;

    for( k1 = v1 = 0; v1 < Valency; k1 += nbits[v1++] )
    {
        v2 = Xindex[v1];
        if( v2 >= 0 )
        {
            for( k2 = i = 0; i < v2; k2 += nbitx[i++] ) ;
            poiss2lhood(Flux1[v1], Flux2[v2], Mock, Data, Acc, Foot,
                        nbits[v1], &ibits[k1], &zbits[k1],
                        nbitx[v2], &ibitx[k2], &zbitx[k2], &t1, &t2);
            dL += t2;
        }
        else
            dL += poiss1lhood(Flux1[v1], Mock, Data, Acc,
                              nbits[v1], &ibits[k1], &zbits[k1]);
    }
    for( k2 = v2 = 0; v2 < Valency; k2 += nbitx[v2++] )
    {
        for( v1 = 0; v1 < Valency; v1++ )
            if( Xindex[v1] == v2 )
                break;
        if( v1 == Valency )
            dL += poiss1lhood(Flux2[v2], Mock, Data, Acc,
                              nbitx[v2], &ibitx[k2], &zbitx[k2]);
    }
    return dL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Poiss1Marginal
//
//                 infinity           -cool*R*z                         Counts
//    (1/cool) log INTEGRAL dz Pr(z) e          ( 1 + z * R/(Mock+Acc) )
//                    z=0
//
//            Counts are integer approximations to cooled data, cool*(Data+Acc)
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int Poiss1Marginal(//   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities R     [nbits]
double*  dLtry)    //   O  (1/cool) INT dprior (L/L0)^cool
{
    double  p;
    double  ProbOFF;
    int     CALLvalue = 0;

    CALL( Poisson1(MassInf, NULL, cool, q, Mock, Data, Acc, Counts,
                   nbits, ibits, zbits, &p) )
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        if( cool * cool * p * p < DBL_EPSILON )
            p *= ProbON;
        else
            p = PLUS(log(ProbOFF), log(ProbON) + cool * p) / cool;
    }
Exit:
    *dLtry = p;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Poiss2Marginal
//
//
//                 infinity           -cool*R*x                         Counts
//    (1/cool) log INTEGRAL dx Pr(x) e          ( 1 + x * R/(Mock+Acc) )
//                    x=0
//   and
//                 infinity                  -cool*u                     Counts
//    (1/cool) log INTEGRAL dxdy Pr(x)Pr(y) e        ( 1 + u/(Mock+Acc) )
//                  x,y = 0
//
//            where  u = R*x + S*y,
//            Counts are integer approximations to cooled data, cool*(Data+Acc)
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int Poiss2Marginal(//   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities R     [nbits]
int      nbitx,    // I    # fragments
int*     ibitx,    // I    Fragment positions        [nbitx]
double*  zbitx,    // I    Fragment quantities S     [nbitx]
double*  dLtry1,   //   O  (1/cool) INT dprior[1] (L/L0)^cool
double*  dLtry2)   //   O  (1/cool) INT dprior[2] (L/L0)^cool
{
    double  p1, p2, p12;
    double  ProbOFF;
    int     CALLvalue = 0;

    CALL( Poisson2(MassInf, NULL, cool, q, Mock, Data, Acc, Counts, Foot,
                   nbits, ibits, zbits, nbitx, ibitx, zbitx, &p1, &p12) )
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        CALL( Poisson1(MassInf, NULL, cool, q, Mock, Data, Acc, Counts,
                       nbitx, ibitx, zbitx, &p2) )
        if( cool * cool * (p1*p1 + p2*p2 + p12*p12) < DBL_EPSILON )
        {
            p12 = ProbON * (ProbOFF * (p1 + p2) + ProbON * p12);
            p1 *= ProbON;
        }
        else
        {
            p1 = PLUS(log(ProbOFF), log(ProbON) + cool*p1);
            p2 = PLUS(log(ProbOFF) + cool*p2, log(ProbON) + cool*p12);
            p12 = PLUS(log(ProbOFF) + p1, log(ProbON) + p2) / cool;
        }
    }
Exit:
    *dLtry1 = p1;
    *dLtry2 = p12;
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Poiss1Sample
//
//            Sample flux z from
//                                -cool*R*z                         Counts
//                         Pr(z) e          ( 1 + z * R/(Mock+Acc) )
//
//            Counts are integer approximations to cooled data, cool*(Data+Acc)
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int Poiss1Sample(  //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
Rand_t   Rand,     // I O  Random generator state
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities R     [nbits]
double*  Flux)     //   O  Sample flux z
{
    double  p, ProbOFF;
    int     k = 1;
    int     CALLvalue = 0;

    *Flux = 0.0;
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        CALL( Poisson1(MassInf, NULL, cool, q, Mock, Data, Acc, Counts,
                       nbits, ibits, zbits, &p) )
        p = log(ProbON) + cool * p;
        if( log(Randouble(Rand)) > p - PLUS(log(ProbOFF), p) )
            k = 0;
    }
    if( k )
        CALL( Poisson1(MassInf, Rand, cool, q, Mock, Data, Acc, Counts,
                       nbits, ibits, zbits, Flux) )
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Poiss2Sample
//
//            Sample single fluxes x,y from
//                           -cool*(R*x+S*y)                            Counts
//               Pr(x)Pr(y) e                ( 1 + (R*x+S*y)/(Mock+Acc) )
//
//            Counts are integer approximations to cooled data, cool*(Data+Acc)
//
// History:   JS   18 Aug 2003, 12 Sep 2003
//-----------------------------------------------------------------------------
static
int Poiss2Sample(  //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive, add 10 for z=0
Rand_t   Rand,     // I O  Random generator state
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double   ProbON,   // I    Pr(individual flux != 0)
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities R     [nbits]
int      nbitx,    // I    # fragments
int*     ibitx,    // I    Fragment positions        [nbitx]
double*  zbitx,    // I    Fragment quantities S     [nbitx]
double*  Flux1,    //   O  Flux x of 1st atom
double*  Flux2)    //   O  Flux y of 2nd atom
{
    double  p, p1, p2, p12, a, b, c, r, ProbOFF;
    int     k = 3;
    int     CALLvalue = 0;

    *Flux1 = *Flux2 = 0.0;
    if( ProbON < 1.0 )
    {
        ProbOFF = 1.0 - ProbON;
        CALL( Poisson2(MassInf, NULL, cool, q, Mock, Data, Acc, Counts,
                       Foot, nbits, ibits, zbits, nbitx, ibitx, zbitx,
                       &p1, &p12) )
        CALL( Poisson1(MassInf, NULL, cool, q, Mock, Data, Acc, Counts,
                       nbits, ibits, zbits, &p2) )
        p   = log(ProbOFF) + log(ProbOFF);
        p1  = log(ProbOFF) + log(ProbON) + cool * p1;
        p2  = log(ProbOFF) + log(ProbON) + cool * p2;
        p12 = log(ProbON) + log(ProbON) + cool * p12;
        a = PLUS(p, p1);
        b = PLUS(p2, p12);
        c = PLUS(a, b);
        p1  = exp(p1 - c);
        p2  = exp(p2 - c);
        p12 = exp(p12 - c);
        r = Randouble(Rand);
        if( r < p12 )
            k = 3;
        else if( r < p12 + p2 )
            k = 2;
        else if( r < p12 + p2 + p1 )
            k = 1;
        else
            k = 0;
    }
    if( k == 1 )
        CALL( Poisson1(MassInf, Rand, cool, q, Mock, Data, Acc, Counts,
                       nbits, ibits, zbits, Flux1) )
    if( k == 2 )
        CALL( Poisson1(MassInf, Rand, cool, q, Mock, Data, Acc, Counts,
                       nbitx, ibitx, zbitx, Flux2) )
    if( k == 3 )
        CALL( Poisson2(MassInf, Rand, cool, q, Mock, Data, Acc, Counts, Foot,
                      nbits, ibits, zbits, nbitx, ibitx, zbitx, Flux1, Flux2) )
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Poisson1
//
// EITHER EVALUATE
//                  infinity          -cool*R*z                         Counts
//     (1/cool) log INTEGRAL dz P(z) e          ( 1 + z * R/(Mock+Acc) )
//                     z=0
// OR SAMPLE FLUX z from
//                                    -cool*R*z                         Counts
//                 Posterior =  P(z) e          ( 1 + z * R/(Mock+Acc) )
//
//            Counts are integer approximations to cooled data, cool*(Data+Acc)
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
int Poisson1(      //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive
Rand_t   Rand,     //(I)   NULL is trial dL, else z for insertion
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities       [nbits]
double*  result)   //   O  if "try", dLtry = (1/cool) INT dprior (L/L0)^cool
                   //      else "insert", z = random flux
{
    int      npoly;         // # Counts covered by "bits"
    double*  poly  = NULL;  // polynomial from "bits"    [0..npoly]
    double   scale;         // dimensional scaling for Footprint
    double   sum;           // polynomial normalisation
    double   r;             // random limit
    double   t;             // temporary
    int      i;             // counter for "bits"
    int      j;             // counter for Counts
    int      k;             // counter for polynomial
    int      CALLvalue = 0;

// MONKEYS
    if( MassInf % 10 == 0 )
    {
// "try"............................
        if( ! Rand )
            *result = poiss1lhood(q, Mock, Data, Acc, nbits, ibits, zbits);
// "insert".........................
        else
            *result = q;
    }

// POSITIVE
    if( MassInf % 10 == 1 )
    {
// Dimensional scaling
        t = 0.0;
        for( i = 0; i < nbits; i++ )
            t += zbits[i];
// Special case of cool = 0
        if( cool * cool < DBL_EPSILON )
        {
// "try"............................
            if( ! Rand )
                *result = - q * t;
// "insert".........................
            else
                *result = - q * log(Randouble(Rand));
        }
        else
// General case
        {
            scale = 1.0 / q + cool * t;
// Get polynomial size from integer-approximated annealed Counts
            npoly = 0;
            for( i = 0; i < nbits; i++ )
                npoly += Counts[ibits[i]];
            CALLOC(poly, 1 + npoly, double)
            for( i = 0; i <= npoly; i++ )
                poly[i] = 0.0;
// Set polynomial, normalised against overflow as exp(s)*poly with SUM(poly)=1
            sum = 0.0;
            npoly = 0;
            poly[0] = 1.0;
            for( i = 0; i < nbits; i++ )
            {
                t = zbits[i] / ((Mock[ibits[i]] + Acc[ibits[i]]) * scale);
                for( j = Counts[ibits[i]]; j > 0; j-- )
                {
                    r = poly[0];
                    for( k = ++npoly; k > 0; k-- )
                        r += poly[k] += t * poly[k-1] * k;
                    sum += log(r);
                    r = 1.0 / r;
                    for( k = 0; k <= npoly; k++ )
                        poly[k] *= r;
                }
            }
// "try"............................
            if( ! Rand )
                *result = (sum - log(q * scale)) / cool;
// "insert".........................
            else
            {
// Sample k from poly
                r = Randouble(Rand);
                t = 0.0;
                for( k = 0; k < npoly; k++ )     // k=npoly is last case
                    if( (t += poly[k]) >= r )
                       break;
// Sample z from Gamma
                *result = Rangamma(Rand, k + 1.0) / scale;
            }
            FREE(poly);
        }
    }
Exit:
    FREE(poly);
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Poisson2
//
// EITHER EVALUATE
//                  infinity          -cool*R*x                         Counts
//     (1/cool) log INTEGRAL dx P(x) e          ( 1 + x * R/(Mock+Acc) )
//                     x=0
//   and
//                  infinity                -cool*u                     Counts
//     (1/cool) log INTEGRAL dxdy P(x)P(y) e        ( 1 + u/(Mock+Acc) )
//                   x,y = 0
// OR SAMPLE FLUX x from
//                                    -cool*R*x                         Counts
//                  Posterior = P(x) e          ( 1 + x * R/(Mock+Acc) )
//   and x,y from
//                                          -cool*u                     Counts
//                  Posterior =   P(x)P(y) e        ( 1 + u/(Mock+Acc) )
// WHERE  u = R*x + S*y.
//            Counts are integer approximations to cooled data, cool*(Data+Acc)
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
int Poisson2(      //   O  0, or -ve error
int      MassInf,  // I    100=monkeys, 101=positive
Rand_t   Rand,     //(I)   NULL is trial dL1 dL2, else (x,y) flux for insertion
double   cool,     // I    Annealing coefficient
double   q,        // I    Flux unit
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int*     Counts,   // I    Annealed data             [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities       [nbits]
int      nbitx,    // I    # fragments
int*     ibitx,    // I    Fragment positions        [nbitx]
double*  zbitx,    // I    Fragment quantities       [nbitx]
double*  result1,  //   O  if "try" dLtry1 = (1/cool) INT dprior[1] (L/L0)^cool
                   //         for "bits" only,        else x = 1st random flux
double*  result2)  //   O  if "try" dLtry2 = (1/cool) INT dprior[2] (L/L0)^cool
                   //         for "bits" and "bitx",  else y = 2nd random flux
{
static const double BIG = DBL_MAX * DBL_EPSILON;
    int      npoly1;        // # Counts covered by "bits" but not "bitx"
    int      npoly2;        // # Counts covered by "bitx" but not "bits"
    int      npoly0;        // # Counts covered by both "bits" and "bitx"
    double*  poly1  = NULL; // polynomial from "bits"              [0...npoly1]
    double*  poly2  = NULL; // polynomial from "bitx"              [0...npoly2]
    double** poly0  = NULL; // polynomial from both     [i][j] (i+j=0...npoly0)
    double   scale1;        // dimensional scaling for "bits" Footprint
    double   scale2;        // dimensional scaling for "bitx" Footprint
    double   s1;            // normalisation of poly1
    double   s2;            // normalisation of poly2
    double   s0;            // normalisation of poly0
    double   r;             // random limit
    double   t1;            // build-up of poly1 and poly0
    double   t2;            // build-up of poly2 and poly0
    double   sum1;          // summation over "bits"
    double   sum2;          // summation over "bitx"
    double   weight;        // row or column weighting of poly0
    double   t;             // temporary
    double   x;             // temporary
    double   y;             // temporary
    int      i;             // index for "bits" (and general)
    int      j;             // index for "bitx"
    int      m;             // 1st index of poly0
    int      n;             // 2nd index of poly0
    int      k;             // temporary
    int      CALLvalue = 0;

// MONKEYS
    if( MassInf % 10 == 0 )
    {
// "try"............................
        if( ! Rand )
            poiss2lhood(q, q, Mock, Data, Acc, Foot,
                        nbits, ibits, zbits, nbitx, ibitx, zbitx,
                        result1, result2);
// "insert".........................
        else
            *result1 = *result2 = q;
    }

// POSITIVE
    if( MassInf % 10 == 1 )
    {
// Dimensional scaling
        t1 = t2 = 0.0;
        for( i = 0; i < nbits; i++ )
            t1 += zbits[i];
        for( j = 0; j < nbitx; j++ )
            t2 += zbitx[j];

// Special case of cool = 0
        if( cool * cool < DBL_EPSILON )
        {
// "try"............................
            if( ! Rand )
            {
                *result1 = -q * t1;
                *result2 = -q * (t1 + t2);
            }
// "insert".........................
            else
            {
                *result1 = - q * log(Randouble(Rand));
                *result2 = - q * log(Randouble(Rand));
             }
        }
        else
// General case
        {
            scale1 = 1.0 / q + cool * t1;
            scale2 = 1.0 / q + cool * t2;
// Get polynomial sizes from integer-approximated annealed Counts
            npoly0 = npoly1 = npoly2 = 0;
            for( j = 0; j < nbitx; j++ )
                Foot[ibitx[j]] = 1.0;
            for( i = 0; i < nbits; i++ )
            {
                k = ibits[i];
                if( Foot[k] != 0.0 )
                {
                    npoly0 += Counts[k];
                    Foot[k] = 0.0;
                }
                else
                    npoly1 += Counts[k];
            }
            for( j = 0; j < nbitx; j++ )
            {
                k = ibitx[j];
                if( Foot[k] != 0.0 )
                {
                    npoly2 += Counts[k];
                    Foot[k] = 0.0;
                }
            }

// Assign polynomials
            CALLOC(poly1, 1 + npoly1, double)
            CALLOC(poly2, 1 + npoly2, double)
            CALLOC(poly0, 1 + npoly0, double*)
            CALLOC(poly0[0], (1 + npoly0) * (2 + npoly0) / 2, double)
            for( i = 1; i <= npoly0; i++ )
                poly0[i] = poly0[i-1] + 2 + npoly0 - i;
            for( i = 0; i <= npoly1; i++ )
                poly1[i] = 0.0;
            for( i = 0; i <= npoly2; i++ )
                poly2[i] = 0.0;
            for( i = 0; i <= npoly0; i++ )
               for( j = 0; i + j <= npoly0; j++ )
                    poly0[i][j] = 0.0;
            poly1[0] = poly2[0] = poly0[0][0] = 1.0;
            s1 = s2 = s0 = 0.0;

// Set polynomials, normalised against overflow as exp(s)*poly with SUM(poly)=1
            npoly0 = npoly1 = npoly2 = 0;
            for( j = 0; j < nbitx; j++ )
               Foot[ibitx[j]] = zbitx[j];
            for( i = 0; i < nbits; i++ )
            {
                k = ibits[i];
                t1 = zbits[i] / ((Mock[k] + Acc[k]) * scale1);
                if( Foot[k] != 0.0 )
                {
                    t2 = Foot[k] / ((Mock[k] + Acc[k]) * scale2);
                    for( j = Counts[k]; j > 0; j-- )
                    {
                        r = poly0[0][0];
                        for( n = ++npoly0; n > 0; n-- )
                        {
                            r += poly0[n][0] += t1 * poly0[n-1][0] * n;
                            for( m = 1; m < n; m++ )
                                r += poly0[n-m][m]
                                  += t1 * poly0[n-m-1][m] * (n-m)
                                   + t2 * poly0[n-m][m-1] * m;
                            r += poly0[0][m] += t2 * poly0[0][m-1] * m;
                        }
                        s0 += log(r);
                        r = 1.0 / r;
                        for( n = 0; n <= npoly0; n++ )
                            for( m = 0; m + n <= npoly0; m++ )
                                poly0[n][m] *= r;
                    }
                    Foot[k] = 0.0;
                }
                else
                {
                    for( j = Counts[k]; j > 0; j-- )
                    {
                        r = poly1[0];
                        for( m = ++npoly1; m > 0; m-- )
                            r += poly1[m] += t1 * poly1[m-1] * m;
                        s1 += log(r);
                        r = 1.0 / r;
                        for( m = 0; m <= npoly1; m++ )
                            poly1[m] *= r;
                    }
                }
            }
            for( j = 0; j < nbitx; j++ )
            {
                k = ibitx[j];
                if( Foot[k] != 0.0 )
                {
                    t2 = zbitx[j] / ((Mock[k] + Acc[k]) * scale2);
                    for( i = Counts[k]; i > 0; i-- )
                    {
                        r = poly2[0];
                        for( m = ++npoly2; m > 0; m-- )
                            r += poly2[m] += t2 * poly2[m-1] * m;
                        s2 += log(r);
                        r = 1.0 / r;
                        for( m = 0; m <= npoly2; m++ )
                            poly2[m] *= r;
                    }
                    Foot[k] = 0.0;
                }
            }

// Weight poly0 rows by SUM[k=0..npoly1] poly1[k] (k+i)!/k!i!
// Get sum1 (usable for dL1)
            sum1 = 0.0;
            for( i = 0; i <= npoly0; i++ )
            {
                weight = 0.0;
                x = 1.0;
                for( k = 0; k <= npoly1; k++ )
                {
                    weight += x * poly1[k];
                    x *= 1.0 + i / (k + 1.0);
                    if( x > BIG )
                        return E_MASSINF_COUNTS;  // incipient overflow
                }
                for( j = 0; i + j <= npoly0; j++ )
                    poly0[i][j] *= weight;
                sum1 += poly0[i][0];
            }

// Weight poly0 columns by SUM[k=0..npoly2] poly2[k] (k+j)!/k!j!
// Get sum2 (usable for dL2)
            sum2 = 0.0;
            for( j = 0; j <= npoly0; j++ )
            {
                weight = 0.0;
                x = 1.0 / sum1;                   // overflow safety
                for( k = 0; k <= npoly2; k++ )
                {
                    weight += x * poly2[k];
                    x *= 1.0 + j / (k + 1.0);
                    if( x > BIG )
                        return E_MASSINF_COUNTS;  // incipient overflow
                }
                for( i = 0; i + j <= npoly0; i++ )
                    sum2 += poly0[i][j] *= weight;
            }

// "try"............................
            if( ! Rand )
            {
                *result1 = (s0 + s1 + log(sum1 / (q * scale1))) / cool;
                *result2 = (s0 + s1 + s2 + log(sum1 / (q * scale1))
                                         + log(sum2 / (q * scale2))) / cool;
            }
// "insert".........................
            else
            {
// Sample i,j from poly0
                r = sum2 * Randouble(Rand);
                t = 0.0;
                for( i = 0; i <= npoly0; i++ )
                {
                    for( j = 0; i + j <= npoly0; j++ )
                        if( (t += poly0[i][j]) >= r )
                            break;
                    if( t >= r )
                        break;
                }
                if( i > npoly0 )  i = j = 0;      // extra protection
// Sample m from poly1
                t = 0.0;
                x = y = 1.0;
                for( m = 0; m <= npoly1; m++ )
                {
                    t += poly1[m] *= x;
                    x *= 1.0 + i / y;
                    y += 1.0;
                }
                t *= Randouble(Rand);
                for( m = 0; m < npoly1; m++ )     // m=npoly1 is last case
                    if( (t -= poly1[m]) <= 0.0 )
                       break;
// Sample n from poly2
                t = 0.0;
                x = y = 1.0;
                for( n = 0; n <= npoly2; n++ )
                {
                    t += poly2[n] *= x;
                    x *= 1.0 + j / y;
                    y += 1.0;
                }
                t *= Randouble(Rand);
                for( n = 0; n < npoly2; n++ )     // n=npoly2 is last case
                    if( (t -= poly2[n]) <= 0.0 )
                        break;
// Sample x,y from Gamma
                *result1 = Rangamma(Rand, m + i + 1.0) / scale1;
                *result2 = Rangamma(Rand, n + j + 1.0) / scale2;
            }
        }
    }
Exit:
    if( poly0 )
        FREE(poly0[0])
    FREE(poly0)
    FREE(poly2)
    FREE(poly1)
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  poiss1lhood
//                               -R*z                   Data+Acc
//            Contribution  log(e    (1 + R*z/(Mock+Acc))       )
//
//            to  DELTA(logL)  from single flux.
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
double poiss1lhood(//   O  DELTA(logL)
double   Flux,     // I    Flux z of atom
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits)    // I    Fragment quantities R     [nbits]
{
    double  dL = 0.0;
    int     i, k;

    for( i = 0; i < nbits; i++ )
    {
        k = ibits[i];
        dL += (Data[k]+Acc[k]) * log(1.0 + Flux * zbits[i] / (Mock[k]+Acc[k]))
             - Flux * zbits[i];
    }
    return dL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  poiss2lhood
//                                -R*x                    Data+Acc
//            Contributions  log(e    (1 + R*x/(Mock+Acc))        )
//
//                          -(R*x+S*y)                          Data+Acc
//            and      log(e          (1 + (R*x+S*y)/(Mock+Acc))        )
//
//            to DELTA(logL) from one, and from two, single fluxes.
//
// History:   JS   18 Aug 2003
//-----------------------------------------------------------------------------
static
void poiss2lhood(
double   Flux1,    // I    Flux x for "bits"
double   Flux2,    // I    Flux y for "bitx"
double*  Mock,     // I    Mock data                 [Ndata]
double*  Data,     // I    Data                      [Ndata]
double*  Acc,      // I    Background                [Ndata]
double*  Foot,     //(I O) Footprint (=0)            [Ndata]
int      nbits,    // I    # fragments
int*     ibits,    // I    Fragment positions        [nbits]
double*  zbits,    // I    Fragment quantities R     [nbits]
int      nbitx,    // I    # fragments
int*     ibitx,    // I    Fragment positions        [nbitx]
double*  zbitx,    // I    Fragment quantities S     [nbitx]
double*  dL1,      //   O  DELTA(logL) for x only
double*  dL2)      //   O  DELTA(logL) for x and y
{
    double sum1, sum2, t;
    int    i, k;

    for( i = 0; i < nbitx; i++ )
       Foot[ibitx[i]] = zbitx[i];
    sum1 = sum2 = 0.0;
    for( i = 0; i < nbits; i++ )
    {
        k = ibits[i];
        sum1 -= Flux1 * zbits[i];
        sum2 -= Flux1 * zbits[i];
        t = (Data[k]+Acc[k]) * log(1.0 + Flux1 * zbits[i] / (Mock[k]+Acc[k]));
        sum1 += t;
        if( Foot[k] == 0.0 )
            sum2 += t;
        else
        {
            t = Flux1 * zbits[i] + Flux2 * Foot[k];
            sum2 += (Data[k] + Acc[k]) * log(1.0 + t / (Mock[k] + Acc[k]));
            Foot[k] = 0.0;
        }
    }
    for( i = 0; i < nbitx; i++ )
    {
        k = ibitx[i];
        sum2 -= Flux2 * zbitx[i];
        if( Foot[k] != 0.0 )
        {
            sum2 += (Data[k] + Acc[k])
                   * log(1.0 + Flux2 * Foot[k] / (Mock[k] + Acc[k]));
            Foot[k] = 0.0;
        }
    }
    *dL1 = sum1;
    *dL2 = sum2;
}

#if PARALLEL
//=============================================================================
//   PARALLEL code replacing DoOperations procedure in BayeSys3
//   History:   Do Kester / John Skilling    4 Nov 2002, 10 Feb 2003
//=============================================================================

//#define THREADS
#define FLOWCHECK 1            // debugging aid
#define E_FLOWCHECK     -230   // parallel flow error

#ifdef THREADS
#include <pthread.h>
#define pthread_create( a, b, c, d )    \
        pthread_create( a, b, (void*(*)(void*))c, (void*)d )
#else
#define pthread_t char                                // not used
#define pthread_create( a, b, c, d )  if( a ) c( d )  // avoid compiler whinge
#define pthread_join( a, b )                          // null statement
#endif

typedef struct             // COLLECTED PARAMETERS
{
    OperStr*   Operations; // all operations
    CommonStr* Common;     // general information
    ObjectStr* Objects;    // all objects
    Node*      Links;      // all linked lists
    int        op;         // current operation #
    int        status;     // 0 = idle, 1 = active, 2 = finished, -ve = error
} WorkStr;

#if FLOWCHECK
#include <stdio.h>
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PrintStart
//        Development diagnostics for start of operation
//-----------------------------------------------------------------------------
static void PrintStart(int slave, WorkStr* Work)
{
    CommonStr* Common  = Work[slave].Common;
    ObjectStr* Objects = Work[slave].Objects;
    OperStr*   Oper    = &Work[slave].Operations[Work[slave].op];
    int        i, flag;

    printf("+ oper%4d on%2d ", Work[slave].op, slave);
    for( i = 0; i < Common->ENSEMBLE; i++ )
    {
        if( (Oper->kType == +1 && Oper->k == i)
         || (Oper->jType == +1 && Oper->j == i)
         || (Oper->iType == +1 && Oper->i == i) )
            printf(" r");
        else if( (Oper->kType == -1 && Oper->k == i)
              || (Oper->jType == -1 && Oper->j == i)
              || (Oper->iType == -1 && Oper->i == i) )
            printf(" w");
        else if( Objects[i].flowcheck )
            printf(" |");
        else
            printf(" .");
    }
    printf("   ");
    flag = 0;
    for( i = 1; i <= PARALLEL; i++ )
    {
        if( Work[i].status == 0 )  printf(".");
        if( Work[i].status == 1 )  printf("#");
        if( Work[i].status == 2 )  printf("-");
        if( Work[i].status )  flag++;
    }
    printf("%2d\n", flag);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  PrintFinish
//        Development diagnostics for finish of operation
//-----------------------------------------------------------------------------
static void PrintFinish(int slave, WorkStr* Work)
{
    CommonStr* Common  = Work[slave].Common;
    ObjectStr* Objects = Work[slave].Objects;
    OperStr*   Oper    = &Work[slave].Operations[Work[slave].op];
    int        i;

    printf("- oper%4d on%2d ", Work[slave].op, slave);
    for( i = 0; i < Common->ENSEMBLE; i++ )
    {
        if( (Oper->kType == +1 && Oper->k == i)
         || (Oper->jType == +1 && Oper->j == i)
         || (Oper->iType == +1 && Oper->i == i) )
            printf(" R");
        else if( (Oper->kType == -1 && Oper->k == i)
              || (Oper->jType == -1 && Oper->j == i)
              || (Oper->iType == -1 && Oper->i == i) )
            printf(" W");
        else if( Objects[i].flowcheck )
            printf(" |");
        else
            printf(" .");
    }
    printf("   ");
    for( i = 1; i <= PARALLEL; i++ )
    {
        if( Work[i].status == 0 )  printf(".");
        if( Work[i].status == 1 )  printf("#");
        if( Work[i].status == 2 )  printf("-");
    }
    printf("\n");
}
#endif

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Do1oper
//        Controller for Do1operation
//        Report return state in  Slave->status
//-----------------------------------------------------------------------------
static void Do1oper( WorkStr*  Slave)
{
    CommonStr* Common     = Slave->Common;        // general information
    ObjectStr* Objects    = Slave->Objects;       // ENSEMBLE of samples
    Node*      Links      = Slave->Links;         // linked lists of labels
    OperStr*   Operations = Slave->Operations;    // all operations
    int        op         = Slave->op;            // this operation
    int        CALLvalue;

    CALLvalue = Do1operation(Operations, Common, Objects, Links, op);
    Slave->status = 2;                            // finished
    if( CALLvalue < 0 )
        Slave->status = CALLvalue;                // error state
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  DoOperations
//        Control operation of batch of parallelisable engine calls
//-----------------------------------------------------------------------------
static int DoOperations( //   O  0, or -ve error
OperStr*   Operations,   // I O  all operations
CommonStr* Common,       // I O  general information
ObjectStr* Objects,      // I O  ENSEMBLE of samples
Node*      Links,        // I O  linked lists of labels
int        nOper)        // I    # operations
{
    int          ENSEMBLE = Common->ENSEMBLE; // # objects
    OperStr*     Oper;                        // & individual operation
    int          op;                          // operation choice
    int          remain;                      // # operations still unfinished
    int          slave;                       // processor id
    int          m;                           // object counter
    WorkStr      Work   [1+PARALLEL];         // Only this many threads...
    pthread_t    threads[1+PARALLEL];         // ...possible at one time
    int*         r     = NULL;
    int*         w     = NULL;
    int          CALLvalue = 0;

// Initialisations
    for( slave = 1; slave <= PARALLEL; slave++ )
    {
        Work[slave].Operations = Operations;
        Work[slave].Common     = Common;
        Work[slave].Objects    = Objects;
        Work[slave].Links      = Links;
        Work[slave].op         = -1;
        Work[slave].status     = 0;
    }
    CALLOC(r, ENSEMBLE, int)
    CALLOC(w, ENSEMBLE, int)

// Set usable ranges Early <= clock <= Late of start/end counts for each object
// clock increments at each start and each finish (start even, finish odd).
// clock= 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17
//        Write---Read----Write---Read----Read----Read----Read----Write---Write
//                               <---any overlapping order--->
// Early= 0       2       4       6       6       6       6      14      16
// Late = 0       2       4      12      12      12      12      14      16
    for( m = 0; m < ENSEMBLE; m++ )
    {
        w[m] = -2;
        r[m] = 0;
    }
    for( op = 0; op < nOper; op++ )
    {
        Oper = &Operations[op];
        if( Oper->iType )
        {
            m = Oper->i;
            Oper->iEarly = (Oper->iType == -1) ? (w[m] = r[m]) : w[m] + 2;
            r[m] += 2;
        }
        if( Oper->jType )
        {
            m = Oper->j;
            Oper->jEarly = (Oper->jType == -1) ? (w[m] = r[m]) : w[m] + 2;
            r[m] += 2;
        }
        if( Oper->kType )
        {
            m = Oper->k;
            Oper->kEarly = (Oper->kType == -1) ? (w[m] = r[m]) : w[m] + 2;
            r[m] += 2;
        }
    }
    for( m = 0; m < ENSEMBLE; m++ )
    {
        w[m] = r[m];
        r[m] -= 2;
    }
    for( op = nOper-1; op >= 0; op-- )
    {
        Oper = &Operations[op];
        if( Oper->iType )
        {
            m = Oper->i;
            Oper->iLate = (Oper->iType == -1) ? (w[m] = r[m]) : w[m] - 2;
            r[m] -= 2;
        }
        if( Oper->jType )
        {
            m = Oper->j;
            Oper->jLate = (Oper->jType == -1) ? (w[m] = r[m]) : w[m] - 2;
            r[m] -= 2;
        }
        if( Oper->kType )
        {
            m = Oper->k;
            Oper->kLate = (Oper->kType == -1) ? (w[m] = r[m]) : w[m] - 2;
            r[m] -= 2;
        }
    }

#if FLOWCHECK
// Optional flowchart checks
    for( m = 0; m < ENSEMBLE; m++ )
        Objects[m].flowcheck = 0;
#endif
// Operation counts
    for( m = 0; m < ENSEMBLE; m++ )
        Objects[m].clock = 0;
// Operations not done (invalid processor #)
    for( op = 0; op < nOper; op++ )
        Operations[op].slave = -1;
// Slaves idle
    for( slave = 1; slave <= PARALLEL; slave++ )
        Work[slave].status = 0;
// Counter
    remain = nOper;

// Cycle round the slaves
    for( slave = 1; remain; slave = 1 + (slave % PARALLEL) )
    {
// Keep monitoring slaves
        if( Work[slave].status < 0 )    // Error state
        {
            CALLvalue = Work[slave].status;
            goto Exit;
        }
        if( Work[slave].status == 2 )   // This slave finished, remove flags
        {
            pthread_join( threads[slave], NULL );           // **JOIN** slave
            Work[slave].status = 0;
            op = Work[slave].op;
            Oper = &Operations[op];
// Update operation count
            if( Oper->iType )    Objects[Oper->i].clock ++;
            if( Oper->jType )    Objects[Oper->j].clock ++;
            if( Oper->kType )    Objects[Oper->k].clock ++;
            remain --;
#if FLOWCHECK
// Optional status checks
            if( Oper->iType == -1 )
            {
                if( Objects[Oper->i].flowcheck != -1 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->i].flowcheck = 0;
            }
            if( Oper->iType == +1 )
            {
                if( Objects[Oper->i].flowcheck <= 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->i].flowcheck --;
            }
            if( Oper->jType == -1 )
            {
                if( Objects[Oper->j].flowcheck != -1 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->j].flowcheck = 0;
            }
            if( Oper->jType == +1 )
            {
                if( Objects[Oper->j].flowcheck <= 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->j].flowcheck --;
            }
            if( Oper->kType == -1 )
            {
                if( Objects[Oper->k].flowcheck != -1 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->k].flowcheck = 0;
            }
            if( Oper->kType == +1 )
            {
                if( Objects[Oper->k].flowcheck <= 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->k].flowcheck --;
            }
            PrintFinish(slave, Work);
#endif
        }
        if( Work[slave].status )
            continue;                   // This slave is not free.  Else...
// Seek usable operation
        for( op = 0; op < nOper; op++ )
        {
            Oper = &Operations[op];
// Already done?
            if( Oper->slave >= 0 )
                continue;
// Out of range?
            if( Oper->iType )
                if( Objects[Oper->i].clock < Oper->iEarly
                 || Objects[Oper->i].clock > Oper->iLate )  continue;
            if( Oper->jType )
                if( Objects[Oper->j].clock < Oper->jEarly
                 || Objects[Oper->j].clock > Oper->jLate )  continue;
            if( Oper->kType )
                if( Objects[Oper->k].clock < Oper->kEarly
                 || Objects[Oper->k].clock > Oper->kLate )  continue;
// Handshake slave<-->operation
            Oper->slave    = slave;
            Work[slave].op = op;
// Slave now active
            Work[slave].status = 1;
// Update operation count
            if( Oper->iType )    Objects[Oper->i].clock ++;
            if( Oper->jType )    Objects[Oper->j].clock ++;
            if( Oper->kType )    Objects[Oper->k].clock ++;
#if FLOWCHECK
// Optional flowchart checks
            if( Oper->iType == -1 )
            {
                if( Objects[Oper->i].flowcheck != 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->i].flowcheck = -1;
            }
            if( Oper->iType == +1 )
            {
                if( Objects[Oper->i].flowcheck < 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->i].flowcheck ++;
            }
            if( Oper->jType == -1 )
            {
                if( Objects[Oper->j].flowcheck != 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->j].flowcheck = -1;
            }
            if( Oper->jType == +1 )
            {
                if( Objects[Oper->j].flowcheck < 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->j].flowcheck ++;
            }
            if( Oper->kType == -1 )
            {
                if( Objects[Oper->k].flowcheck != 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->k].flowcheck = -1;
            }
            if( Oper->kType == +1 )
            {
                if( Objects[Oper->k].flowcheck < 0 )
                { CALLvalue = E_FLOWCHECK;    goto Exit; }
                Objects[Oper->k].flowcheck ++;
            }
            PrintStart(slave, Work);
#endif
// Do usable operation. A free processor is guaranteed.     **CREATE** slave.
            pthread_create( &threads[slave], NULL, Do1oper, &Work[slave] );
            break;
        }
    }
Exit:
    FREE(w)
    FREE(r)
    return CALLvalue;
}
#endif
