//=============================================================================
//               BAYESYS3/MASSINF  BEGINNER PROGRAM  toy.c
//
// If BayeSys operation (bayestoy), link
//    bayestoy.c = {main,UserMonitor} , userstr.h                        *USER*
//       |
//      bayesys3.[ch] --- {random.[ch], hilbert.[ch]}                  *SYSTEM*
//         |
//        bayesapp.c = {UserEmpty,UserInsert1,UserInsert2,
//                      UserDelete1,UserTry1,UserTry2}                   *USER*
//                      with dummy UserFoot
//
//
// If MassInf operation (masstoy or poisstoy), link
//    masstoy.c = {main,UserMonitor} , userstr.h                         *USER*
//       |
//      bayesys3.[ch] --- random.[ch], hilbert.[ch]                    *SYSTEM*
//         |
//        massapp.c = {UserFoot}                                         *USER*
//                    with dummy UserEmpty,UserInsert1,UserInsert2,
//                               UserDelete1,UserTry1,UserTry2
//
// Programs bayestoy.c, masstoy.c, poisstoy differ only in the MASSINF setting.
//=============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bayesys3.h"
#include "userstr.h"

#define MASSINF   1    // -ve is BayeSys,  >= 0 is MassInf

int main()
{
// Simulation data are {SUM(f), SUM(x.f)}, so that
// Data[0]  =  zeroth moment over x  =  flux
// Data[1]  =  first  moment over x  =  flux * mean location
    int Ndata = 2;
    double Data[2] = {40.0, 12.0}; // flux = 40.0,   location = 12/40 = 0.3
    double Acc [2] = { 1.0,  1.0}; // noise is plus or minus 1/1.0

    int Ncell = 4;                 // collect results with 4-cell resolution

// Addressing
    int           ENSEMBLE = 10;         //  >= 1
    ObjectStr     Objects   [10];        // [ENSEMBLE]
    CommonStr     Common[1];
    UserCommonStr UserCommon[1];
    double        Mockbar[2];            // [Ndata]
    double        PrPos[4];              // [Ncell]
    double        Objbar[4];             // [Ncell]
    double        Objdev[4];             // [Ncell]
    double        Z;
    int           k;
    int           code;

// Initialise BayeSys
    Common->MinAtoms    = 1;          // >= 1
    Common->MaxAtoms    = 0;          // >= MinAtoms, or 0 = infinity
    Common->Alpha       = -1.0;       // +ve for Poisson, -ve for geometric
    Common->ENSEMBLE    = ENSEMBLE;   // # objects in ensemble
    Common->Method      = -1;         // Algorithm method
    Common->Rate        = 0.1;        // Speed of calculation (dimensionless)
    Common->Iseed       = 4321;       // Random seed, -ve is time seed
#if MASSINF < 0                       // BayeSys...
    Common->Ndim        = 2;          // dimension = # coordinates
    Common->Valency     = 0;          // # MassInf fluxes per atom
    for( k = 0; k < ENSEMBLE; k++ )   // UserBuild uses mock data
        Objects[k].Mock = (double*)malloc(Ndata * sizeof(double));
#else                                 // MassInf...
    Common->Ndim        = 1;          // dimension = # coordinates
    Common->Valency     = 1;          // # MassInf fluxes per atom
    Common->MassInf     = MASSINF;    // 0=monkeys,1=pos,2=posneg,3=gauss,
                                      // Add 100 for Poisson data
    Common->ProbON      = 1.0;        // Prob(individual flux != 0)
    Common->FluxUnit0   = 10.0;       // +ve assigned, 0 auto, -ve flexible
#endif
// Initialise Likelihood
    Common->Ndata    = Ndata;         // # data
    Common->Data     = Data;          // data                         [Ndata]
    Common->Acc      = Acc;           // accuracies = 1/sigma         [Ndata]
// Initialise statistics
    Common->UserCommon  = (void*)UserCommon;
    UserCommon->Ncell   = Ncell;      // # cells
    UserCommon->Mockbar = Mockbar;    // <mock data>                  [Ndata]
    UserCommon->PrPos   = PrPos;      // Pr(cell occupied>            [Ncell]
    UserCommon->Objbar  = Objbar;     // <intensity per cell>         [Ncell]
    UserCommon->Objdev  = Objdev;     // std.dev. intensity per cell  [Ncell]
    UserCommon->Nsample = 0;
    UserCommon->atoms   = 0.0;
    for( k = 0; k < Ndata; k++ )
        Mockbar[k] = 0.0;
    for( k = 0; k < Ncell; k++ )
        Objbar[k] = Objdev[k] = PrPos[k] = 0.0;

// Compute ...................................................................
    code = BayeSys3(Common, Objects);
    printf("\nBayeSys3 return code %d\n", code);
// ................................................................... compute

// Normalise all desired statistics
    Z = Common->ENSEMBLE * UserCommon->Nsample;
    UserCommon->atoms /= Z;                                     // <# atoms>
    for( k = 0; k < Ndata; k++ )
        Mockbar[k] /= Z;                                        // <mock data>
    for( k = 0; k < Ncell; k++ )
    {
        PrPos[k] = 0.01*(double)(int)(0.5 + 100.0*PrPos[k]/Z);  // Pr(positive)
        Objbar[k] /= Z;                                         // <intensity>
        Objdev[k] /= Z;
        Objdev[k] -= Objbar[k] * Objbar[k];
        Objdev[k]  = (Objdev[k] > 0.0) ? sqrt(Objdev[k]) : 0.0; // std.dev.
    }

// Display results
    printf("    Random seed was   %d\n"  , Common->Iseed);
    printf("    Success / CPU = %g / %g = %6.4f\n",
                    Common->Success, Common->CPU,
                    Common->Success / Common->CPU);
    printf("    # ENSEMBLE samples was %d\n"  ,  UserCommon->Nsample);
    printf("    < Atoms >   =%10.2f\n",          UserCommon->atoms);
    printf("    Evidence    =%10.2f (log[e])\n", Common->Evidence);
    printf("    Information =%10.2f (log[e])\n", Common->Information);

    printf("\n  Cell    Total +- StdDev    Pr(+)%%\n");
    for( k = 0; k < Ncell; ++k )
        printf("%5d %9.2f %9.2f %7.0f\n", k,Objbar[k],Objdev[k],100.*PrPos[k]);
    printf("\n   Data    <Mock>\n");
    for( k = 0; k < Ndata; ++k )
        printf("%8.2f %8.2f\n", Data[k], Mockbar[k]);
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserMonitor
//
// Purpose:   I have provided a new ensemble of objects.
//
//        1.  For probabilistic exploration, restrict Common->cool to <= 1.
//                    cool = 0 is prior, seen at the first call.
//                0 < cool < 1 is "burn-in"
//                    cool = 1 is evolve posterior,
//                    cool -> infinity is maximum likelihood
//            Return with a POSITIVE return code to exit BayeSys,
//            usually after the end of an evolution stage with cool = 1.
//            You can use Common->Nsystem, which counts calls to UserMonitor.
//
//        2.  You should reset any nuisance parameters x in UserCommon and/or
//            each of your UserObject structures, preferably by sampling from
//            their conditonal probabilities Pr(x|Objects) and/or Pr(x|Object).
//
//        3.  Collect your statistics and display your diagnostics. 
//-----------------------------------------------------------------------------
//
int UserMonitor(            //   O  0 = continue, +ve = finish, -ve = abort
      CommonStr* Common,    // I    Full general information
      ObjectStr* Objects)   // I    Full ensemble of sample objects  [ENSEMBLE]
{
// control
    int     Nrem;    // # to go, as (2 * Nsample++) catches up with Nsystem++
// desired statistics
    int     Ndata   = Common->Ndata;
    UserCommonStr* UserCommon = (UserCommonStr*)Common->UserCommon;
    int     Ncell   = UserCommon->Ncell;
    double* Mockbar = UserCommon->Mockbar; // < mock data >             [Ndata]
    double* PrPos   = UserCommon->PrPos;   // Pr(cell occupied)         [Ncell]
    double* Objbar  = UserCommon->Objbar;  // < intensity per cell >    [Ncell]
    double* Objdev  = UserCommon->Objdev;  // std.dev.intensity per cell[Ncell]
// local
    ObjectStr* Object;         // individual object
    double*    Cube;           // individual atom coordinates in cube
    double     z;              // gives individual atom flux
    double     flux;           // cell-based flux from Cube[1]
    int        cell;           // cell number from Cube[0]
    int        k;              // object counter
    int        r;              // atom counter
    int        j;              // general counter

//.............................................................................
// Re-set any user parameters in UserCommon or UserObject (none in this code).
// You will not need to reset the loglikelihoods.

//.............................................................................
// USER IMPOSES FINAL TEMPERATURE AND TERMINATION CONDITION !
    if( Common->cool >= 1.0 )
        Common->cool = 1.0;

    Nrem = Common->Nsystem - 2 * UserCommon->Nsample;  // or any other setting

// Always print a line of diagnostics (cool, Nrem, logLikelihood[0])
    fprintf(stderr,                     // (stderr flushes output immediately)
      "%10.6f %5d %14.4f\r", Common->cool, Nrem, Objects[0].Lhood);
//.............................................................................
// Set any NUISANCE PARAMETERS you may have here
// Could over-ride system's choices of FluxUnit (both Common and Objects) here
// Could CALIBRATE data here (corrupting Evidence and Information)
//.............................................................................
// BURN-IN
    if( Common->cool < 1.0 )  // final temperature not reached..
        return 0;             // ...normal return, keep going
//.............................................................................
// EVOLUTION: Accumulate any desired statistics
    for( k = 0; k < Common->ENSEMBLE; k++ )
    {
        Object = &Objects[k];
  // # atoms
        UserCommon->atoms += Object->Natoms;
  // mock data
        for( j = 0; j < Ndata; j++ )
            Mockbar[j] += Object->Mock[j];
  // object calculated on Ncell cells
        for( j = 0; j < Ncell; j++ )
        {
            flux = 0.0;
            for( r = 0; r < Object->Natoms; r++ )
            {
                Cube = Object->Cubes[r];
                cell = (int)(Ncell * Cube[0]);
                if( cell == j )
                {
                     z = Cube[1];
                     if( Common->Valency )
                         flux += z;                // pure MassInf operation
                     else
                         flux += -10.0 * log(z);   // direct BayeSys, q=10
                }
            }
            if( flux > 0.0 )
                PrPos[j]  += 1.0;                  // sum(occupied)
            Objbar[j] += flux;                     // sum intensity
            Objdev[j] += flux * flux;              // sum intensity squared
        }
    }
  // # samples of full ensemble
    UserCommon->Nsample ++;

    if( Nrem > 0 )                                 // not finished...
        return 0;                                  // ..normal return
//.............................................................................
// EXIT WITH POSITIVE FINISH CODE WHEN EVOLUTION IS DEEMED COMPLETE
    return 1;
}
