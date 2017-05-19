//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//            Bayesian Inference
// 
// Filename:  bayesapp.c
//
// Purpose:   Test operation of BayeSys with
//            UserEmpty,UserTry1,UserTry2,UserInsert1,UserInsert2,UserDelete1.
//=============================================================================
//
// BayeSys3 can operate with arbitrary likelihood functions.
// In this example, an atom has 2 coordinates interpreted as
//            location in [0,1] = Cube[0]  (restricted to [0.0.75])
// and
//         flux in [0,infinity) = -10.0*log(Cube[1])
// to correspond with testing the operation of MassInf with FluxUnit0=10
//
// These are only example programs, involving inefficient calculations that
// repeatedly build an object from scratch with "UserBuild" instructions.
//
//=============================================================================
// History:   JS       2 Jan 2002, 4 Jan 2003, 10 Feb 2003, 20 Aug 2003
//            Copyright (c) 2002,2003 Maximum Entropy Data Consultants Ltd.
//-----------------------------------------------------------------------------

#include <math.h>
#include "bayesys3.h"

static int UserBuild(double*, CommonStr*, ObjectStr*, int, int);

//=============================================================================
//                           SIMPLE CODE
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserEmpty
//
// Purpose:   Set Lhood = logLikelihood(coords) for empty sample
//            with 0 atoms in Object, and initialise its other information.
//
//            Lhood := log L(empty)
//
//            I have already put the number 0 in Object->Natoms,
//            so you can use a simple "UserBuild" instruction.
//-----------------------------------------------------------------------------
// 
int UserEmpty(        //   O  >=0, or -ve error code
double*    Lhood,     //   O  loglikelihood
CommonStr* Common,    // I O  general information
ObjectStr* Object)    // I O  sample object, output new Lhood
{
    UserBuild(Lhood, Common, Object, 0, 1);
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserTry1
//
// Purpose:   d(logLikelihood(coords)) after supposedly adding one new atom
//            to Object, without any side effects on it.
//
//            If Valency = 0,
//              dLtry := d logL(...,x)
//            else mimic the above with                                cool
//              dLtry := (1/cool) log INTEGRAL dPrior(z) dlogL(...,x,z)
//
//            I have already put the new atom x after the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms ], so
//            if Valency=0 you can use simple "UserBuild" instructions.
//-----------------------------------------------------------------------------
// 
int UserTry1(         //   O  +ve = OK, 0 = DO NOT USE, -ve = error
double*    dLtry,     //   O  trial d(logLikelihood) value
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I    sample object (DO NOT UPDATE Lhood)
{
    int    Natoms = Object->Natoms;
    double Lhood  = Object->Lhood;                      // existing Lhood
    double Ltry;
    int    OK;

    *dLtry = 0.0;
    OK = UserBuild(&Ltry, Common, Object, Natoms+1, 0); // trial
    if( OK > 0 )
        *dLtry = Ltry - Lhood;                          // increment
    return  OK;                                         // OK?
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserTry2
//
// Purpose:   d(logLikelihood(coords)) after supposedly adding one,
//            then two new atoms to Object, without any side effects on it.
//
//   If Valency = 0,
//     dLtry1 := d logL(...,x1)
//
//     dLtry2 := d logL(...,x1,x2)
//
//   else mimic the above with                                    cool
//     dLtry1 := (1/cool) log INTEGRAL dPrior(z1) dlogL(...,x1,z1)
//                                                                         cool
//     dLtry2 := (1/cool) log INTEGRAL dPrior(z1,z2) dlogL(...,x1,z1,x2,z2)
//
//            I have already put the new atoms x1,x2 after the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms ] and
//                               at Object->Cubes[ Object->Natoms + 1 ]
//            so if Valency=0 you can use simple "UserBuild" instructions.
//-----------------------------------------------------------------------------
// 
int UserTry2(         //   O  +ve = OK, 0 = DO NOT USE, -ve = error
double*    dLtry1,    //   O  trial d(logLikelihood) value for 1st atom
double*    dLtry2,    //   O  trial d(logLikelihood) value for both atoms
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I    sample object (DO NOT UPDATE Lhood)
{
    int    Natoms = Object->Natoms;
    double Lhood  = Object->Lhood;                           // existing Lhood
    double Ltry1, Ltry2;
    int    OK;

    *dLtry1 = *dLtry2 = 0.0;
    OK = UserBuild(&Ltry1, Common, Object, Natoms+1, 0);     //trial for 1 more
    if( OK > 0 )
    {
        *dLtry1 = Ltry1 - Lhood;                             // increment for 1
        OK = UserBuild(&Ltry2, Common, Object, Natoms+2, 0); //trial for 2 more
        if( OK > 0 )
            *dLtry2 = Ltry2 - Lhood;                         // increment for 2
    }
    return  OK;                             // return OK only if both trials OK
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserInsert1
//
// Purpose:   Insert 1 new atom into Object, keeping it up-to-date, and
//            set d(loglikelihood(coords)). 
//
//            If Valency = 0,
//              dL := d logL(...,x)
//
//            else
//                                                cool
//              sample z from  Prior(z) L(...,x,z)
//
//              and set    dL := d logL(...,x,z)  at the sampled z.
//
//            I have already put the new atom x at the last location of
//            the updated Atoms list, at Object->Cubes[ Object->Natoms - 1 ],
//-----------------------------------------------------------------------------
// 
int UserInsert1(      //   O  >=0, or -ve error code
double*    dL,        //   O  d(loglikelihood)
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I O  sample object
{
    int    Natoms = Object->Natoms;              // new number
    double Lold   = Object->Lhood;               // not yet updated
    double Lnew;

    UserBuild(&Lnew, Common, Object, Natoms, 1); // new updated state
    *dL = Lnew - Lold;                           // I will update Object->Lhood
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserInsert2
//
// Purpose:   Insert 2 new atoms into Object, keeping it up-to-date, and
//            set d(loglikelihood(fluxes)). 
//
//            If Valency = 0,
//              dL := d logL(...,x1,x2)
//
//            else
//                                                                cool
//              sample z1,z2 from  Prior(z1,z2) L(...,x1,z1,x2,z2)
//
//              and set  dL := d logL(...,x1,z1,x2,z2)  at the sampled z1,z2.
//
//            I have already put the new atoms x1,x2 at the last location
//            of the Atoms list, at Object->Cubes[ Object->Natoms - 2 ] and
//                               at Object->Cubes[ Object->Natoms - 1 ]
//            so if Valency=0 you can use a simple "UserBuild" instruction.
//-----------------------------------------------------------------------------
// 
int UserInsert2(      //   O  >=0, or -ve error code
double*    dL,        //   O  d(loglikelihood)
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I O  sample object
{
    int    Natoms = Object->Natoms;              // new number
    double Lold   = Object->Lhood;               // not yet updated
    double Lnew;

    UserBuild(&Lnew, Common, Object, Natoms, 1); // new updated state
    *dL = Lnew - Lold;                           // I will update Object->Lhood
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserDelete1
//
// Purpose:   Delete 1 old atom from Object, keeping it up-to-date, and
//            set d(loglikelihood(fluxes)). 
//  
//            dL := d logL(...)
//
//            I have already put the old atom after the last location of
//            the updated Atoms list, at Object->Cubes[ Object->Natoms ],
//            so if Valency=0 you can use a simple "UserBuild" instruction.
//-----------------------------------------------------------------------------
// 
int UserDelete1(      //   O  >=0, or -ve error code
double*    dL,        //   O  d(loglikelihood)
CommonStr* Common,    // I    general information
ObjectStr* Object)    // I O  sample object
{
    int    Natoms = Object->Natoms;              // new number
    double Lold   = Object->Lhood;               // not yet updated
    double Lnew;

    UserBuild(&Lnew, Common, Object, Natoms, 1); // new updated state
    *dL = Lnew - Lold;                           // I will update Object->Lhood
    return 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserBuild
//
// Purpose:   Build Lhood, and any other data info such as Mock from scratch.
//            Can be called with any number of atoms
//
//            This example implementation has
//            Mock[0] = SUM[atoms] flux ,             flux = - 10*log(Cube[1])
//            Mock[1] = SUM[atoms] flux * location ,  location = Cube[0]
//-----------------------------------------------------------------------------
// 
static int UserBuild( //   O  +ve = OK, 0 = DO NOT USE, -ve = error
double*    Lhood,     //   O  loglikelihood 
CommonStr* Common,    // I    General information
ObjectStr* Object,    // I(O) Cubes in, perhaps Mock out
int        Natoms,    // I    # atoms
int        output)    // I    +ve if Mock to be written out, -ve if not
{
const double Sqrt2Pi = 2.50662827463100050240;
    double** Cubes  = Object->Cubes;  // I    Cubes in [0,1)  [Natoms][Ndim]
    int      Ndata  = Common->Ndata;  // I    # data
    double*  Data   = Common->Data;   // I    Data            [Ndata]
    double*  Acc    = Common->Acc;    // I    Accuracies      [Ndata]
    double   Mock[2];                 // Local mock data      [Ndata]
    double   Lnorm  = 0.0;            // Likelihood normalisation
    double   C      = 0.0;            // chisquared
    double   flux;
    double*  Cube;
    int      k;

// Accumulate mock data
    Mock[0] = Mock[1] = 0.0;
    for( k = 0; k < Natoms; k++ )
    {
        Cube = Cubes[k];
        if( Cube[0] > 0.75 )                    // restrict x to range [0,0.75]
            return 0;                           // or force Lhood = -LARGE
        flux = -10.0 * log(Cube[1]);            // [0,infinity) with q=10
        Mock[0] += flux;                        // Mock[0] = SUM( flux )
        Mock[1] += flux * Cube[0];              // Mock[1] = SUM( flux * x )
    }
// Accumulate logLikelihood
    for( k = 0; k < Ndata; k++ )
    {
        if( Acc[k] > 0.0 )
        {
            C += (Mock[k] - Data[k]) * Acc[k] * Acc[k] * (Mock[k] - Data[k]);
            Lnorm += log(Acc[k] / Sqrt2Pi);     // (fussy to include this)
        }
    }
    *Lhood = Lnorm - C / 2.0;
    if( output )
        for( k = 0; k < Ndata; k++ )
            Object->Mock[k] = Mock[k];
    return 1;
}

//=============================================================================
//       Dummy procedure to link to BayeSys3 when not running MassInf
//=============================================================================
int UserFoot(
double*    Cube,
CommonStr* Common,
int*       ibits,
double*    zbits,
int*       nbits)
{
    return (Cube && Common && ibits && zbits && nbits) ? 1 : 1;
}
