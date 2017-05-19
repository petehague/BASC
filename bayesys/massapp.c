//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//            Massive Inference
// 
// Filename:  massapp.c
// 
// Purpose:   Test operation of MassInf with UserFoot procedure
//=============================================================================
//
// MassInf operates with linear data.
// In this example, an atom has just one (Ndim = 1) spatial dimension
//            location in [0,1] = Cube[0]  (restricted to [0.0.75])
// augmented by a single flux (Valency = 1) stored in Cube[.][1]
//
//=============================================================================
// History:   JS       2 Jan 2002, 16 Oct 2002, 5 Feb 2003, 20 Aug 2003
//            Copyright (c) 2002,2003 Maximum Entropy Data Consultants Ltd.
//-----------------------------------------------------------------------------

#include "bayesys3.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  UserFoot
// 
// Purpose:   Mock-data footprint from unit flux of single atom.
//            No more than Ndata fragments must be returned.
//                nbits[ Valency ]       elements obey nbits[.] >= 0
//                ibits[ SUM(nbits[.]) ] elements obey 0 <= ibits[.] < Ndata
//                zbits[ SUM(nbits[.]) ]
// 
// History:   JS          12 Jan 1998, 4 Feb 1998, 9 Feb 1999, 2 Jan 2002
//-----------------------------------------------------------------------------
// 
int UserFoot(      //   O  +ve = OK
                   //       0  = DO NOT USE this position
                   //      -ve = ERROR
double*    Cube,   // I   Atom position, here only 1 dimension in Cube[0]
CommonStr* Common, // I   Definition of footprint
int*       ibits,  //   O Fragment coords on data, serially per valency
double*    zbits,  //   O Fragment quantities,     serially per valency
int*       nbits)  //   O # fragments >= 0, for each valency, here 1 only
{
    if( Cube[0] > 0.75 )                     // Restrict x to range [0, 0.75]
        return 0;
    nbits[0] = Common->Ndata;                // Valency #0 sees 2 data
    ibits[0] = 0;    zbits[0] = 1.0;         // Data[0] is SUM( flux )
    ibits[1] = 1;    zbits[1] = Cube[0];     // Data[1] is SUM( x * flux )
                                             // No further valencies
    return 1;
}

//=============================================================================
//     Dummy procedures to link to BayeSys3 when running pure MassInf
//=============================================================================
int UserEmpty( double* Lhood, CommonStr* Common, ObjectStr* Object)
{
    *Lhood = 0.0;
    return (Common && Object) ? 1 : 1;
}
int UserTry1(double* Ltry, CommonStr* Common, ObjectStr* Object)
{
    *Ltry = 0.0;
    return (Common && Object) ? 1 : 1;
}
int UserTry2(double* L1, double* L2, CommonStr* Common, ObjectStr* Object)
{
    *L1 = *L2 = 0.0;
    return (Common && Object) ? 1 : 1;
}
int UserInsert1(double* dL, CommonStr* Common, ObjectStr* Object)
{
    *dL = 0.0;
    return (Common && Object) ? 1 : 1;
}
int UserInsert2(double* dL, CommonStr* Common, ObjectStr* Object)
{
    *dL = 0.0;
    return (Common && Object) ? 1 : 1;
}
int UserDelete1(double* dL, CommonStr* Common, ObjectStr* Object)
{
    *dL = 0.0;
    return (Common && Object) ? 1 : 1;
}

