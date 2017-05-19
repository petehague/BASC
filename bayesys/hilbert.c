//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Filename:  hilbert.c
// 
// Purpose:   Hilbert and Linked-list utility procedures for BayeSys3.
// 
// History:   TreeSys.c   17 Apr 1996 - 31 Dec 2002
//            Peano.c     10 Apr 2001 - 11 Jan 2003
//            merged       1 Feb 2003
//            Arith debug 28 Aug 2003
//            Hilbert.c   14 Oct 2003
//                         2 Dec 2003
//-----------------------------------------------------------------------------
/*
    Copyright (c) 1996-2003 Maximum Entropy Data Consultants Ltd,
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
#include <stdlib.h>
#include "hilbert.h"

typedef Atom*  pAtom;          // Atom address
typedef Node*  pNode;          // Node address

// Internal prototypes
static void   LinetoTranspose(coord_t*, coord_t*, int, int);
static void   TransposetoLine(coord_t*, coord_t*, int, int);
static void   TransposetoAxes(coord_t*, int, int);
static void   AxestoTranspose(coord_t*, int, int);
static int    ResetLink   (pNode);
static void   Balance     (pNode);

// Internal macros
#undef  CALLOC    // allocates vector p[0:n-1] of type t
#define CALLOC(p,n,t) {p=NULL;\
if((n)>0&&!(p=(t*)calloc((size_t)(n),sizeof(t))))\
{CALLvalue=E_MALLOC;goto Exit;}/*printf("%p %d\n",p,(size_t)(n)*sizeof(t));*/}

#undef  FREE      // frees CALLOC or REALLOC or NULL vector p[0:*], sets p=NULL
#define FREE(p) {if(p){/*printf("%p -1\n",p);*/(void)free((void*)p);} p=NULL;}

#undef  CALL      // catches negative error codes
#define CALL(x) {if( (CALLvalue = (x)) < 0 ) goto Exit;}

//=============================================================================
//                    Composite-integer arithmetic library
//=============================================================================
//
//  A composite-integer is a multi-word unsigned integer "Label" stored
//  "big endian" in N conventional unsigned integers with [0] high.
//        ___________________________________________________
//       |            |            |            |            |
//       |  Label[0]  |  Label[1]  |    ....    | Label[N-1] |
//       |____________|____________|____________|____________|
//            high                                   low

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  CmpLabel
//
// Purpose:   Compare labels by generating
//                                    +1 if u > v,
//                    sign(u - v)  =   0 if u = v,
//                                    -1 if u < v.
//
// History:   John Skilling   12 Apr 2001
//-----------------------------------------------------------------------------
// 
int  CmpLabel(       //   O  comparison
coord_t* u,          // I    composite integer ([0] high)
coord_t* v,          // I    composite integer ([0] high)
int      Ndim)       // I    dimension
{
    int  j;
    for( j = 0; j < Ndim; j++ )
        if( u[j] < v[j] )
            return -1;
        else if( u[j] > v[j] )
            return 1;
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  AddLabel
//
// Purpose:   Set     w = u + v
//
// History:   JS            28 Jan 2002, 31 Dec 2002
//            Julian Center 28 Aug 2003 debug
//-----------------------------------------------------------------------------
// 
void  AddLabel(
coord_t*  w,      //   O  w = u + v               [Ndim]
coord_t*  u,      // I    can be overwritten      [Ndim]
coord_t*  v,      // I    must not be overwritten [Ndim]
int       Ndim)   // I    dimension
{ 
    int  carry = 0;
    int  i;
    for( i = Ndim-1; i >= 0; i-- )
    {
        w[i] = u[i] + v[i];
        carry = carry ? (++w[i] <= v[i]) : (w[i] < v[i]);
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SubLabel
//
// Purpose:   Set     w = u - v
//
// History:   JS            28 Jan 2002, 31 Dec 2002
//            Julian Center 28 Aug 2003 debug
//-----------------------------------------------------------------------------
// 
void  SubLabel(
coord_t*  w,      //   O  w = u - v               [Ndim]
coord_t*  u,      // I    can be overwritten      [Ndim]
coord_t*  v,      // I    must not be overwritten [Ndim]
int       Ndim)   // I    dimension
{
    int  carry = 0;
    int  i;
    for( i = Ndim-1; i >= 0; i-- )
    {
        w[i] = u[i] - v[i];
        carry = carry ? (--w[i] >= u[i]) : (w[i] > u[i]);
    }
}

//=============================================================================
//              Hilbert-curve (a space-filling Peano curve) library
//=============================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Functions: LinetoAxes
//            AxestoLine
//
// Purpose:   Serial Hilbert length  <---->   multidimensional Axes position.
//
//   Space  = n-dimensional hypercube of side R = 2^b
//            Number of cells = N = R^n = 2^(n*b)
//
//   Line   = serial number of cell along Hilbert curve through hypercube
//          = extended integer of n*b bits ranging from 0 to N-1,
//            stored as vector of n unsigned b-bit integers with [0] high.
//
//   Axes   = Geometrical position of cell
//          = n b-bit integers representing coordinates.
//
// Example:   side R = 16, dimension n = 2, number of cells = N = 256.
//            Line = 9, stored in base-16 words as
//                   Line[0] = 0 (high),   Line[1] = 9 (low),
//            corresponds to position (2,3) as in diagram, stored as
//                   Axes[0] = 2,   Axes[1] = 3.
// 
//        |
//     15 |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
//        |    @   @---@   @   @   @---@   @   @   @---@   @   @   @---@   @
//        |    |           |   |           |   |           |   |           |
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |        |   |           |   |           |   |           |   |    
//        |    @---@   @---@---@---@   @---@   @---@   @---@---@---@   @---@
//        |    |                           |   |                           |
//        |    @   @---@---@   @---@---@   @   @   @---@---@   @---@---@   @
//        |    |   |       |   |       |   |   |   |       |   |       |   |
// Axes[1]|    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |            |           |                   |           |        
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |   |       |   |       |   |   |   |       |   |       |   |
//        |    @   @---@---@   @---@---@   @---@   @---@---@   @---@---@   @
//        |    |                                                           |
//        |    @---@   @---@---@   @---@---@   @---@---@   @---@---@   @---@
//        |        |   |       |   |       |   |       |   |       |   |    
//        |    @---@   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |    |           |           |           |           |           |
//        |    @   @---@   @   @---@   @---@   @---@   @---@   @   @---@   @
//        |    |   |   |   |   |   |       |   |       |   |   |   |   |   |
//        |    @---@   @---@   @   @---@---@   @---@---@   @   @---@   @---@
//        |                    |                           |                
//      3 |    5---6   9---@   @   @---@---@   @---@---@   @   @---@   @---@
//        |    |   |   |   |   |   |       |   |       |   |   |   |   |   |
//      2 |    4   7---8   @   @---@   @---@   @---@   @---@   @   @---@   @
//        |    |           |           |           |           |           |
//      1 |    3---2   @---@   @---@   @---@   @---@   @---@   @---@   @---@
//        |        |   |       |   |       |   |       |   |       |   |    
//      0 |    0---1   @---@---@   @---@---@   @---@---@   @---@---@   @--255
//        |
//         -------------------------------------------------------------------
//             0   1   2   3          ---> Axes[0]                         15
//
// Notes: (1) Unit change in Line yields single unit change in Axes position:
//            the Hilbert curve is maximally local.
//        (2) CPU proportional to total number of bits, = b * n.
//
// History:   John Skilling  20 Apr 2001, 11 Jan 2003, 3 Sep 2003
//-----------------------------------------------------------------------------
// 
void LinetoAxes(
coord_t* Axes,    //   O  multidimensional geometrical axes   [n]
coord_t* Line,    // I    linear serial number, stored as     [n] 
int      b,       // I    # bits used in each word
int      n)       // I    dimension
{
    if( n <= 1 )            // trivial case
        *Axes = *Line;
    else
    {
        LinetoTranspose(Axes, Line, b, n);
        TransposetoAxes(Axes,       b, n);
    }
}

void AxestoLine(
coord_t* Line,    //   O  linear serial number, stored as     [n] 
coord_t* Axes,    // I    multidimensional geometrical axes   [n]
int      b,       // I    # bits used in each word
int      n)       // I    dimension
{
    coord_t  store[1024];   // avoid overwriting Axes
    int      i;             // counter

    if( n <= 1 )            // trivial case
        *Line = *Axes;
    else if( n <= 1024 )    // surely the usual case
    {
        for( i = 0; i < n; ++i )
           store[i] = Axes[i];
        AxestoTranspose(      store, b, n);
        TransposetoLine(Line, store, b, n);
    }
    else                    // must do in place at greater cost
    {
        AxestoTranspose(      Axes, b, n);
        TransposetoLine(Line, Axes, b, n);
        TransposetoAxes(      Axes, b, n);
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Functions: LinetoTranspose
//            TransposetoLine
//
// Purpose:   Recover Hilbert integer by bit-transposition
//
// Example:   b=5 bits for each of n=3 coordinates
//               15-bit Hilbert integer = A B C D E a b c d e 1 2 3 4 5
//                                        X[0]..... X[1]..... X[2].....
//            transposed to
//               X[0](high) = A D b e 3
//               X[1]       = B E c 1 4
//               X[2](low)  = C a d 2 5
//                            high  low
//
// History:   John Skilling  20 Apr 2001, 3 Sep 2003, 14 Oct 2003
//-----------------------------------------------------------------------------
// 
static void LinetoTranspose(
coord_t* X,            //   O  Transpose        [n]
coord_t* Line,         // I    Hilbert integer  [n] 
int      b,            // I    # bits
int      n)            // I    dimension
{
    coord_t  j, p, M;
    int      i, q;

    M = 1 << (b - 1);
    for( i = 0; i < n; i++ )
        X[i] = 0;
    q = 0;    p = M;
    for( i = 0; i < n; i++ )
    {
        for( j = M; j; j >>= 1 )
        {
            if( Line[i] & j )
                X[q] |= p;
            if( ++q == n )
            {
                q = 0;    p >>= 1;
            }
        }
    }
}

static void TransposetoLine(
coord_t* Line,         //   O  Hilbert integer  [n] 
coord_t* X,            // I    Transpose        [n]
int      b,            // I    # bits
int      n)            // I    dimension
{
    coord_t  j, p, M;
    int      i, q;

    M = 1 << (b - 1);
    q = 0;    p = M;
    for( i = 0; i < n; i++ )
    {
        Line[i] = 0;
        for( j = M; j; j >>= 1 )
        {
            if( X[q] & p )
                Line[i] |= j;
            if( ++q == n )
            {
                q = 0;    p >>= 1;
            }
        }
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Functions: TransposetoAxes
//            AxestoTranspose
//
// Purpose:   Transform between Hilbert transpose and geometrical axes
//
// Example:   b=5 bits for each of n=3 coordinates
//            Hilbert transpose
//             X[0] = A D b e 3                  X[1]|  
//             X[1] = B E c 1 4    <------->         |  /X[2]
//             X[2] = C a d 2 5                axes  | /
//                    high  low                      |/______
//                                                         X[0]
//            Axes are stored conventially as b-bit integers.
//         
// History:   John Skilling  20 Apr 2001, 3 Sep 2003, 14 Oct 2003
//-----------------------------------------------------------------------------
// 
static void TransposetoAxes(
coord_t* X,            // I O  position   [n]
int      b,            // I    # bits
int      n)            // I    dimension
{
    coord_t  M, P, Q, t;
    int      i;

// Gray decode by  H ^ (H/2)
    t = X[n-1] >> 1;
    for( i = n-1; i; i-- )
        X[i] ^= X[i-1];
    X[0] ^= t;

// Undo excess work
    M = 2 << (b - 1);
    for( Q = 2; Q != M; Q <<= 1 )
    {
        P = Q - 1;
        for( i = n-1; i; i-- )
            if( X[i] & Q ) X[0] ^= P;                              // invert
            else{ t = (X[0] ^ X[i]) & P;  X[0] ^= t;  X[i] ^= t; } // exchange
        if( X[0] & Q ) X[0] ^= P;                                  // invert
    }
} 
static void AxestoTranspose(
coord_t* X,            // I O  position   [n]
int      b,            // I    # bits
int      n)            // I    dimension
{
    coord_t  P, Q, t;
    int      i;

// Inverse undo
    for( Q = 1 << (b - 1); Q > 1; Q >>= 1 )
    {
        P = Q - 1;
        if( X[0] & Q ) X[0] ^= P;                                  // invert
        for( i = 1; i < n; i++ )
            if( X[i] & Q ) X[0] ^= P;                              // invert
            else{ t = (X[0] ^ X[i]) & P;  X[0] ^= t;  X[i] ^= t; } // exchange
    }

// Gray encode (inverse of decode)
    for( i = 1; i < n; i++ )
        X[i] ^= X[i-1];
    t = X[n-1];
    for( i = 1; i < b; i <<= 1 )
        X[n-1] ^= X[n-1] >> i;
    t ^= X[n-1];
    for( i = n-2; i >= 0; i-- )
        X[i] ^= t;
}

//=============================================================================
//                        Linked-list library
//=============================================================================
//                                    atom___________________________________
//                                   Base|     |     | Ptr | Ptr | Ptr | Ptr |
//                                   Free| Ptr | ... | Ptr | Ptr | Ptr | Ptr |
//                                       |_____|_____|_____|_____|_____|_____|
//                                        /           /  /  /  /  /  /  /  /
//                                0      /           /  /  /  /  /  /  /  / 
//                                :     /           /  /  /  /  /  /  /  /
//                                :    /                vacant nodes
//                              ______/
//                             | Node |
//                             |depth3|
//                             |______|
//                             /:     \.
//                            /:       \.
//                           /:         \.
//                          /:           \.
//                         /:             \.
//                        /:               \.
//                       /:                 \.
//                      /:                   ______
//                     /:                   | Node |
//                    /:                    |depth2|
//                   /:                     |______|
//                  /:                       /:   \.
//                 /:                       /:     \.
//              ______                     /:      ______
//             | Node |                   /:      | Node |
//             |depth1|                  /:       |depth1|
//             |______|                 /:        |______|
//              /:   \                 /:          /:   \.
//             /:     \               /:          /:     \.
//            /:       \             /:          /:       \.
//    Base______      ______      ______      ______      ______
//       | Node |    | Node |    | Node |    | Node |    | Node |
//  0----|depth0|----|depth0|----|depth0|----|depth0|----|depth0|----0
//       |______|    |______|    |______|    |______|    |______|
//          :           :           :           :           :
//          :           :           :           :           :
//          :           :           :           :           :
//    atom_____       _____       _____       _____       _____
//       |extra|     |Atom2|     |Atom0|     |Atom3|     |Atom1|
//       |  0  |     |Label|     |Label|     |Label|     |Label|
//       |atom |     | etc.|     | etc.|     | etc.|     | etc.|
//       |_____|     |_____|     |_____|     |_____|     |_____|
//                   <---------- random-order storage --------->
// 
//          0  <=    <  <  <  <  <  <   Label   <  <  <  <  <     
//
//  Inserted atoms have serial numbers  0, 1, 2, ..., N-1  with labels x >= 0,
//  and occupy consecutive storage. To maintain this efficiently under random
//  insertion and deletion demands that they be stored in dynamically changing
//  order.  Hence coupling to their neighbours in x must be by a linked list,
//  and to avoid ambiguity all labels must be different.
//  
//  The binary tree allows insertion of a new atom with known label x,
//  and deletion of a random atom, at a cost that is asymptotically
//  logarithmic O(log N) because of the search.
//  After each insertion, the newly inserted atom is at the end of storage.
//  and after each deletion, the newly deleted atom is just beyond the end.  
//  After each insertion or deletion, the tree is re-balanced to keep its
//  maximum depth close to log2(N).
//
//  There are 2N+1 nodes of the tree, which only occupy consecutive storage
//  when N is at a historic maximum: when N falls away through deletions, 
//  holes appear in memory.  These are tracked and re-used by a consecutive
//  vector of controlling pointers (which is conveniently stored in the
//  unused atoms to save memory). 
//  
//  Memory = sizeof(Atom) + 2 * sizeof(Node),  per atom.
//  
//  CPU = asymptotically O(log N) per insertion or deletion.
//        Tests with a range of practical N give roughly
//  =======================================
//  CPU = 10 N (log10(N) + 5) multiply-adds per (insertion and deletion) cycle.
//  =======================================
//  (Compare   10 N log10(N)  multiply-adds for a single FFT.)
// 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  SetLink
//
// Purpose:   Allocate and initialise memory,
//            using guess for maximum number of atoms needed.
//            
//            See also ResetLink and FreeLink
//
// History:   John Skilling   9 Mar 1996, 6 Feb 1998, 12 Apr 2001, 10 Sep 2001,
//                           22 Jan 2003
//=============================================================================
//                         STRUCTURE OF ATOMS
//    Node = apex of linked list
//        \.
//         \. 
//          \.
//   List    \______
//       0   | Label|  (unused)
//           | flux |  (unused)
//           | base |  Link up to origin base of tree
//           | free |  Point to apex of tree
//           |______|
//       1   | Label|  x                            ./|\.
//           | flux |  extra properties of atom       |
//           | base |  Link up to base of tree        |
//           | free |  (unused)                       |
//           |______|                                 |  Surviving atoms
//           |      |                                 |
//           ........                                 |  in scrambled
//           |______|                                 |
//   Natoms  |      |                                 |  storage order
//           | .... |                                 |
//           |      |                                 |
//           | .... |                                 |
//           |______|                                \|/
//  Natoms+1 | Label|  (unused)
//           | flux |  (unused)
//           | base |  Point to vacant node address
//           | free |  Point to vacant node address
//           |______|
//           |      |
//           ........
//           |______|
//    NMAX   |      |
//           | .... |
//           |      |
//           | .... |
//           |______|  Reallocate memory if Natoms attempts to exceed NMAX
//
//=============================================================================
//                         STRUCTURE OF NODES
//
//  Link node ______
//           | depth|  >= 1
//           |number|  >= 2
//           |parent|  Link up to parent (NULL if apex)
//           |  lo  |  Link down to "next" child
//           |  hi  |  Link down to "previous" child
//           | atom |  Point down "previous" branches to atom
//           |______|
//
//  Base node ______
//           | depth|  0
//           |number|  1
//           |parent|  Link up to parent (NULL if apex of empty tree)
//           |  lo  |  Link across to  "next"  neighbour   (NULL if last)
//           |  hi  |  Link across to "previous" neighbour (NULL if first)
//           | atom |  Link down to atom
//           |______|
//
//-----------------------------------------------------------------------------
// 
int SetLink(       //   O  0, or -ve error code
int       Ndim,    // I    # words per label
Node*     psLink)  // I O  Linked list
{
    int      NMAX  = 100;   // Initial number of allowed insertions of an atom
    pNode    Nodes = NULL;  // Begin tree storage
    pAtom    List  = NULL;  // List of atoms
    int      i;             // Local counter
    int      CALLvalue = 0;

    psLink->atom   = NULL;
    psLink->parent = NULL;
// Allocate storage for tree
    CALLOC( Nodes, NMAX+NMAX+1, Node )
    CALLOC( List, NMAX+1, Atom )
    CALLOC( List->Label, Ndim * (NMAX+1), coord_t )

// Initialise artificial header atom
    List->Base = Nodes;
    List->Free = Nodes;
// Initialise pointers to available nodes
    for( i = 1; i <= NMAX; ++i )
    {
        List[i].Base  = &Nodes[i+i-1];
        List[i].Free  = &Nodes[i+i];
        List[i].Label = List[i-1].Label + Ndim;
    }
// Initialise top of tree
    Nodes->depth   = 0;
    Nodes->number  = 1;
    Nodes->parent  = NULL;
    Nodes->hi      = NULL;
    Nodes->lo      = NULL;
    Nodes->atom    = List;
    psLink->depth   = NMAX;
    psLink->number   = Ndim;
    psLink->atom   = List;
    psLink->parent = Nodes;
    return 0;
Exit:
    FREE( List->Label )
    FREE( List )
    FREE( Nodes )
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FreeLink
//
// Purpose:   Free memory allocated by SetLink.
//
// History:   John Skilling   9 Mar 1996, 20 Apr 2001
//-----------------------------------------------------------------------------
// 
int  FreeLink(             //   O  0 (or -ve debug code)
Node*  psLink)             // I O  Linked list
{ 
    if( psLink )
    {
        if( psLink->atom )    FREE( psLink->atom->Label )
        FREE( psLink->atom )
        FREE( psLink->parent )
        psLink->parent = NULL;
        psLink->atom   = NULL;
        psLink->depth  = 0;
    }
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  ResetLink
//
// Purpose:   Rellocate memory, if initial guess in SetLink was inadequate
//
// History:   John Skilling   9 Mar 1996, 6 Feb 1998, 12 Apr 2001
//-----------------------------------------------------------------------------
// 
static int  ResetLink(             //   O  0, or -ve error code
pNode  psLink)                     // I O  Linked list
{
    int    NMAX    = psLink->depth;     // Initial # allowed atoms
    pAtom  List    = psLink->atom;      // Old list of atoms
    pNode  Nodes   = psLink->parent;    // Old nodes of tree
    int    Ndim    = psLink->number;    // # words per label

    int    MORE    = NMAX + NMAX/2;     // Extended # allowed atoms
    pAtom  ListX   = NULL;              // Replacement list of atoms
    pNode  NodesX  = NULL;              // Replacement nodes of tree

    int    i;                           // Local counter
    int    j;                           // Counter
    int    CALLvalue = 0;

// Allocate new tree
    CALLOC( NodesX, MORE+MORE+1, Node )
    CALLOC( ListX, MORE+1, Atom )
    CALLOC( ListX->Label, Ndim * (MORE+1), coord_t )
    for( i = 1; i <= MORE; ++i )
        ListX[i].Label = ListX[i-1].Label + Ndim;

// Copy across
    for( i = 0; i <= NMAX + NMAX; ++i )
        NodesX[i] = Nodes[i];
    for( j = 0; j < Ndim; j++ )
        ListX[0].Label[j] = List[0].Label[j];
    for( i = 1; i <= NMAX; ++i )
        for( j = 0; j < Ndim; j++ )
            ListX[i].Label[j] = List[i].Label[j];

// Correct all internal pointers
    for( i = 0; i <= NMAX; ++i )
    {
        ListX[i].Base = NodesX + (List[i].Base - Nodes);
        ListX[i].Free = NodesX + (List[i].Free - Nodes);
    }
    for( i = NMAX + 1; i <= MORE; ++i )
    {
        ListX[i].Base = &NodesX[i+i-1];
        ListX[i].Free = &NodesX[i+i];
    }
    for( i = 0; i <= NMAX + NMAX; ++i )
    {
        if( NodesX[i].parent )
            NodesX[i].parent = NodesX + (Nodes[i].parent - Nodes);
        if( NodesX[i].lo )
            NodesX[i].lo = NodesX + (Nodes[i].lo - Nodes);
        if( NodesX[i].hi )
            NodesX[i].hi = NodesX + (Nodes[i].hi - Nodes);
        NodesX[i].atom = ListX + (Nodes[i].atom - List);
    }
    psLink->atom   = ListX;
    psLink->parent = NodesX;
    psLink->depth  = MORE;

// Free old tree
    FREE( List->Label )
    FREE( List )
    FREE( Nodes )
    return 0;
Exit:
    FREE( ListX->Label )
    FREE( ListX )
    FREE( NodesX )
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  InsAtom
//
// Purpose:   Update list to insert new atom, if new label is distinct.
//            Split node to left, forking to new and previous nodes,
//            then re-balance the tree.
//
// Notes:     After InsAtom, the serial number of the inserted atom is Natoms-1
//
// History:   John Skilling   6 Mar 1996,  7 Oct 1996, 15 Dec 1997, 21 Jan 1998
//                            6 Feb 1998, 12 Apr 2001, 10 Sep 2001, 22 Jan 2003
//=============================================================================
//        Original configuration          Final configuration before balancing
//
//         \         /:       /:                 \          /:       /:
//        ____     ____     ____                ____      ____     ____   
//    ___|Node|___|Node|___|Node|___        ___|Node|    |Node|   |Node|___  
//       |____|   |____|   |____|              |____|    |____|   |____|  
//           :        :                            \     /:  \     / :
//            :        :                            \   /:    \   /  :
//           ____     ____                           ____     ____   :
//          |Atom|   |Atom|                         |New2|___|New1|  :
//          | x- |   | x+ |                         |____|   |____|  :
//          |____|   |____|                          :         :     :
//                                                  :         :      :
//                ____                          ____      ____      ____ 
//               |    |                        |Atom|    |    |    |Atom|
//      InsAtom  | x  | > x-                   | x- |    | x  |    | x+ |
//               |____|                        |____|    |____|    |____|
//
//-----------------------------------------------------------------------------
// 
int InsAtom(         //   O  1 = accept,
                     //      0 = reject because same label,
                     //      -ve = error
pAtom    insertion,  // I    New atom to be inserted into list
pNode    psLink)     // I O  Linked list
{
    int       Natoms;     // # inserted atoms
    pAtom     patom;      // Copy insertion atom into list
    pAtom     List;       // List of atoms
    pNode     node;       // Node of tree
    pNode     New1;       // Available node for inclusion in tree
    pNode     New2;       // Available node for inclusion in tree
    int       Ndim;       // # words per label
    int       j;          // Counter
    int       CALLvalue = 0;

    if( psLink->number <= 0 )
        return CALLvalue;
// Ensure memory is available
    if( psLink->depth < 2 )
        return E_TREEDATA;
    if( NumAtoms(psLink) >= psLink->depth )
        CALL( ResetLink(psLink) )
// Goto base node at or just left of insertion
    List = psLink->atom;
    node = List->Free;
    Ndim = psLink->number;
// log(N) loop
    while( node->depth )
        node = (CmpLabel(insertion->Label, node->hi->atom->Label, Ndim) < 0)
              ? node->lo : node->hi;
    if( CmpLabel(insertion->Label, node->atom->Label, Ndim) == 0 )
    {
        if( CmpLabel(insertion->Label, psLink->atom->Label, Ndim) > 0 )
            return 0;               // Avoid overlapping existing user atom
        if( node->lo )              // but allow one user atom at Label=0
            return 0;
    }
    Natoms = NumAtoms(psLink) + 1;

// Copy atom to vacancy at top of list 
    patom = List + Natoms;
    New1  = patom->Free;
    New2  = patom->Base;
    for( j = 0; j < psLink->number; j++ )
        patom->Label[j] = insertion->Label[j];

// Reconnect base with additional node and twig
    patom->Free  = NULL;
    patom->Base  = New1;
    New1->atom   = patom;
    New1->parent = node;
    New1->hi     = node->hi;
    if( New1->hi )
        New1->hi->lo = New1;
    New1->lo     = New2;
    New1->depth  = 0;
    New1->number = 1;
    node->hi = New1;
    node->atom->Base = New2;
    New2->atom   = node->atom;
    New2->parent = node;
    New2->hi     = New1;
    New2->lo     = node->lo;
    if( New2->lo )
        New2->lo->hi = New2;
    New2->depth  = 0;
    New2->number = 1;
    node->lo = New2;

// Update depth information back up tree, and equilibrate accordingly
    Balance(New1);
    return 1;                        // accepted
Exit:
    return CALLvalue;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  DelAtom
//
// Purpose:   Update list to delete an atom.
//            Cut out specified node and its parent, then re-balance the tree.
//
// History:   John Skilling   6 Mar 1996, 6 Feb 1998, 22 Jan 2003, 2 Dec 2003
//=============================================================================
//        Original configuration          Final configuration before balancing
//  (sibling can be to left or right)      (sibling can be same as lo or hi)
//
//            /                                      \.
//          ____                                      \.
//         | del|                                      \.
//         |____|                                       \.
//          /  \.                                        \.
//         /    \.                                        \.
//       ____    \.                                      ____ 
//      |sibl|    \.                                    |sibl|
//      |____|     \.                                   |____|
//      /    \      \.                                 /      \.
//                   \.    
//          \         \.                                   \           /
//          ____     ____     ____                         ____     ____ 
//      ___| lo |___|kill|___| hi |___                 ___| lo |___| hi |___
//         |____|   |____|   |____|                       |____|   |____|
//           :        :        :                            :        :
//           :        :        :                            :        :
//          ____     ____     ____                         ____     ____ 
//         |Atom|   |    |   |Atom|                       |Atom|   |Atom|
//         | x- |   | x  |   | x+ |                       | x- |   | x+ |
//         |____|   |____|   |____|                       |____|   |____|
//                 deletion
//-----------------------------------------------------------------------------
// 
int DelAtom(              //   O  0 (or -ve error code)
pAtom     deletion,       // I    New atom to be deleted from list
pNode     psLink)         // I O  Linked list
{
    pAtom    xatom;            // Interacting atom
    pAtom    List;             // List of atoms
    pAtom    Top;              // Top of list of atoms
    pNode    node;             // Current node
    pNode    kill;             // Node on tree base, to be killed
    pNode    del;              // kill->parent, also to be killed
    pNode    sibling;          // del->otherchild
    int      Natoms;           // Decremented # atoms
    int      j;                // Counter

    if( psLink->number <= 0 )
        return 0;
    Natoms = NumAtoms(psLink);
// Ensure deletion is from list
    List = psLink->atom;
    Top = List + Natoms;
    if( deletion <= List || deletion > Top )
        return E_TREEDATA;

// Cut node out of base of tree
    kill = deletion->Base;
    if( kill->hi )
        kill->hi->lo = kill->lo;
    kill->lo->hi = kill->hi;

// Cut kill and del out of the tree
    del     = kill->parent;
    sibling = ( del->hi == kill ) ? del->lo : del->hi;
    xatom   = sibling->atom;
    sibling->parent = del->parent;
    if( del->parent )
    {
        if( del->parent->hi == del )
            del->parent->hi = sibling;
        else
            del->parent->lo = sibling;
// Replace all references to killed atom with refs to its sibling
        for( node = del; node->atom == deletion; node = node->parent )
            node->atom = xatom;
    }
    else
    {
        List->Free = sibling;
    }

// Copy atom at top of list into freed slot to keep list storage consecutive
    if( deletion < Top )
    {
        xatom = List + Natoms;
        for( node = xatom->Base; node->atom == xatom; node = node->parent )
            node->atom = deletion;
        deletion->Base = xatom->Base;
        deletion->Free = xatom->Free;
        for( j = 0; j < psLink->number; ++j )
            deletion->Label[j] = xatom->Label[j];
    }

// Collect addresses of two freed leaves
    List[Natoms].Base = del;
    List[Natoms].Free = kill;

// Update depth information back up tree, and equilibrate accordingly
    Balance(sibling);
    return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FindOrder
//
// Purpose:   Find n'th atom in ordered list, or NULL if out of range.
//            Inverse of OrderNum.
//
// History:   John Skilling   10 Sep 2001
//-----------------------------------------------------------------------------
// 
pAtom FindOrder(          //   O  & required atom
      int       n,        // I    order number (0,1,...,NumAtoms-1)
const Node*     psLink)   // I    Linked list
{
    pNode  node;
    if( n < 0 || n >= NumAtoms(psLink) )
        return NULL;
    n++;
    node = psLink->atom->Free;
    while( node->depth )
    {
        if( n >= node->lo->number )
        {
            n -= node->lo->number;
            node = node->hi;
        }
        else
            node = node->lo;
    }
    return node->atom;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FindLeft
//
// Purpose:   Find atom at or left of label.
//            If absent,  return NULL.
//            If present, return right-hand such atom.
//
// History:   John Skilling   6 Mar 1996, 30 Dec 1997, 21 Jan 1998, 12 Apr 2001
//-----------------------------------------------------------------------------
// 
pAtom FindLeft(            //   O  & required atom
      coord_t*  Label,     // I    Label limit
const Node*     psLink)    // I    Linked list
{
    int    Ndim = psLink->number;
    pNode  node = psLink->atom->Free;
    if( Ndim <= 0 )
        return NULL;
    while( node->depth )
        node = ( CmpLabel(Label, node->hi->atom->Label, Ndim) < 0 )
              ? node->lo : node->hi;
    return node->lo ? node->atom : NULL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FindRight
//
// Purpose:   Find atom at or right of label.
//            If absent,  return NULL.
//            If present, return left-hand such atom.
//
// History:   John Skilling   12 Apr 2001
//-----------------------------------------------------------------------------
// 
pAtom FindRight(           //   O  & required atom
      coord_t*  Label,     // I    Label limit
const Node*     psLink)    // I    Linked list
{
    int    Ndim = psLink->number;
    pNode  node = psLink->atom->Free;
    if( Ndim <= 0 )
        return NULL;
    while( node->depth )
        node = ( CmpLabel(Label, node->hi->atom->Label, Ndim) <= 0 )
              ? node->lo : node->hi;
    return node->hi ? node->hi->atom : NULL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  FindHere
//
// Purpose:   Find atom exactly at label.
//            If absent,  return NULL.
//            If present, return address.
//
// History:   John Skilling   2 Dec 2001
//-----------------------------------------------------------------------------
// 
pAtom FindHere(            //   O  & required atom
      coord_t*  Label,     // I    Label limit
const Node*     psLink)    // I    Linked list
{
    int    Ndim = psLink->number;
    pNode  node = psLink->atom->Free;
    if( Ndim <= 0 )
        return NULL;
    while( node->depth )
        node = ( CmpLabel(Label, node->hi->atom->Label, Ndim) < 0 )
              ? node->lo : node->hi;
    if( node->lo && CmpLabel(Label, node->atom->Label, Ndim) == 0 )
        return node->atom;
    else
        return NULL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  EndLink
//
// Purpose:   Find right-most atom.
//
// History:   John Skilling   12 Apr 2001
//-----------------------------------------------------------------------------
// 
pAtom EndLink(             //   O  & atom in required label range
const Node*    psLink)     // I    Linked list
{
    pNode  node;
    if( NumAtoms(psLink) )
    {
        node = psLink->atom->Free;
        while( node->depth )
            node = node->hi;
        return node->atom;
    }
    else
        return NULL;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Storage
//
// Purpose:   Locate atom in current storage list.  Inverse of FindAtom.
//
// History:   John Skilling   22 Jan 2003
//-----------------------------------------------------------------------------
// 
int Storage(    //   O  location in current storage list 0,,..,NumAtoms-1
pAtom   atom)   // I    & atom in list
{
    pNode  node;
    for( node = atom->Base; node->parent; node = node->parent );
    return  atom - node->atom - 1;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  OrderNum
//
// Purpose:   Order number of atom in linked list.  Inverse of FindOrder.
//
// History:   John Skilling   22 Jan 2003
//-----------------------------------------------------------------------------
// 
int OrderNum(   //   O  location in ordered list 0,1,..,NumAtoms-1
pAtom   atom)   // I    & atom in list
{
    pNode  node;
    int    number = -1;
    for( node = atom->Base; node->parent; node = node->parent )
        if( node == node->parent->hi )
            number += node->parent->lo->number;
    return  number;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Function:  Balance
//
// Purpose:   Re-balance the tree depths, after insertion or deletion.
//
// History:   John Skilling   25 Jan 1996, 10 Sep 2001
//=============================================================================
//           Original configurations      Final configurations
//  
//  Case A               /                        /
//                     Node                     Node
//                    /    \                    /    \.
//                   /     Node              Node     \.
//                  /     /    \            /    \     \.
//  
//
//  Case B               /                        /
//                     Node                     Node
//                    /    \                   /    \.
//                   /     Node             Node    Node
//                  /     /   \             /  \    /  \.
//                 /   Node    \. 
//                /    /   \    \.
//
//  
//  Case C               /                        /             
//                     Node                     Node
//                    /    \                   /    \.
//                 Node     \                 /     Node
//                /    \     \               /     /    \.
//  
//  
//  Case D               /                        /
//                     Node                     Node
//                    /    \                   /    \.
//                 Node     \               Node    Node
//                 /   \     \              /  \    /  \.
//                /    Node   \. 
//               /    /   \    \.
//  
//-----------------------------------------------------------------------------
// 
static void Balance(
pNode   node)          // I    Base node of revised tree
{
#define DEPTH(x,y)    (((x) > (y) ? (x) : (y)) + 1)
    pNode t;                // Temporary pointer for swapping
    int   asymmetry;        // Depth difference between branches
    int   deep;             // Recalculated depth of node

    for( node = node->parent; node; node = node->parent )
    {
        asymmetry = node->hi->depth - node->lo->depth;
        if( asymmetry > 1 )
        {
            if( node->hi->hi->depth >= node->hi->lo->depth )
            {
// Case A
                t         = node->hi;
                node->hi  = t->hi;       node->hi->parent = node;
                t->hi     = t->lo;
                t->lo     = node->lo;    t->lo->parent    = t;
                node->lo  = t;
                t->atom   = node->atom;
                t->depth  = DEPTH(t->hi->depth, t->lo->depth);
                t->number = t->hi->number + t->lo->number;
            }
            else
            {
// Case B
                t = node->hi->lo;
                node->hi->atom = t->hi->atom;
                t->atom = node->atom;
                node->hi->lo = t->hi;       t->hi->parent = node->hi;
                t->hi = t->lo;
                t->lo = node->lo;           t->lo->parent = t;
                node->lo = t;               t->parent = node;
                t->depth  = DEPTH(t->hi->depth, t->lo->depth);
                t->number = t->hi->number + t->lo->number;
                t = node->hi;
                t->depth  = DEPTH(t->hi->depth, t->lo->depth);
                t->number = t->hi->number + t->lo->number;
            }
        }
        if( asymmetry < -1 )
        {
            if( node->lo->lo->depth >= node->lo->hi->depth )
            {
// Case C
                t         = node->lo;
                node->lo  = t->lo;       node->lo->parent = node;
                t->lo     = t->hi;
                t->hi     = node->hi;    t->hi->parent    = t;
                node->hi  = t;
                t->atom   = t->lo->atom;
                t->depth  = DEPTH(t->hi->depth, t->lo->depth);
                t->number = t->hi->number + t->lo->number;
            }
            else
            {
// Case D
                t = node->lo->hi;
                t->atom = t->hi->atom;
                node->lo->hi = t->lo;       t->lo->parent = node->lo;
                t->lo = t->hi;
                t->hi = node->hi;           t->hi->parent = t;
                node->hi  = t;              t->parent = node;
                t->depth  = DEPTH(t->hi->depth, t->lo->depth);
                t->number = t->hi->number + t->lo->number;
                t = node->lo;
                t->depth  = DEPTH(t->hi->depth, t->lo->depth);
                t->number = t->hi->number + t->lo->number;
            }
        }
        deep = node->depth;
        node->depth  = DEPTH(node->hi->depth, node->lo->depth);
        node->number = node->hi->number + node->lo->number;
        if( node->depth == deep )
            break;
    }
    for( ; node; node = node->parent )
        node->number = node->hi->number + node->lo->number;
}
