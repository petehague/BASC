//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Filename:  hilbert.h
// 
// Purpose:   Header for hilbert.c
// 
// History:   Incorporates Hilbert and Linked-list codes.
//            John Skilling   28 Jan 2002 - 11 Oct 2003
//-----------------------------------------------------------------------------
//
#ifndef HILBERTH
#define HILBERTH

typedef unsigned int coord_t; // char,short,int for 8,16,32 bits per word

struct Node_;
struct Atom_;

typedef struct Atom_
{
    coord_t*       Label;   // label as composite integer ([0] high)
    struct Node_*  Base;
    struct Node_*  Free;
} Atom;                     // Element of linked list

typedef struct Node_
{
    int            depth;
    int            number;
    struct Node_*  parent;
    struct Node_*  lo;
    struct Node_*  hi;
    struct Atom_*  atom;
} Node;                     // Linkage in linked list

#ifdef __cplusplus
extern "C" {
#endif
// Composite-integer arithmetic ...............................................

extern int CmpLabel(        // Output sign(u-v) for comparison
       coord_t*  u,         // Input composite integer ([0] high)
       coord_t*  v,         // Input composite integer ([0] high)
       int       Ndim);     // Dimension

extern void AddLabel(       // Add labels
       coord_t*  w,         // Output sum u+v
       coord_t*  u,         // Input composite integer
       coord_t*  v,         // Input composite integer
       int       Ndim);     // Dimension

extern void SubLabel(       // Subtract labels
       coord_t*  w,         // Output difference u-v
       coord_t*  u,         // Input composite integer
       coord_t*  v,         // Input composite integer
       int       Ndim);     // Dimension

// Hilbert curve ..............................................................

extern void LinetoAxes(     // multidimensional posn <--- serial Hilbert length
       coord_t*  Axes,      //   O  multidimensional position         [Ndim]
       coord_t*  Line,      // I    linear serial number, stored as   [Ndim] 
       int       Nbits,     // I    # bits used for each coordinate
       int       Ndim);     // I    dimension

extern void AxestoLine(     // serial Hilbert length <--- multidimensional posn
       coord_t*  Line,      //   O  linear serial number, stored as   [Ndim] 
       coord_t*  Axes,      // I    multidimensional position         [Ndim]
       int       Nbits,     // I    # bits used for each coordinate
       int       Ndim);     // I    dimension

// Linked list ................................................................

extern int SetLink(         // Setup Link
       int       Ndim,      // Dimension
       Node*     Link);     // Linked list

extern int FreeLink(        // Free memory allocated by SetLink
       Node*     Link);     // Linked list

extern int InsAtom(         // Insert atom (1=accept, 0=reject, -ve=error)
       Atom*     insertion, // &(atom to be inserted)
       Node*     Link);     // Linked list

extern int DelAtom(         // Delete atom
       Atom*     deletion,  // &(atom to be deleted)
       Node*     Link);     // Linked list

extern Atom* EndLink(       // Output &(last atom in list): see BeginLink
 const Node*     Link);     // Linked list

extern Atom* FindHere(      // Output &(atom at Label), or NULL
       coord_t*  Label,     // Input composite integer
 const Node*     Link);     // Linked list

extern Atom* FindLeft(      // Output &(atom at or before Label), or NULL
       coord_t*  Label,     // Input composite integer
 const Node*     Link);     // Linked list

extern Atom* FindRight(     // Output &(atom at or after Label), or NULL
       coord_t*  Label,     // Input composite integer
 const Node*     Link);     // Linked list

extern Atom* FindOrder(     // Output &(n'th atom in ordered list)
       int       n,         // order number in list
 const Node*     Link);     // Linked list

extern int OrderNum(        // Output order number in list: see FindOrder
       Atom*     atom);     // Input &(atom)

extern int Storage(         // Output position in current storage: see FindAtom
       Atom*     atom);     // Input &(atom)

#undef  FindAtom            // Output &(n'th atom in current storage of Link)
#define FindAtom(n,Link)     (NumAtoms(Link) ? (Link)->atom + (n) + 1 : NULL)

#undef  BeginLink           // Output &(1st atom in list of Link): see EndLink 
#define BeginLink(Link)      (NumAtoms(Link)?(Link)->atom->Base->hi->atom:NULL)

#undef  NextAtom            // Output &(rightward neighbour atom to A), or NULL
#define NextAtom(A)          ((A)&&(A)->Base->hi ? (A)->Base->hi->atom : NULL)

#undef  PrevAtom            // Output &(leftward neighbour atom to A), or NULL
#define PrevAtom(A)          ((A)&&(A)->Base->lo->lo?(A)->Base->lo->atom:NULL)

#undef  NumAtoms            // Output number of atoms in Link
#define NumAtoms(Link)       ((Link)->atom->Free->number - 1)

//.............................................................................
#ifdef __cplusplus
};
#endif

#undef  E_TREEDATA
#define E_TREEDATA        -261  // Wrong input sizes
#undef  E_MALLOC
#define E_MALLOC          -130  // Can't allocate memory

#endif

