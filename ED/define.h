#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

#include <vector>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"
#include <bitset>

#define DIM 2
#define PI 2*asin(1.0)
#define tn 1.0
#define Nlink 72
#define Nbit 108
extern int *next[2*DIM+1], *parity;
extern int LX, LY, VOL, GL, iter, FLAG;
extern double Vn;


extern uint64_t NST;
//extern int t,V;
extern int alignment;
extern double *psiBpsi, *Sx, *Sy, *cdw, *plaqOp1;

#endif
