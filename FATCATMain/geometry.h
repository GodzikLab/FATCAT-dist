//---------------------------------------------------------------------------
//This template includes some common geometry related functions,
//eg. distance of two points, angle of three points, dihedral of three points et al
//Y.Y 10/22/02 modified from MyTool.C (MyTool.h)
//---------------------------------------------------------------------------

#ifndef MYGEOMETRY_H
#define MYGEOMETRY_H

#include "basic.h"
#define	PI 3.1415926
#define pow1(x) ((x) * (x))

namespace GEOMETRY	{
	//---------------------------------------------------------------------------
	//return the distance between two points x and y
	//---------------------------------------------------------------------------
	double distance(double *x, double *y);
	//---------------------------------------------------------------------------
	// translates the point X0 to LOC_COORD system , then rotates it
	//----------------------------------------------------------------------------
	void tran_ord(double rotmtx[][3],double *orig,double *x0);
	void tran_ord(double *rotmtx,double *orig,double *x0);
	//---------------------------------------------------------------------------
	//superimpose two sets of coordinations: c1 and c2
	//record the transformation information in r and t
	//return rmsd
	//----------------------------------------------------------------------------
	double kearsay(int ln, double *c1, double *c2, double *r, double *t);
	//---------------------------------------------------------------------------
	// The program carries on the Jacobi translation
	//---------------------------------------------------------------------------
	int jacobi(double mat[][4], int n, double *lamda, double vec[][4]);
	//---------------------------------------------------------------------------
	// The program return the root mean square of two distance matrix 
	//---------------------------------------------------------------------------
	double dismatrixrms(int ln, double *dis1, double *dis2);
	//---------------------------------------------------------------------------
	void srteig(double *lamda, int n, double vec[][4]);
	void rot2euler(double *rot, double *euler_angle);
};

#endif
