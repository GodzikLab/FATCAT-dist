#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "basic.h"

namespace MATRIX {
	void 	invmatrix(int nrow, double **a); 
	//inverse matrix by gauss-jordan elimination (modified from gaussj() in numeric recipe) 
	double** mulmatrix(int row1, int col1, double **a, int row2, int col2, double **b);
	//multiply two a and b, return the multiplier
	double	tracematrix(int nrow, double **a);
	void	rot(double **a, const double s, const double tau, const int i, const int j, const int k, const int l);
	void	jacobi(int row, double **a, double *d, double **v, int &nrot);
	void	eigsrt(int row, double *d, double **v);

}
;
#endif
