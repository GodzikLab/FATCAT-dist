#include <cstdlib>

#include "matrix.h"

//inverse a matrix (a)
//input matrix a, the row number
//output: the inversed matrix in a
void MATRIX::invmatrix(int nrow, double **a)
{
	int 	i,icol,irow,j,k,l,ll;
	double	big,dum,pivinv;
	double	sw;
	icol = irow = 0;
	int 	n = nrow;
	int	*indxc = new int[n];
	int	*indxr = new int[n];
	int	*ipiv = new int[n];
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<n;l++) {
				sw = a[irow][l];
				a[irow][l] = a[icol][l];
				a[icol][l] = sw;
			}
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) {
			printf("InverseMatrix: Singular Matrix");
			exit(1);
		}
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l]) {
			for (k=0;k<n;k++)	{
				sw = a[k][indxr[l]];
				a[k][indxr[l]] = a[k][indxc[l]];
				a[k][indxc[l]] = sw;
			}
		}
	}
	delete[] indxc;
	indxc = NULL;
	delete[] indxr;
	indxr = NULL;
	delete[] ipiv;
	ipiv = NULL;
}

//return a multiplier of two matrices
double** MATRIX::mulmatrix(int row1, int col1, double **a, int row2, int col2, double **b)
{
	if(col1 != row2)	return NULL;
	double	**c = new double*[row1];
	int	i, j, m;
	double	cell;
	for(i = 0; i < row1; i ++)	c[i] = new double[col2];
	for(i = 0; i < row1; i ++)	{
		for(j = 0; j < col2; j ++)	{
			cell = 0.0;
			for(m = 0; m < col1; m ++)	{
				cell += a[i][m] * b[m][j];
			}
			c[i][j] = cell;
		}	
	}
	return c;
}

//return the trace of a matrix: the sum of the diagonal elements
double	MATRIX::tracematrix(int row, double **a)
{
	double	sum = 0;
	for(int i = 0; i < row; i ++)	{
		sum += a[i][i];
	}
	return sum;
}


void MATRIX::rot(double **a, const double s, const double tau, const int i, const int j, const int k, const int l)
{
	double g,h;

	g=a[i][j];
	h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}

//jacobi function to diagonal the input matrix **a,
//eigenvalue is in *d, eigenvector is in **v
void MATRIX::jacobi(int row, double **a, double *d, double **v, int &nrot)
{
	int 	i, j, ip, iq;
	double 	tresh, theta, tau, t, sm, s, h, g, c;

	int 	n = row;
	double	*b = new double[n];
	double	*z = new double[n];
	for(ip = 0; ip < n; ip ++) {
		for(iq = 0; iq < n; iq ++) v[ip][iq] = 0.0;
		v[ip][ip] = 1.0;
	}
	for(ip = 0; ip < n; ip ++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	nrot = 0;
	for(i = 1; i <= 50; i ++) {
		sm = 0.0;
		for(ip = 0; ip < n - 1; ip ++) {
			for(iq = ip + 1; iq < n; iq ++)
				sm += fabs(a[ip][iq]);
		}
		if(sm == 0.0)	{
			delete[] b;
			delete[] z;
			return;
		}
		if(i < 4)
			tresh = 0.2 * sm / (n * n);
		else
			tresh = 0.0;
		for(ip = 0; ip < n - 1; ip ++) {
			for(iq = ip + 1; iq < n; iq ++) {
				g = 100.0 * fabs(a[ip][iq]);
				if(i > 4 && (fabs(d[ip]) + g) == fabs(d[ip])
					&& (fabs(d[iq]) + g) == fabs(d[iq]))
						a[ip][iq] = 0.0;
				else if(fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if((fabs(h) + g) == fabs(h))
						t = (a[ip][iq])/h;
					else {
						theta = 0.5*h / (a[ip][iq]);
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta*theta));
						if(theta < 0.0) t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for(j = 0; j < ip; j ++)
						rot(a, s, tau, j, ip, j, iq);
					for(j = ip + 1; j < iq; j ++)
						rot(a, s, tau, ip, j, j, iq);
					for(j = iq + 1; j < n; j ++)
						rot(a, s, tau, ip, j, iq, j);
					for(j = 0; j < n; j ++)
						rot(v, s, tau, j, ip, j, iq);
					++ nrot;
				}
			}
		}
		for(ip = 0; ip < n; ip ++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	printf("Too many iterations in routine jacobi");
	delete[] b;
	delete[] z;
	exit(1);
}

//sort the eigenvalues
void MATRIX::eigsrt(int row, double *d, double **v)
{
	int 	i,j,k;
	double 	p;

	int 	n = row;
	for(i = 0; i < n - 1; i ++) {
		p = d[k=i];
		for(j = i; j < n; j ++)
			if (d[j] >= p) p = d[k=j];
		if(k != i) {
			d[k] = d[i];
			d[i] = p;
			for(j = 0;j < n; j ++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}
