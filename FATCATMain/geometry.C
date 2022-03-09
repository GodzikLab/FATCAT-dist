//---------------------------------------------------------------------------
//This template includes some common geometry related functions,
//eg. distance of two points, angle of three points, dihedral of three points et al
//Y.Y 10/22/02 modified from MyTool.C (MyTool.h)
//---------------------------------------------------------------------------

#include "geometry.h"

//---------------------------------------------------------------------------
//return the distance between two points x and y
//---------------------------------------------------------------------------
double GEOMETRY::distance(double *x, double *y)
{
	int    	i;
	double  v;
	double	v0;
        v = 0.0;
	for (i = 0; i < 3; i++) {
		v0 = x[i] - y[i];
	        v += v0 * v0;
	}
	return (sqrt(v));
}

//---------------------------------------------------------------------------
// translates the point X0 to LOC_COORD system , then rotates it
//----------------------------------------------------------------------------
void GEOMETRY::tran_ord(double rotmtx[][3],double *orig,double *x0)
{
	int	i;
	double	x[3];

	for (i=0;i<3;i++)
		x[i]=x0[i]-orig[i];
	for (i=0;i<3;i++)
		x0[i]=rotmtx[i][0]*x[0]+rotmtx[i][1]*x[1]+rotmtx[i][2]*x[2];
        return;
}

void GEOMETRY::tran_ord(double *rotmtx,double *orig,double *x0)
{
	int	i;
	double	x[3];

	for (i=0;i<3;i++)
		x[i]=x0[i]-orig[i];
	for (i=0;i<3;i++)
		x0[i]=rotmtx[3 * i]*x[0]+rotmtx[3 * i + 1]*x[1]+rotmtx[3 * i + 2]*x[2];
        return;
}

//---------------------------------------------------------------------------
//superimpose two sets of coordinations: c1 and c2
//record the transformation information in r and t
//return rmsd
//----------------------------------------------------------------------------
double GEOMETRY::kearsay(int ln, double *c1, double *c2, double *r, double *t)
{
	int           i, j, k;
	double        cen1[3], cen2[3], mat[4][4], temp;
	double        vec[4][4], lamda[4], d;
        double        *xm, *ym, *zm, *xp, *yp, *zp;

        xm = new double[ln];
        ym = new double[ln];
        zm = new double[ln];
        xp = new double[ln];
        yp = new double[ln];
        zp = new double[ln];

	for(j = 0; j < 3; j ++)
		cen1[j] = cen2[j] = 0.0;
        for(i = 0; i < 4; i ++) {
                lamda[i] = 0.0;
                for(j = 0; j < 4; j ++) {
                        mat[i][j] = vec[i][j] = 0.0;
                }
        }
	for(i = 0; i < ln; i ++)     {
		for(j = 0; j < 3; j ++)       {
			cen1[j] += c1[3 * i + j];
			cen2[j] += c2[3 * i + j];
		}
	}
	for(j = 0; j < 3; j ++)       {
		cen1[j] /= ln;
		cen2[j] /= ln;
	}

	for(i = 0; i < ln; i ++)     {
		xm[i] = (c1[3 * i] - cen1[0]) - (c2[3 * i] - cen2[0]);
		ym[i] = (c1[3 * i + 1] - cen1[1]) - (c2[3 * i + 1] - cen2[1]);
		zm[i] = (c1[3 * i + 2] - cen1[2]) - (c2[3 * i + 2] - cen2[2]);
		xp[i] = (c1[3 * i] - cen1[0]) + (c2[3 * i] - cen2[0]);
		yp[i] = (c1[3 * i + 1] - cen1[1]) + (c2[3 * i + 1] - cen2[1]);
		zp[i] = (c1[3 * i + 2] - cen1[2]) + (c2[3 * i + 2] - cen2[2]);
	}
	for(i = 0; i < ln; i ++)     {
		mat[0][0] += pow1(xm[i]) + pow1(ym[i]) + pow1(zm[i]);
		mat[1][1] += pow1(yp[i]) + pow1(zp[i]) + pow1(xm[i]);
		mat[2][2] += pow1(xp[i]) + pow1(zp[i]) + pow1(ym[i]);
		mat[3][3] += pow1(xp[i]) + pow1(yp[i]) + pow1(zm[i]);
		mat[0][1] += yp[i] * zm[i] - ym[i] * zp[i];
		mat[0][2] += xm[i] * zp[i] - xp[i] * zm[i];
		mat[0][3] += xp[i] * ym[i] - xm[i] * yp[i];
		mat[1][2] += xm[i] * ym[i] - xp[i] * yp[i];
		mat[1][3] += xm[i] * zm[i] - xp[i] * zp[i];
		mat[2][3] += ym[i] * zm[i] - yp[i] * zp[i];
	}
	for(i = 0; i < 3; i ++)
		for(j = i + 1; j < 4; j ++)  {
			mat[j][i] = mat[i][j];
		}
	jacobi(mat, 4, lamda, vec);
	srteig(lamda, 4, vec);
	d = lamda[3] / ln;
	d = sqrt(fabs(d));
	r[0] = pow1(vec[0][3]) + pow1(vec[1][3]) - pow1(vec[2][3])
			- pow1(vec[3][3]);
	r[1] = 2.0 * (vec[1][3] * vec[2][3] + vec[0][3] * vec[3][3]);
	r[2] = 2.0 * (vec[1][3] * vec[3][3] - vec[0][3] * vec[2][3]);
	r[3] = 2.0 * (vec[1][3] * vec[2][3] - vec[0][3] * vec[3][3]);
	r[4] = pow1(vec[0][3]) - pow1(vec[1][3]) + pow1(vec[2][3])
			- pow1(vec[3][3]);
	r[5] = 2.0 * (vec[2][3] * vec[3][3] + vec[0][3] * vec[1][3]);
	r[6] = 2.0 * (vec[1][3] * vec[3][3] + vec[0][3] * vec[2][3]);
	r[7] = 2.0 * (vec[2][3] * vec[3][3] - vec[0][3] * vec[1][3]);
	r[8] = pow1(vec[0][3]) - pow1(vec[1][3]) - pow1(vec[2][3])
			+ pow1(vec[3][3]);
	for(j = 0; j < 3; j ++)    {
		temp = 0.0;
		for(k = 0; k < 3; k ++)
		    temp += cen1[k] * r[3 * k + j];
		t[j] = cen2[j] - temp;
	}

        delete[] xm;
        xm = NULL;
        delete[] ym;
        ym = NULL;
        delete[] zm;
        zm = NULL;
        delete[] xp;
        xp = NULL;
        delete[] yp;
        yp = NULL;
        delete[] zp;
        zp = NULL;

	return d ;
}

// The program carries on the Jacobi translation
//---------------------------------------------------------------------------
int GEOMETRY::jacobi(double mat[][4], int n, double *lamda, double vec[][4])
{
	int         i, j, k, m, iter;
	double      tresh, tau, sm, g, C, H, S, T, theta, *B, *Z;

        B = new double[n];
        Z = new double[n];
	for(i = 0; i < n; i ++)   {
		for(j = 0; j < n; j ++)
			vec[i][j] = 0.0;
		vec[i][i] = 1.0;
	}
	for(i = 0; i < n; i ++)    {
		lamda[i] = B[i] = mat[i][i];
		Z[i] = 0.0;
	}
	iter = 0;
	for(i = 0; i < 50; i ++)     {
		sm = 0.0;
		for(j = 0; j < n - 1; j ++)
			for(k = j + 1; k < n; k ++)
				 sm += fabs(mat[j][k]);
		if(sm == 0.0)     {
			delete[] B;
                        B = NULL;
                        delete[] Z;
                        Z = NULL;
			return iter;
		}
		if(i < 3)
			tresh = 0.2 * sm / pow1(n);
		else
			tresh = 0.0;
		for(j = 0; j < n - 1; j ++)    {
			for(k = j + 1; k < n; k ++)   {
				g = 100.0 * fabs(mat[j][k]);
				if((i > 3) &&
				      (fabs(lamda[j]) + g == fabs(lamda[j])) &&
				      (fabs(lamda[j]) + g == fabs(lamda[k])))
						mat[j][k] = 0.0;
				else if(fabs(mat[j][k]) > tresh)  {
					H = lamda[k] - lamda[j];
					if(fabs(H) + g == fabs(H))
						T = mat[j][k] / H;
					else   {
						theta = 0.5 * H / mat[j][k];
						T = 1.0 / (fabs(theta) +
						   sqrt(1.0 + pow1(theta)));
						if(theta < 0.0)  T = -T;
					}
					C = 1.0 / sqrt(1.0 + pow1(T));
					S = T * C;
					tau = S / (1.0 + C);
					H = T * mat[j][k];
					Z[j] -= H;
					Z[k] += H;
					lamda[j] -= H;
					lamda[k] += H;
					mat[j][k] = 0.0;
					for(m = 0; m < j; m ++)   {
						g = mat[m][j];
						H = mat[m][k];
						mat[m][j]=g-S*(H+g*tau);
						mat[m][k]=H+S*(g-H*tau);
					}
					for(m = j + 1; m < k; m ++)   {
						g = mat[j][m];
						H = mat[m][k];
						mat[j][m]=g-S*(H+g*tau);
						mat[m][k]=H+S*(g-H*tau);
					}
					for(m = k + 1; m < n; m ++)   {
						g = mat[j][m];
						H = mat[k][m];
						mat[j][m]=g-S*(H+g*tau);
						mat[k][m]=H+S*(g-H*tau);
					}
					for(m = 0; m < n; m ++)    {
						g = vec[m][j];
						H = vec[m][k];
						vec[m][j]=g-S*(H+g*tau);
						vec[m][k]=H+S*(g-H*tau);
					}
					iter ++;
				}
			}
		}
		for(j = 0; j < n; j ++)   {
			B[j] += Z[j];
			lamda[j] = B[j];
			Z[j] = 0.0;
		}
	}
	cout<<"50 iteration never happen!"<<endl;

        delete[] B;
        B = NULL;
        delete[] Z;
        Z = NULL;
	return iter;
}

//---------------------------------------------------------------------------
double GEOMETRY::dismatrixrms(int ln, double *dis1, double *dis2)
{
	double	rms = 0;
	for(int i = 0; i < ln; i ++)	{
		rms += (dis1[i] - dis2[i]) * (dis1[i] - dis2[i]);
	}
	return (sqrt(rms / ln));
}

//---------------------------------------------------------------------------
void GEOMETRY::srteig(double *lamda, int n, double vec[][4])
{
	int           i, j, k;
	double        p;

	for(i = 0; i < n - 1; i ++)     {
		k = i;
		p = lamda[i];
		for(j = i + 1; j < n; j ++)    {
			if(lamda[j] >= p)   {
				k = j;
				p = lamda[j];
			}
		}
		if(k != i)   {
			lamda[k] = lamda[i];
			lamda[i] = p;
			for(j = 0; j < n; j ++)    {
				p = vec[j][i];
				vec[j][i] = vec[j][k];
				vec[j][k] = p;
			}
		}
	}
}

//---------------------------------------------------------------------------
//rotation-eulerangle transformation
//---------------------------------------------------------------------------
void GEOMETRY::rot2euler(double *rot, double *euler_angle)
{
        int     i;
        double  costheta, tanfai, tanpsi;

        for(i = 0; i < 3; i ++) euler_angle[i] = 0.0;
        costheta = rot[8];
        euler_angle[1] = acos(costheta);

        if(rot[6] == 0.0)       euler_angle[2] = 0.0;
        else if(rot[7] == 0.0)  euler_angle[2] = PI / 2;
        else    {
                tanpsi = rot[6] / rot[7];
                euler_angle[2] = atan(tanpsi);
        }
        if(rot[2] == 0.0)       euler_angle[0] = 0.0;
        else if(rot[5] == 0.0)  euler_angle[0] = -PI / 2;
        else    {
                tanfai = -rot[2] / rot[5];
                euler_angle[0] = atan(tanfai);
        }

        if(sin(euler_angle[1]) * rot[7] < 0)    {
                if(euler_angle[2] > 0)  euler_angle[2] -= PI;
                else                    euler_angle[2] += PI;
        }
        if(sin(euler_angle[1]) * rot[5] > 0)    {
                if(euler_angle[0] > 0)  euler_angle[0] -= PI;
                else                    euler_angle[0] += PI;
        }

        for(i = 0; i < 3; i ++)
                euler_angle[i] = 180.0 * euler_angle[i] / PI;
}

