#include <cstdlib> 

#include "SAlnOpt.h"
#include "Align0.h"

using namespace ARRAY;
//latest update, 5/1/03, Y.Y

//--------------------------------------------------------------------------------------------------------
//optimize the structural alignment by update the equivalent residues 
//and then run dynamic programming
//input: len1 the length of structure 1; c1: the structure information of 1
//       len2 the length of structure 2; c2: the structure information of 2
//       iniLen and iniSet is the length and list of initial equivalent residues
//       maxi: the maximum iteration
//--------------------------------------------------------------------------------------------------------
SALNOPT::SALNOPT(int len1, double *c1, int len2, double *c2, int iniLen, int **iniSet, int maxi)
{
	//initial structures
	int	i, k;
	pro1Len = len1;
	pro2Len = len2;
	cod1 = NewArray <double> (len1 * 3); 
	cod2 = NewArray <double> (len2 * 3);
	for(i = 0; i < len1; i ++)	{
		for(k = 0; k < 3; k ++)	{
			cod1[3 * i + k] = c1[3 * i + k];
		}
	}
	for(i = 0; i < len2; i ++)	{
		for(k = 0; k < 3; k ++)	{
			cod2[3 * i + k] = c2[3 * i + k];
		}
	}

	//initial equivalent sets
	maxLen = (len1 < len2)?len1:len2;
	equSet = NewMatrix <int> (2, maxLen); 
	for(i = 0; i < iniLen; i ++)	{
		equSet[0][i] = iniSet[0][i];
		equSet[1][i] = iniSet[1][i];
		if(iniSet[0][i] > len1 || iniSet[1][i] > len2)	{
			printf("Warn: focus exceeds the protein 1 or 2 length\n");
			exit(1);
		}
	}	
	equLen = iniLen;
	equLen0 = equLen;

	SetParameters();

	//calloc the matrix
	sij = NewMatrix <double> (pro1Len, pro2Len);

	SuperimposeBySet();
	//printf("   initial rmsd %f\n", rmsd);

	maxKeepStep = 4;
	keepStep = 0;

	Optimize(maxi);
}

//decalloc the class
SALNOPT::~SALNOPT(void)
{
	if(equSet != NULL)	DelMatrix <int> (equSet, 2);	
	if(cod1 != NULL)	DelArray <double> (cod1);
	if(cod2 != NULL)	DelArray <double> (cod2);
	if(sij != NULL)		DelMatrix <double> (sij, pro1Len);
}

//--------------------------------------------------------------------------------------------------------
//refer CE, similarity = Dc - dij, Dc is increased by 0.5 each cycle,
//optimization continues until either
//i)alignment length is less than 95% of alignment length before optimization
//ii)rmsd is less than 110% of rmsd at the cycle when condition i) was first satisfied
//--------------------------------------------------------------------------------------------------------
void SALNOPT::SetParameters(void)
{
	Dc = 3.0; //Dc = 2.0 
	increase = 0.5;
	stopLenPer = 0.95;
	stopRmsdPer = 1.1;
	stopRmsd = -1.0;
	rmsdCut = 3.0;
	gapIni = 5.0;
	gapExt = 0.5;
}

//--------------------------------------------------------------------------------------------------------
void SALNOPT::Optimize(int maxi)
{
	int	i, ifstop, alnLen;
	int	**alnList = NewMatrix <int> (2, maxLen);
	for(i = 0; i < maxi; i ++)	{
		CalMatrix();
		ALIGN0 *aln = new ALIGN0(sij, pro1Len, pro2Len, gapIni, gapExt);
		alnLen = aln->GetAlignPos(alnList);
		if(alnLen < 3)	ifstop = 1; //very rare, mark by Y.Y on 5/1/03
		else	ifstop = DefineEquPos(alnLen, alnList);
		delete aln;
		if(ifstop)	break;
		Dc += increase;
	}
	DelMatrix <int> (alnList, 2);
	/*
	if(i < maxi)	{
		printf("   optimize converged at %d iterations\n", i);
	}
	else	printf("   optimize stop without convergence\n");
	*/
}

//--------------------------------------------------------------------------------------------------------
//the definition of matrix between residues:
//             Sij = Dc^2 - Dij^2 if Dij <= Dc
//                   0            else
//--------------------------------------------------------------------------------------------------------
void SALNOPT::CalMatrix(void)
{
	int	i, j;
	double	dis;
	for(i = 0; i < pro1Len; i ++)	{
		for(j = 0; j < pro2Len; j ++)	{
			dis = GEOMETRY::distance(&cod1[3 * i], &cod2[3 * j]);
			if(dis < Dc)	sij[i][j] = Dc - dis;
			else	sij[i][j] = 0;
		}
	}
}

//--------------------------------------------------------------------------------------------------------
//superimpose two structures according to the equivalent residues
//--------------------------------------------------------------------------------------------------------
void SALNOPT::SuperimposeBySet(void)
{
	//extract the coordinations of equivalent residues
	double *tmp1 = NewArray <double> (3 * equLen);
	double *tmp2 = NewArray <double> (3 * equLen);
	int	i, k, r1, r2, n1, n2;
	n1 = n2 = 0;
	for(i = 0; i < equLen; i ++)	{
		r1 = equSet[0][i];
		r2 = equSet[1][i];
		for(k = 0; k < 3; k ++)	{
			tmp1[n1 ++] = cod1[3 * r1 + k]; 
			tmp2[n2 ++] = cod2[3 * r2 + k]; 
		}
	}

	//superimpose the equivalent residues
	double	r[9], t[3];
	rmsd = GEOMETRY::kearsay(equLen, tmp1, tmp2, r, t);

	//transformed structure 2 accroding to the superimposing of equivalent residues
	for(i = 0; i < pro2Len; i ++)	{
		GEOMETRY::tran_ord(r, t, cod2 + 3 * i);
	}

	DelArray <double> (tmp1);
	DelArray <double> (tmp2);
}

//--------------------------------------------------------------------------------------------------------
//the equivalent residues: residues where Dij <= Dc and i,j is an aligned pair
//use the previous superimposing
//--------------------------------------------------------------------------------------------------------
int SALNOPT::DefineEquPos(int alnLen, int **alnList)
{
	int	i, r1, r2;
	int 	equLenOld = equLen;
	int	**equSetOld = NewMatrix <int> (2, equLenOld);
	for(i = 0; i < equLen; i ++)	{
		equSetOld[0][i] = equSet[0][i];
		equSetOld[1][i] = equSet[1][i];
	}
	double	rmsdOld = rmsd;	
	double	dis;
	equLen = 0;
	//printf("Dc %f, equLenOld %d, rmsdOld %f\n", Dc, equLenOld, rmsdOld);
	for(i = 0; i < alnLen; i ++)	{
		r1 = alnList[0][i];
		r2 = alnList[1][i];
		dis = GEOMETRY::distance(cod1 + 3 * r1, cod2 + 3 * r2);  
		if(dis <= Dc)	{
			equSet[0][equLen] = r1;
			equSet[1][equLen] = r2;
			equLen ++;
		}
	}

	SuperimposeBySet();
	//printf("new equLen %d rmsd %f\n", equLen, rmsd);

	int	ifstop = 0;
	if(fabs(rmsd - rmsdOld) < 1e-10 && equLenOld == equLen)	keepStep ++;
	else	keepStep = 0;

	if(keepStep > maxKeepStep)	{
		ifstop = 1; //converge
	} //allowing up to maxKeepStep instead of 1 is essential for some special cases
	else if(stopRmsd == -1)	{
		ifstop = 0; //condition 1, continue
	}
	else if(rmsd <= stopRmsd * stopRmsdPer || rmsd < rmsdCut) {
		ifstop = 0; //condition 2, continue
	} //rmsdCut is adopted or not? to be tuned
	else	{
		ifstop = 1; //get worse	
	}

	DelMatrix <int> (equSetOld, 2);
	
	if(stopRmsd == -1 && equLen >= stopLenPer * equLen0)	{
		stopRmsd = rmsd; //condition 1
	}

	return ifstop;
}

//--------------------------------------------------------------------------------------------------------
double SALNOPT::OptimizeResult(int *len, int **list)
{
	*len = equLen; 
	for(int i = 0; i < equLen; i ++)	{
		list[0][i] = equSet[0][i];
		list[1][i] = equSet[1][i];
	}
	return rmsd;
}
