//this class is for structural alignment optimization
//by dynamic programming
#ifndef SALNOPT_H
#define SALNOPT_H

#include "Align0.h"

class SALNOPT	
{
	private:
		int	maxLen;
		int	equLen;
		int	equLen0;
		int	**equSet; //set of equivalent residues
		int	pro1Len;
		int	pro2Len;
		double	*cod1;//coordination of structure 1
		double	*cod2;//coordination of structure 2
		double	rmsd; //the rmsd of superimposing equivalent residues
		double	**sij;//the structure based similarity matrix
		double	Dc;   //the criteria for structural equivalent residues, eg. 3.0 (CE), 6.0(ProSup)
		double	rmsdCut;//the criteria for stoping optimization
		double	increase;
		double	stopLenPer;
		double	stopRmsdPer;
		double	stopRmsd;
		int	maxKeepStep;
		int	keepStep;
		double	gapIni;	
		double	gapExt;
		void	Optimize(int maxi);
		void	SuperimposeBySet(void);
		void	CalMatrix(void);
		int 	DefineEquPos(int alnLen, int **alnList);
		void	SetParameters(void);
	public:
		SALNOPT(int len1, double *cd1, int len2, double *cd2, int iniLen, int **iniSet, int maxi);
		~SALNOPT(void);
		double	OptimizeResult(int *len, int **optSet);
}
;

#endif
