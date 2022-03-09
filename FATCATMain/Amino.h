#ifndef _AMINO_H_
#define _AMINO_H_

#include "basic.h"

namespace AMINO	{

	extern	char 	AA3[AANUM][5];
	extern	char 	AA1[AANUM + 1];
	extern	double	AaMatrix[AANUM][AANUM];
	extern	double	GapExtend;
	extern	double	GapCreate;
	//NOTE: since Amino.h will appear in several files, variables appear in must be declared as "extern"
	//otherwise the compilier will report multiple definition errors.
	//see the actual definitions in Amino.C 

	int 	aa3Index(char *res);
	int 	aa1Index(char res);
	void	GetSstMatrix(char *filename);
	void	GetLinearMatrix(char *sstname, double sw);
	void	GetCombiMatrix(char *filename);
	void	ChangeGapPenalty(double gapcreate, double gapextend);
	double	ppScore(double *prof1, double *prof2);
	double	apScore(int s1, double *prof2);
	double	aaScore(char s1, char s2);	
	void	AssignSSTSymbol(char *type);
	int	sstIndex(char sst);
	double	sstScore(char sst1, char sst2);
	char	aa1Get(char *aa3);
	int	ifThisAProt(char *str);
};

#endif
