#include "Amino.h"

char AMINO::AA3[AANUM][5] = {"GLY", "ALA", "PRO", "CYS", "THR", "SER",
                           "ASP", "ASN", "GLU", "GLN", "LYS", "HIS", "ARG",
                           "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP", "XXX"};
char AMINO::AA1[AANUM + 1] = {  'G',   'A',   'P',   'C',   'T',   'S',
                              'D',   'N',   'E',   'Q',   'K',   'H',   'R',
                              'V',   'I',   'L',   'M',   'F',   'Y',   'W', '-'};
double	AMINO::GapExtend = -1;
double	AMINO::GapCreate = -11; //default parameter
double	AMINO::AaMatrix[AANUM][AANUM] = {
	{6,0,-2,-3,-2,0,-1,0,-2,-2,-2,-2,-2,-3,-4,-4,-3,-3,-3,-2,-4},
	{0,4,-1,0,0,1,-2,-2,-1,-1,-1,-2,-1,0,-1,-1,-1,-2,-2,-3,-4},
	{-2,-1,7,-3,-1,-1,-1,-2,-1,-1,-1,-2,-2,-2,-3,-3,-2,-4,-3,-4,-4},
	{-3,0,-3,9,-1,-1,-3,-3,-4,-3,-3,-3,-3,-1,-1,-1,-1,-2,-2,-2,-4},
	{-2,0,-1,-1,5,1,-1,0,-1,-1,-1,-2,-1,0,-1,-1,-1,-2,-2,-2,-4},
	{0,1,-1,-1,1,4,0,1,0,0,0,-1,-1,-2,-2,-2,-1,-2,-2,-3,-4},
	{-1,-2,-1,-3,-1,0,6,1,2,0,-1,-1,-2,-3,-3,-4,-3,-3,-3,-4,-4},
	{0,-2,-2,-3,0,1,1,6,0,0,0,1,0,-3,-3,-3,-2,-3,-2,-4,-4},
	{-2,-1,-1,-4,-1,0,2,0,5,2,1,0,0,-2,-3,-3,-2,-3,-2,-3,-4},
	{-2,-1,-1,-3,-1,0,0,0,2,5,1,0,1,-2,-3,-2,0,-3,-1,-2,-4},
	{-2,-1,-1,-3,-1,0,-1,0,1,1,5,-1,2,-2,-3,-2,-1,-3,-2,-3,-4},
	{-2,-2,-2,-3,-2,-1,-1,1,0,0,-1,8,0,-3,-3,-3,-2,-1,2,-2,-4},
	{-2,-1,-2,-3,-1,-1,-2,0,0,1,2,0,5,-3,-3,-2,-1,-3,-2,-3,-4},
	{-3,0,-2,-1,0,-2,-3,-3,-2,-2,-2,-3,-3,4,3,1,1,-1,-1,-3,-4},
	{-4,-1,-3,-1,-1,-2,-3,-3,-3,-3,-3,-3,-3,3,4,2,1,0,-1,-3,-4},
	{-4,-1,-3,-1,-1,-2,-4,-3,-3,-2,-2,-3,-2,1,2,4,2,0,-1,-2,-4},
	{-3,-1,-2,-1,-1,-1,-3,-2,-2,0,-1,-2,-1,1,1,2,5,0,-1,-1,-4},
	{-3,-2,-4,-2,-2,-2,-3,-3,-3,-3,-3,-1,-3,-1,0,0,0,6,3,1,-4},
	{-3,-2,-3,-2,-2,-2,-3,-2,-2,-1,-2,2,-2,-1,-1,-1,-1,3,7,2,-4},
	{-2,-3,-4,-2,-2,-3,-4,-4,-3,-2,-3,-2,-3,-3,-3,-2,-1,1,2,11,-4},
	{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1}};
//this is blosum62 matrix (so GetAaMatrix function can be removed);
//Y.Y Feb 1,2005

/*---------------------Residue Index----------------------------------------*/
int AMINO::aa3Index(char *res)
{
	int	i;
	for(i = 0; i < AANUM; i ++)	{
		if(!strncmp(AA3[i], res, 3))	return i;
	}	
	return (AANUM - 1);
}

/*--------------------Residue Index----------------------------------------*/
int AMINO::aa1Index(char res)
{
	int	i;
	for(i = 0; i < AANUM - 1; i ++)	{
		if(AA1[i] == res)	return i;
	}	
	if(res == '*' || res == '-' || res == '_' || res == ' ')	
		return (AANUM - 1);	/*deletion*/
	else	return -1; /* unknown character as delete */
}

/*--------------------Amino acid mutation matrix----------------------------*/
void AMINO::ChangeGapPenalty(double create, double extend)
{
	if(create > 0)	GapCreate = -create;	
	else	GapCreate = create;
	if(extend > 0)	GapExtend = -extend;
	else	GapExtend = extend;
}

/*--------------------Profile-Profile mutation score-----------------------*/
double	AMINO::ppScore(double *prof1, double *prof2)
{
	int	i, j;
	double	s = 0;
	for(i = 0; i < AANUM; i ++)	{
		for(j = 0; j < AANUM; j ++)	{
			s += prof1[i] * prof2[i] * AaMatrix[i][j];
		}
	}
	return s;
}

/*--------------------Profile-Profile mutation score-----------------------*/
double	AMINO::apScore(int s1, double *prof2)
{
	int	j;
	double	s = 0;
	for(j = 0; j < AANUM; j ++)	{
		s += prof2[j] * AaMatrix[s1][j];
	}
	return s;
}

/*--------------------amino-amino mutation score-----------------------*/
double	AMINO::aaScore(char s1, char s2)
{
	return AaMatrix[aa1Index(s1)][aa1Index(s2)];
}

char	AMINO::aa1Get(char *aa3)
{
	int	idx = aa3Index(aa3);
	return AA1[idx];
}	

int	AMINO::ifThisAProt(char *str)
{
	int	i, j;
	for(i = 0; i < int(strlen(str)); i ++)	{
		j = aa1Index(str[i]);
		if(j < 0)	return 0;
	}
	return 1;
}
