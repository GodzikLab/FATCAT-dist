#include "FatDom.h"

int main(int argc, char *argv[])
{
	char	reference[100], alignfile[100], pdbfile[100], curvetxt[100], curvefile[100];
	double	pvaluecut = 0.05;
	reference[0] = alignfile[0] = pdbfile[0] = curvefile[0] = curvetxt[0] = 0;
	int	i;
	for(i = 0; i < argc; i ++)	{
		if(!strcmp(argv[i], "-reference") && argc > i + 1)	{
			strcpy(reference, argv[++i]);	
		}
		else if(!strcmp(argv[i], "-align") && argc > i + 1)	{
			strcpy(alignfile, argv[++i]);	
		}
		else if(!strcmp(argv[i], "-refpdb") && argc > i + 1)	{
			strcpy(pdbfile, argv[++i]);	
		}
		else if(!strcmp(argv[i], "-curve") && argc > i + 1)	{
			strcpy(curvefile, argv[++i]);
		}
		else if(!strcmp(argv[i], "-curvetxt") && argc > i + 1)	{
			strcpy(curvetxt, argv[++i]);
		}
		else if(!strcmp(argv[i], "-pvalue") && argc > i + 1)	{
			sscanf(argv[++i], "%lf", &pvaluecut);
		}
	}
	if(reference[0] == 0 || alignfile[0] == 0 || pdbfile[0] == 0)	{
		printf("FATPAT <-reference name> <-refpdb pdbfile> <-align align-file> [-curve curve-ps-file] [-curvetxt curve-txt-file] [-pvalue value]\n");
		exit(1);
	}

	FATDOM	*fatdom = new FATDOM(reference, pdbfile, alignfile, pvaluecut); 

	if(curvefile[0] != 0)	fatdom -> DrawTwistCurve(curvefile);
	if(curvetxt[0] != 0)	fatdom -> PrintTwistFreq(curvetxt);
}
