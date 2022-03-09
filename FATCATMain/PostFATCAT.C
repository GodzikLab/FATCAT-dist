#include "PostAlign.h"

//this program is relaced by Align2Struct

//this program is for generating the complex pdb with the structures been compared in 
//different chains, and the chaining postscript graph, and the rasmol script file
//GIVEN an alignment of two structures
//writen by Y.Y 6/10/03
//maintained by Y.Y. 11/20/03

int main(int argc, char *argv[])
{
	char	alignfile[100], pdb1file[100], pdb2file[100], name1[100], name2[100], matrixfile[100];
	char	outfile[100], cpdb[100], psfile[100], rasmol[100];
	alignfile[0] = pdb1file[0] = pdb2file[0] = name1[0] = name2[0] = matrixfile[0] = 0;
	outfile[0] = cpdb[0] = psfile[0] = rasmol[0] = 0;
	char	chain1[10] = {"A"}; 
	char	chain2[10] = {"B"};
	int	ref = 1;

	int	i;
	for(i = 0; i < argc; i ++)	{
		if(!strcmp(argv[i], "-a") && argc > i + 1)	strcpy(alignfile, argv[++i]);
		else if(!strcmp(argv[i], "-p1") && argc > i + 1)	strcpy(pdb1file, argv[++i]);
		else if(!strcmp(argv[i], "-p2") && argc > i + 1)	strcpy(pdb2file, argv[++i]);
		else if(!strcmp(argv[i], "-n1") && argc > i + 1)	strcpy(name1, argv[++i]);
		else if(!strcmp(argv[i], "-n2") && argc > i + 1)	strcpy(name2, argv[++i]);
		else if(!strcmp(argv[i], "-o") && argc > i + 1)	strcpy(outfile, argv[++i]);
		else if(!strcmp(argv[i], "-c") && argc > i + 1)	strcpy(cpdb, argv[++i]);
		else if(!strcmp(argv[i], "-s") && argc > i + 1)	strcpy(psfile, argv[++i]);
		else if(!strcmp(argv[i], "-r") && argc > i + 1)	strcpy(rasmol, argv[++i]);
		else if(!strcmp(argv[i], "-m") && argc > i + 1) strcpy(matrixfile, argv[++i]);
		else if(!strcmp(argv[i], "-ref") && argc > i + 1)	sscanf(argv[++i], "%d", &ref);
		else if(!strcmp(argv[i], "-c1") && argc > i + 1) strcpy(chain1, argv[++i]);
		else if(!strcmp(argv[i], "-c2") && argc > i + 1) strcpy(chain2, argv[++i]);
	}
	if(alignfile[0] == 0 || pdb1file[0] == 0 || pdb2file[0] == 0 || name1[0] == 0 || name2[0] == 0)	{
		printf("usage: PostFATCAT <-a alignment> <-p1 pdb1> <-p2 pdb2> <-n1 name1> <-n2 name2> <-c output-pdb>\n");
		printf("  options:\n");
		printf("  [-c1 chain-label-for-protein1-in-the complex] #default 'A'\n");
		printf("  [-c2 chain-label-for-protein2-in-the complex] #default 'B'\n");
		printf("  [-ref 1/2] # 1: transform the second structure (default); otherwise the first one\n");
		printf("  [-s chaining-ps-file]  # write the chaining result in .ps file\n");
		printf("  [-r rasmol-script-file]# write the rasmol-script file\n"); 
		printf("  [-m transformation-matrix] #output the transformation matrix\n");
		printf("  [-o output-align-file] # write the alignment file\n");
		exit(1);
	}
	if(cpdb[0] == 0 && psfile[0] == 0 && rasmol[0] == 0 && matrixfile[0] == 0)	{
		printf("no function is defined\n");
		exit(1);
	}
	POSTALIGN	*palign;
	if(ref == 1) 	palign = new POSTALIGN(pdb1file, pdb2file, name1, name2, alignfile);
	else 		palign = new POSTALIGN(pdb2file, pdb1file, name2, name1, alignfile);
	if(cpdb[0] != 0 && rasmol[0] != 0)	palign -> WriteTwistPdb(cpdb, chain1[0], chain2[0], rasmol);  
	else if(cpdb[0] != 0)	palign -> WriteTwistPdb(cpdb, chain1[0], chain2[0]);
	if(outfile[0] != 0)	{
		if(ref == 1)	palign -> WriteAln(alignfile, name1, name2, outfile);
		else		palign -> WriteAln(alignfile, name2, name1, outfile);
	}
	if(psfile[0] != 0)	palign -> ShowAfp(psfile, 1);
	printf("matrixfile %s\n", matrixfile);
	if(matrixfile[0] != 0)	palign -> WriteMatrix(matrixfile);

	delete palign;

	printf("PostFATCAT done\n");
	return 0;
}
