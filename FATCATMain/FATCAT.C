#include <cstring> 
#include <cstdlib> 

#include "AFPchain.h"

int main(int argc, char *argv[])
{
	time_t	timebegin = time(NULL);
	int	i;
	char	pdb1[100], pdb2[100], file0[100], dir1[100], dir2[100]; 
	pdb1[0] = pdb2[0] = 0;
	strcpy(file0, "tmp");
	strcpy(dir1, "./");
	strcpy(dir2, "./");
	int	rgd, ct, cg, afp_bk, afp_cl, sup, tws, aln, basic, report, s1, s2, l1, l2, showt, ft, sparse;
	double	probcut = 0;
	rgd = ct = cg = afp_bk = afp_cl = sup = tws = aln = basic = report = showt = ft = sparse = 0;
	s1 = s2 = 0;
	l1 = l2 = 100000;
	for(i = 0; i < argc; i ++)	{
		if(!strcmp(argv[i], "-p1") && argc > i + 1)		strcpy(pdb1, argv[++i]);
		else if(!strcmp(argv[i], "-p2") && argc > i + 1)	strcpy(pdb2, argv[++i]);
		else if(!strcmp(argv[i], "-o")&& argc > i + 1)	strcpy(file0, argv[++i]);
		else if(!strcmp(argv[i], "-b"))	basic = 1;
		else if(!strcmp(argv[i], "-f"))	report = 1;
		else if(!strcmp(argv[i], "-i") && argc > i + 1)	 {
			strcpy(dir1, argv[++i]);	
			strcpy(dir2, dir1);	
		}
		else if(!strcmp(argv[i], "-i1") && argc > i + 1) strcpy(dir1, argv[++i]);	
		else if(!strcmp(argv[i], "-i2") && argc > i + 1) strcpy(dir2, argv[++i]);	
		else if(!strcmp(argv[i], "-s1") && argc > i + 1) sscanf(argv[++i], "%d", &s1); 
		else if(!strcmp(argv[i], "-s2") && argc > i + 1) sscanf(argv[++i], "%d", &s2);
		else if(!strcmp(argv[i], "-l1") && argc > i + 1) sscanf(argv[++i], "%d", &l1); 
		else if(!strcmp(argv[i], "-l2") && argc > i + 1) sscanf(argv[++i], "%d", &l2);
		else if(!strcmp(argv[i], "-r"))	rgd = 1; 
		else if(!strcmp(argv[i], "-c"))	ct = 1; 
		else if(!strcmp(argv[i], "-g"))	cg = 1;	
		else if(!strcmp(argv[i], "-ab"))afp_bk = 1;	
		else if(!strcmp(argv[i], "-ac"))afp_cl = 1;	
		else if(!strcmp(argv[i], "-q"))	aln = 1;	
		else if(!strcmp(argv[i], "-m"))	aln = 2;	
		else if(!strcmp(argv[i], "-s"))	sup = 1;	
		else if(!strcmp(argv[i], "-t"))	tws = 1;	
		else if(!strcmp(argv[i], "-time"))	showt = 1;	
		else if(!strcmp(argv[i], "-filter") && argc > i + 1)	{
			sscanf(argv[++i], "%lf", &probcut);
			ft = 1;	
		}
		else if(!strcmp(argv[i], "-sparse") && argc > i + 1)	{
			sscanf(argv[++i], "%d", &sparse);
		}
	}
	if(pdb1[0] == 0 || pdb2[0] == 0 || sparse < 0 || sparse > 3)	{
		printf("FATCAT <-p1 file> <-p2 file> (the input pdb files)\n");
		printf("  [-o output-initial] (default tmp)\n");
		printf("  [-i string] (data directory for both structures, default ./)\n");
		printf("  [-i1 string] (data directory for 1st structures, default ./)\n");
		printf("  [-i2 string] (data directory for 2st structures, default ./)\n");
		printf("  [-s1 num] (read start position of protein 1, default from the begin)\n");
		printf("  [-s2 num] (read start position of protein 2, default from the begin)\n");
		printf("  [-l1 num] (read length of protein 1, default whole protein)\n");
		printf("  [-l2 num] (read length of protein 2, default whole protein)\n");
		printf("  [-r] (force program run rigid structural alignment, default off)\n");
		printf("  [-filter probcut] (filter the alignment quickly, set a big probcut, eg 0.2, useful in database searching, default off)\n");
		printf("  [-sparse number[0-3]] (sparsely fragment sampling, for speeding up the calculation, default off)\n");
		printf("  [-b] (print a basic report to stdout)\n");
		printf("  [-f] (print a full report to stdout. When -b or -f is on, following options are all automatically off)\n");
		printf("  [-m] (print alignment to a file)\n");
		printf("  [-q] (print alignment to stdout, useful in database-search in queue)\n");
		printf("  [-ab] (print the postscript graph of all AFPs and final AFP chain in black-white to a file)\n");
		printf("  [-ac] (print the postscript graph of all AFPs and final AFP chain in color to a file)\n");
		printf("  [-c] (print AFP chaining result to file.chain.txt)\n");
		printf("  [-t] (print the files of transformed pdbs and corresponding rasmol scripts)\n");
		printf("  [-s] (print the files of superimposed pdbs and corresponding rasmol scripts)\n");
		printf("  [-time] (print the total running time, default off)\n");
		exit(1);
	}
	if(basic || report)	ct = cg = afp_bk = afp_cl = sup = tws = 0;

	char	pdb1f[100], pdb2f[100];
	sprintf(pdb1f, "%s/%s", dir1, pdb1);
	sprintf(pdb2f, "%s/%s", dir2, pdb2);
	AFPCHAIN	*afpc = new AFPCHAIN(pdb1f, s1, l1, pdb2f, s2, l2, sparse, file0);

	if(showt)	afpc->ShowTimeSet();

	if(ft)	{
		int	cont = afpc->QuickFilter(probcut);
		if(cont == 0)	{
			delete afpc;
			return 0;
		}
	}

	if(rgd)	afpc->RChainAfp();
	else	afpc->ChainAfp();

	if(ct)	afpc->ShowAfpChainText();
	//if(cg)	afpc->ShowAfpChainPs();
	if(afp_bk)	afpc->ShowAfp(0);
	if(afp_cl)	afpc->ShowAfp(1);
	if(sup)	afpc->SuperPdb();
	if(tws)	afpc->WriteTwistPdb();
	if(aln)	afpc->Display(aln);

	if(report)	afpc->Report();
	if(basic)	afpc->ShtReport();

	delete afpc;

	if(showt)	{
		time_t	timeend = time(NULL);
		double	diff = difftime(timeend, timebegin);
		printf("total time %.1f\n", diff);
	}

	return 1;
}
