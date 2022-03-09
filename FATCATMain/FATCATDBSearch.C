#include "FATCATDB.h"

/*
 the post-processing of the structure database searching by FATCAT
 1. collect the alignment-result 
 2. draw a score-distribution graph in postscript
 3. output potential candidates with high score
*/

void printusage(void);

int main(int argc, char *argv[])
{
	// read parameters
	char	alnfile[100], reportfile[100], outfile[100], target[100], function[100], pro1[100], pro2[100]; 
	char	cgi[100], dir[100];
	alnfile[0] = reportfile[0] = outfile[0] = function[0] = pro1[0] = pro2[0] = cgi[0] = dir[0] = '\0';
	double	part = 1.0;
	int	graph = 0;
	int	sort = 0;
	int	cut = 0; //big enough number
	strcpy(target, "score");
	strcpy(function, "report");
	int	i;
	for(i = 1; i < argc; i ++)	{
		if(!strcmp(argv[i], "-a") && argc > i + 1)	strcpy(alnfile, argv[++i]);	
		else if(!strcmp(argv[i], "-r") && argc > i + 1)	strcpy(reportfile, argv[++i]);	
		else if(!strcmp(argv[i], "-o") && argc > i + 1)	strcpy(outfile, argv[++i]);	
		else if(!strcmp(argv[i], "-t") && argc > i + 1)	strcpy(target, argv[++i]);
		else if(!strcmp(argv[i], "-f") && argc > i + 1)	strcpy(function, argv[++i]);
		else if(!strcmp(argv[i], "-c1") && argc > i + 1)strcpy(pro1, argv[++i]);
		else if(!strcmp(argv[i], "-c2") && argc > i + 1)strcpy(pro2, argv[++i]);
		else if(!strcmp(argv[i], "-p") && argc > i + 1) sscanf(argv[++i], "%lf", &part);
		else if(!strcmp(argv[i], "-g"))	graph = 1;	
		else if(!strcmp(argv[i], "-s"))	sort = 1;
		else if(!strcmp(argv[i], "-c") && argc > i + 1)	sscanf(argv[++i], "%d", &cut);
		else if(!strcmp(argv[i], "-cgi") && argc > i + 1)	sscanf(argv[++i], "%s", cgi);
		else if(!strcmp(argv[i], "-dir") && argc > i + 1)	sscanf(argv[++i], "%s", dir);
	}
	if(!(alnfile[0] != 0 || reportfile[0] != 0) || outfile[0] == 0 || function[0] == 0)	{
		printusage();
		exit(1);
	}
	if(!strcmp(function, "htmlreport") && (cgi[0] == 0 || dir[0] == 0))	{
		printusage();
		exit(1);		
	}	

	FATCATDB	*db = new FATCATDB(sort, cut);
	int	step = 100;

	//printf("function %s, pro1 %s#, pro2 %s#\n", function, pro1, pro2);

	if(alnfile[0] != 0)	{	
		db -> ReadAlign(alnfile);
	}
	else if(reportfile[0] != 0)	{
		db -> ReadReport(reportfile);
	}
	else	{
		printf("no alnfile or reportfile is available\n");
		exit(1);
	}
	printf("gettarget\n");
	db -> GetTarget(target);


	if(!strcmp(function, "report"))	{
		db -> WriteReport(1.0, outfile);
	}
	else if(!strcmp(function, "sortalign"))	{
		db -> WriteAlign(outfile);
	}
	else if(!strcmp(function, "htmlreport"))	{
		db -> WriteHtml(part, cgi, dir, outfile);
	}
	else if(!strcmp(function, "distribute"))	{
		db -> Distribute(step);
		if(graph)	db -> PsDistribute(target, outfile);
		else		db -> WriteDistribute(target, outfile);
	}
	else if(!strcmp(function, "extalign") && alnfile[0] != 0 && pro1[0] != 0 && pro2[0] != 0)	{
		db -> ExtractAlign(pro1, pro2, alnfile, outfile);
	}
	else if(!strcmp(function, "extprob") && (alnfile[0] != 0 || reportfile[0] != 0))	{
		db -> WriteProb(outfile);
	}
	else	{
		printf("undefined function %s\n", function);
		printusage();
	}

	delete db;

	printf("Success\n");

	return 0;
}

void printusage(void)
{
	printf("FATCATDB [-a alignment-file]/[-r report-file] [-t target] [-f functions] [-g] [-o output]\n");
	printf("  Note: both alignment-file or report-file are allowed as input\n");
	printf("      alignment-file: the detailed alignment result of the database searching\n");
	printf("      report-file: the one line one pair format of the searching result\n");
	printf("  Possible targets: \n");
	printf("     -t score       (fatcat chaining score, default)\n");
	printf("     -t normscore   (score*sqrt(optlen/rmsd*block)\n");
	printf("     -t oprobability (probability only related with alignment-score and query length (and twist))\n");
	printf("     -t rprobability (probability related with alignment-score, rmsd and alignment-length(rigid-fatcat))\n");
	printf("     -t fprobability (probability related with alignment-score, rmsd, alignment-length and twist(flexible-fatcat))\n");
	printf("     -t rprobabilitys2 (probability related with alignment-score, rmsd, alignment-length and twist(flexible-fatcat)),sparsesampling=2\n");
	printf("     -t fprobabilitys2 (probability related with alignment-score, rmsd, alignment-length and twist(flexible-fatcat)),sparsesampling=2\n");
	printf("     -t rprobabilitys1 (probability related with alignment-score, rmsd, alignment-length and twist(flexible-fatcat)),sparsesampling=1\n");
	printf("     -t fprobabilitys1 (probability related with alignment-score, rmsd, alignment-length and twist(flexible-fatcat)),sparsesampling=1\n");
	printf("     -t alignlen    (length of equivilant positions)\n");
	printf("     -t rmsd        (rmsd)\n");
	printf("     -t gap         (gaps)\n");
	printf("     -t twist       (twists)\n");
	printf("  Possible functions: \n");
	printf("     -f report      (report the alignment results in a one line one pair format in text, default)\n");
	printf("     -f htmlreport  (report the alignment results in a one line one pair format in html, for website)\n");
	printf("     -f sortalign   (output the alignment results sorted by probability)\n");
	printf("     -f distribute  (analyze the distribution of database searching result)\n");
	printf("     -f extprob     (extract the probabilities for all given protein pairs)\n");
	printf("     -f extalign    (extract the alignment result for a given protein pairs)\n");
	printf("  -g  to draw a postscript file for the distribution, else output a text file\n");
	printf("  -s  to output report after sorting by score\n");
	printf("  -c  to limit the minimum length of proteins, default no limit\n");
	printf("  -p  to define the probability threshold for output list (default 0.1), work when -f report or -f htmlreport is on\n");
	printf("  -cgi to define the cgi, work only with -f htmlreport is on\n");
	printf("  -dir to define the data dir, work only with -f htmlreport is on\n");
	printf("  -c1 to specify code of protein 1\n");
	printf("  -c2 to give code of protein 2, when use -f extalign, -c1 and -c2 are required\n");
}
