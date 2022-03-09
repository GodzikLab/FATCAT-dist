//this program is for output the distance matrix based on
//twists and scores, respectively in one file
#include "AFPchain.h"
#include "hiCluster.h"
#include "tempkit.h"

#define MAXITM 1000

using namespace ARRAY;
using namespace STATIS;

int main(int argc, char *argv[])
{
	int	i, j, phylip;
	double	prob = 1;
	char	inputfile[100], outfile0[100], matrixfile[100], treefile[100], directory[100], str[1000], alnfile[100], oldinput[100];
	inputfile[0] = outfile0[0] = matrixfile[0] = treefile[0] = alnfile[0] = oldinput[0] = 0;
	strcpy(directory, "./");
	phylip = 0;
	int rgd=0;
	for(i = 0; i < argc; i ++)	{
		if(!strcmp(argv[i], "-i") && argc > i + 1)	strcpy(inputfile, argv[++i]); 
		else if(!strcmp(argv[i], "-f") && argc > i + 1)	strcpy(oldinput, argv[++i]);
		else if(!strcmp(argv[i], "-d") && argc > i + 1)	strcpy(directory, argv[++i]);
		else if(!strcmp(argv[i], "-o") && argc > i + 1)	strcpy(outfile0, argv[++i]);
		else if(!strcmp(argv[i], "-c") && argc > i + 1)	sscanf(argv[++i], "%lf", &prob);	
		else if(!strcmp(argv[i], "-r")) rgd = 1;
		else if(!strcmp(argv[i], "-m"))	phylip = 1;
		
	}	
	if(inputfile[0] == 0 || outfile0[0] == 0)	{
		printf("FATCAT2Tree -i list-file <-f oldinput> -o output-prefix <-m> <-c probability-cut(default 1)> <-d directory>\n");
		printf("            [-r]: if given, then run rigid FATCAT, otherwise run flexible FATCAT\n");
		printf("             program will run FATCAT if -f is not provided\n");
		exit(1);
	}

	char	**list = NewMatrix <char> (MAXITM, 200);
	//20 -> 200 YY Apr 8, 2008
	int	item = 0;
	ifstream in(inputfile);
	if(!in)	{
		printf("open file %s error\n", inputfile);
		exit(1);
	}
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 1000);
		if(str[0] == '#' || str[0] == 0)	continue;
		if(item >= MAXITM)	{ printf("please increase list\n"); exit(1); }
		sscanf(str, "%s", list[item ++]);
	}
	in.close();

	double	**dis = NewMatrix <double> (item, item);
	double	**dis2 = NewMatrix <double> (item, item);
	char	pdb1f[100], pdb2f[100], file0[100], code0[10];
	int	s1, s2, l1, l2, optlen, maxtwist, alnlen, gaplen;
	double	ns;
        s1 = s2 = 0;
        l1 = l2 = 100000;
	maxtwist = 5;


	//use fatcat probability as distance
	if(oldinput[0] != 0)	{
		char	tmp1[100], tmp2[100];
		for(i = 0; i < item; i ++)	{
			for(j = 0; j < item; j ++)	dis[i][j] = -1;
		}
		ifstream in2(oldinput);
		if(!in2) { printf("open %s error\n", oldinput); exit(1); }
		while(!in2.eof())	{
			in2.getline(str, 1000);
			sscanf(str, "%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%d%d%lf", tmp1, tmp2, &alnlen, &gaplen, &ns);
			for(i = 0; i < item; i ++)	{
				if(!strncmp(tmp1, list[i], strlen(list[i])))	break;
			}
			if(i == item)	continue;	
			for(j = 0; j < item; j ++)	{
				if(!strncmp(tmp2, list[j], strlen(list[j])))	break;
			}
			if(j == item)	continue;	
			dis[i][j] = dis[j][i] = ns; 
			if(alnlen == 0)	dis2[i][j] = dis2[j][i] = 1.0;
			else	dis2[i][j] = dis2[j][i] = double(gaplen) / double(alnlen);
		}
		in2.close();
		for(i = 0; i < item; i ++)	{
			for(j = 0; j < item; j ++)	{
				if(dis[i][j] == -1)	{
					printf("Distance between item %s %s is not defined\n", list[i], list[j]);
					exit(1);
				}
			}
		}
	}
	else	{	
		sprintf(alnfile, "%s.fatcat", outfile0);
		ofstream alnout(alnfile);
		if(!alnout) { printf("open %s error\n", alnfile); exit(1); }
		printf("now run FATCAT\n");
		for(i = 0; i < item; i ++)	{
			printf("compare item %d (total %d)\n", i, item);
			for(j = i; j < item; j ++)	{
				printf("   to %d\n", j);
				sprintf(pdb1f, "%s/%s.pdb", directory, list[i]);
				sprintf(pdb2f, "%s/%s.pdb", directory, list[j]);
				sprintf(file0, "%s_%s", list[i], list[j]);

				AFPCHAIN        *afpc = new AFPCHAIN(pdb1f, s1, l1, pdb2f, s2, l2, 0, file0);

				if(rgd) afpc->RChainAfp();
				else    afpc->ChainAfp();

				ns = afpc -> GetProb();

				dis[i][j] = dis[j][i] = ns;

				optlen = afpc -> GetOptLen();
				gaplen = afpc -> GetGapLen();
				alnlen = optlen + gaplen;
				if(alnlen == 0)	dis2[i][j] = dis2[j][i] = 1.0;
				else	dis2[i][j] = dis2[j][i] = double(gaplen) / double(alnlen);

				if(alnfile[0] != 0)	{
					char	*tmp = afpc -> ShtReportStr();
					alnout<<tmp;
					delete[] tmp;
				}

				delete afpc;
			}
		}
		alnout.close();
	}

	//output matrix in phylip format
	if(phylip)	{
		sprintf(matrixfile, "%s.phylip", outfile0);
		ofstream out(matrixfile);
		if(!out)	{
			printf("open file %s error\n", argv[2]);
			exit(1);
		}
		out<<" "<<item<<endl;
        	for(i = 0; i < item; i ++)   {
		        code0[0] = 0;
			strncpy(code0, list[i], 50);
			code0[50] = 0;
			out<<code0;
			for(j = 0; j < item; j ++)	{
				if(j == item - 1)    sprintf(str, " %.4f\n", dis[i][j]);
				else                 sprintf(str, " %.4f", dis[i][j]);
				out<<str;
			}
		}
	}

	printf("now clustering..\n");
	hiCluster *hi = new hiCluster(item, list, dis);
	hi -> assignSecondDis(dis2);

	hi -> singleLinkageClust();
	int	mult;
	int	clnum = 0;	
	int	*clsize = NewArray <int> (item);
	int	**cllist = NewMatrix <int> (item, item);
	if(prob < 1.0)	{
		clust *root = hi -> getroot();
		printf("prob-cut %.2f\n", prob);
		hi -> extclust(root, prob, &clnum, clsize, cllist);
		printf("clnum %d\n", clnum);
		if(clnum > 1)	mult = 1;
		else	mult = 0;
	}
	else	mult = 0;
	if(mult == 0)	{
		char	*tree = hi -> exttree();

		sprintf(treefile, "%s.tre", outfile0);
		ofstream treeout(treefile);
		if(!treeout)	{ printf("open file %s error\n", treefile); exit(1); }
		treeout<<tree;
		treeout.close();
		delete[] tree;
	}
	else	{
		char	*subtree, subtreefile[100], sublistfile[100];
		int	d = 0;
		for(i  = 0; i < clnum; i ++)	{
			printf("now cluster %d\n", i);
			if(clsize[i] <= 1)	continue; //siglenton
			d ++;
			subtree = hi -> subtreestr(clsize[i], cllist[i]);
			printf("subtree %s\n", subtree);
			//write subtree
			sprintf(subtreefile, "%s_c%d.tree", outfile0, d);
			ofstream subtreeout(subtreefile);
			if(!subtreeout)	{ printf("open file %s error\n", subtreefile); exit(1); }
			subtreeout<<subtree;
			subtreeout.close();

			//write sublist
			sprintf(sublistfile, "%s_c%d.list", outfile0, d);
			ofstream sublistout(sublistfile);
			if(!sublistout)	{ printf("open file %s error\n", sublistfile); exit(1); }
			sprintf(str, "# probability-cut %.3e\n", prob);
			sublistout<<str;
			for(j = 0; j < clsize[i]; j ++)	{
				sublistout<<list[cllist[i][j]]<<endl;	
			}
			sublistout.close();
			
			delete[] subtree;
		}
		printf("Warning: %d clusters %d non-singletons\n", clnum, d);
	}

	delete hi;

	DelMatrix <char> (list, MAXITM);
	DelMatrix <double> (dis, item);
	DelMatrix <int> (cllist, item);
	delete[] clsize;
}
