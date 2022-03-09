//------------------------------------------------------------------
//Purpose: refer AFPchain.h
//By Y.Y, 10/24/02
//Latest update: 5/1/03
//Key functions: ExtractAfp()   extract the AFPs between two proteins
//Key functions: ChainAfp()	chain (assembly) the AFPs to form a complete alignment
//------------------------------------------------------------------

#include <cstring>
#include <cstdlib>

#include "AFPchain.h"
#include "SAlnOpt.h"
#include "SigEva.h"

using namespace ARRAY;
using namespace MATRIX;
using namespace GEOMETRY;

int debug = 0;
int colorn = 7;
int colorlist[] = {255,0,0,  205,0,205,  255,255,0, 0,255,0,  0,0,255,  0,255,255, 139,101,8};
                  //red      magenta     yellow     green     blue      cyan       darkgolden

//------------------------------------------------------------------
//Constructors: assign two proteins and extract their AFPs
//------------------------------------------------------------------
AFPCHAIN::AFPCHAIN(char *file1, char *file2, char *output)
{
	AFPCHAIN(file1, 0, 100000, file2, 0, 100000, 0, output); 
		          //big enough number to cover the whole length of possible longest protein
}

//------------------------------------------------------------------
//Constructors: assign two proteins and extract their AFPs
//------------------------------------------------------------------
AFPCHAIN::AFPCHAIN(char *file1, char *file2, int sps, char *output)
{
	AFPCHAIN(file1, 0, 100000, file2, 0, 100000, sps, output); 
		          //big enough number to cover the whole length of possible longest protein
}

//------------------------------------------------------------------
//Constructors: assign two proteins and extract their AFPs
//------------------------------------------------------------------
AFPCHAIN::AFPCHAIN(char *file1, int s1, int l1, char *file2, int s2, int l2, int sps, char *output)
{
	bgtime = time(NULL);
	Initial();
	SetParameters();
	sparse = sps;
	int	i, p1, p2;
	for(i = strlen(file1) - 1; i >= 0; i --)	{
		if(file1[i] == '/')	break;	
	}	
	p1 = i + 1; //not include the directory
	for(i = strlen(file2) - 1; i >= 0; i --)	{
		if(file2[i] == '/')	break;	
	}	
	p2 = i + 1; 
	strcpy(pro1Name, file1 + p1);
	strcpy(pro2Name, file2 + p2);

	pro1 = new PROT(file1, s1, l1);
	pro1Len = pro1->length;
	if(debug)	printf("protein 1 %s len %d\n", pro1Name, pro1Len);

	pro2 = new PROT(file2, s2, l2);
	pro2Len = pro2->length;
	if(debug)	printf("protein 2 %s len %d\n", pro2Name, pro2Len);

	if(pro1Len <= 0 || pro2Len <= 0)	{
		printf("Wrong pdb read %d %d\n", pro1Len, pro2Len);
		exit(1);
	}

	if(pro1Len < pro2Len)	minLen = pro1Len;
	else			minLen = pro2Len;

	strcpy(output0, output);

	Calloc();

	ExtractAfp();

	//MergeAfp(); //othwerwise, it will cause problem in chainning the AFPs 

	SortAfp();
}

//------------------------------------------------------------------
//Initialize the class
//------------------------------------------------------------------
void AFPCHAIN::Initial(void)
{
	pro1Len = pro2Len = afpNum = 0;
	pro1 = pro2 = NULL;
	iniTwistPdb = optTwistPdb = NULL;

	twi = NULL; //the number of twists making the best score ending at each AFP
	afpChainList = NULL;
	afpChainTwiBin = NULL;
	afpChainTwiList = NULL;
	afpChainLen = 0;
	afpChainTwiNum = 0;
	chainRmsd = 0;
	chainLen = misLen = gapLen = 0;

	blockNum = blockNumIni = blockNumClu = blockNumSpt = 0;
	blockSize = NULL;
	blockRmsd = NULL;

	afpIndex = afpAftIndex = afpBefIndex = NULL;

	optAln = NULL;
	optLen = NULL;
	optRmsd = NULL;
	optLength = 0;

	focusResn = 0;
	focusRes1 = NULL;
	focusRes2 = NULL;
	focusAfpn = 0;
	focusAfpList = NULL;

	totalRmsdIni = totalRmsdOpt = 0.0;
	totalLenIni = totalLenOpt = 0;
	blockResList = NULL;
	blockResSize = NULL;
	blockScore = NULL;
	blockGap = NULL;

	alnseq1 = NULL;
	alnseq2 = NULL;
	alnsymb = NULL;
	alnbeg1 = alnbeg2 = alnLength = 0;

	alignScore = 0.0;
	alignScoreUpdate = 0.0; //the chaining score after delete and merge process
	normAlignScore = 0.0;
	probability = 0.0;

	sparse = 0;

	shortAln = 0; //if the alignment is too short (meanless)
	showtime = 0;

	disTable1 = disTable2 = NULL;
}

//------------------------------------------------------------------
//setup the parameters
//------------------------------------------------------------------
void AFPCHAIN::SetParameters(void)
{
	fragLen = 8;
	fragLenSq = fragLen * fragLen;
	rmsdCut = 3.0; //cutoff for AFP detection
	disCut = 5.0; //for AFPs connection, to be tuned, 4.0
	afpDisCut = fragLenSq * disCut * disCut;
	afpDisCut0 = fragLenSq * disCut;
	disSmooth = 4.0; //for smoothly calculation of twist penalty calculation
	misCut = 2 * fragLen; //structural-dismilar ranges allowed between AFPs
	maxGap = 40; //try-1 30
	maxGapFrag = fragLen + maxGap;
	disFilter = 2.0 * rmsdCut; //for single AFP denifition to be tuned! //two CA-dis is 3.6
	badRmsd = 4.0; //very important paramerter for twists detection
	maxTra = 5;
	gapCreate = -5.0;
	gapExtend = -0.5;
	misScore = gapExtend; //comparable to gapExtend 
	torsionPenalty = 5 * gapCreate; //to be tuned 
	maxPenalty = 1 * gapCreate; //to be tuned
	resScore = 3.0; //on average, the score for each well-matched residue pair
	fragScore = resScore * fragLen; //the score for each well-matched fragment
}

void AFPCHAIN::ShowTimeSet(void)
{
	showtime = 1;
}

//------------------------------------------------------------------
void AFPCHAIN::ChangeDisCut(double value)
{
	disCut = value;
}

//------------------------------------------------------------------
void AFPCHAIN::ChangeBadRmsd(double value)
{
	badRmsd = value;
}
//------------------------------------------------------------------
void AFPCHAIN::ChangeMaxTra(int value)
{
	maxTra = value;
}
//------------------------------------------------------------------
void AFPCHAIN::ChangeTorsionPenalty(double value)
{
	torsionPenalty = value;
}
//------------------------------------------------------------------
void AFPCHAIN::ChangeMaxPenalty(double value)
{
	maxPenalty = value;
}

//-----------------------------------------------------------------
//calloc the major variants: know pro1Len, pro2Len, and maxTra
//-----------------------------------------------------------------
void AFPCHAIN::Calloc(void)
{
	if(afpIndex == NULL)	{
		afpIndex = NewMatrix <int> (pro1Len, pro2Len);//the index of (i,j) pair in AFP list, otherwise -1
		afpAftIndex = NewMatrix <int> (pro1Len, pro2Len); 
		//the index of AFP (i,j*) nearest to (i,j), j*<j. if a AFP exits for (i,j), it equals to afpIndex  
		afpBefIndex = NewMatrix <int> (pro1Len, pro2Len);
		//the index of AFP (i,j*) nearest to (i,j), j*>j. if a AFP exits for (i,j), it equals to afpIndex  
	}
	if(blockResList == NULL)	{
		blockResSize = NewArray <int> (maxTra + 1);
		blockResList = NewArray3 <int> (maxTra + 1, 2, minLen);
	}
	if(blockSize == NULL)	{ //n twists result in n + 1 blocks
		blockSize = NewArray <int> (maxTra + 1);
		block2Afp = NewArray <int> (maxTra + 1);
		blockRmsd = NewArray <double> (maxTra + 1);
		blockScore = NewArray <double> (maxTra + 1);
		blockGap = NewArray <int> (maxTra + 1);
	}
	if(focusAfpList == NULL) focusAfpList = NewArray <int> (minLen); 
	if(focusRes1 == NULL)	focusRes1 = NewArray <int> (minLen); 
	if(focusRes2 == NULL)	focusRes2 = NewArray <int> (minLen); 
}

//------------------------------------------------------------------
//desconstructor
//------------------------------------------------------------------
AFPCHAIN::~AFPCHAIN(void)
{
	if(pro1 != NULL) delete pro1;
	if(pro2 != NULL) delete pro2;
	if(iniTwistPdb != NULL)	delete iniTwistPdb; 
	if(optTwistPdb != NULL)	delete optTwistPdb;
	if(twi != NULL)		DelArray <int> (twi);
	if(afpChainList != NULL)	{
		DelArray <int> (afpChainList);
		DelArray <int> (afpChainTwiBin);
		DelArray <double> (afpChainTwiList);
	}
	if(afpIndex != NULL)	DelMatrix <int> (afpIndex, pro1Len);
	if(afpAftIndex != NULL)	DelMatrix <int> (afpAftIndex, pro1Len);
	if(afpBefIndex != NULL)	DelMatrix <int> (afpBefIndex, pro1Len);

	if(blockRmsd != NULL)	DelArray <double> (blockRmsd);
	if(block2Afp != NULL)	DelArray <int> (block2Afp);
	if(blockSize != NULL)	DelArray <int> (blockSize); //calloc size maxTra + 1
	if(blockScore != NULL)	DelArray <double> (blockScore);
	if(blockGap != NULL)	DelArray <int> (blockGap);
	if(blockResSize != NULL) DelArray <int> (blockResSize); //calloc size maxTra + 1
	if(blockResList != NULL) DelArray3 <int> (blockResList, maxTra + 1, 2); //calloc size maxTra + 1

	if(optAln != NULL)	DelArray3 <int> (optAln, maxTra + 1, 2);
	if(optLen != NULL)	DelArray <int> (optLen); //calloc size maxTra + 1
	if(optRmsd != NULL)	DelArray <double> (optRmsd);

	if(focusAfpList != NULL)	DelArray <int> (focusAfpList);
	if(focusRes1 != NULL)	DelArray <int> (focusRes1);
	if(focusRes2 != NULL)	DelArray <int> (focusRes2);
	if(alnseq1 != NULL)	DelArray <char> (alnseq1);
	if(alnseq2 != NULL)	DelArray <char> (alnseq2);
	if(alnsymb != NULL)	DelArray <char> (alnsymb);

	if(disTable1 != NULL)	DelMatrix <double> (disTable1, pro1Len); 
	if(disTable2 != NULL)	DelMatrix <double> (disTable2, pro2Len); 

	if(debug)	printf("FATCAT completed\n");
}

//------------------------------------------------------------------
//Key function:
//Extract the AFPs of the two proteins
//two segments with RMSD < rmsdCut are considered as a potential AFP
//possible improvements:
//1. consider the amino residues composition (might be good for homology detection)
//2. consider the secondary structure (eg. penalize strand->helix, helix->strand replacements)
//------------------------------------------------------------------
void AFPCHAIN::ExtractAfp(void)
{
	int	p1, p2, filter2, n0, n, n1, n2;
	double	rmsd, r[9], t[3], score, filter1;
	int	add = sparse + 1; //if add > 1, use sparse sampling
	n0 = n = n1 = n2 = 0;
	AFP	afptmp;
	for(p1 = 0; p1 < pro1->length - fragLen; p1 += add )	{
		for(p2 = 0; p2 < pro2->length - fragLen; p2 += add)	{
			n0 ++;
			filter1 = End2End(p1, p1 + fragLen - 1, p2, p2 + fragLen - 1);
				//difference bewteen end-to-end distances
			if(filter1 > disFilter)	{ n1 ++; continue; }
			filter2 = TermFilter(p1, p1 + fragLen - 1, p2, p2 + fragLen - 1);
			if(filter2)	{
				n2 ++;
				continue;
			} //be cautious to use this filter !!
			rmsd = kearsay(fragLen, &(pro1->caCod[3 * p1]), &(pro2->caCod[3 * p2]), r, t);
			//printf("afp %d: p1 %d p2 %d rmsd %f end-to-end dis %f\n", afpSet.size(), p1, p2, rmsd, filter1);
			if(rmsd < rmsdCut)	{
				afptmp.Set(p1, p2, fragLen);
				afptmp.Set(rmsd, r, t);
				score = ScoreAfp(afptmp);
				afptmp.Set(score);
				afpSet.insert(afpSet.end(), afptmp); 
				n ++;
			}
		}
	}	
	afpNum = afpSet.size();
	if(debug)	printf("possible AFP-pairs %d, remain %d after filter 1 remove %d; filter 2 remove %d\n", 
			n0, afpNum, n1, n2); 
}

//------------------------------------------------------------------
//filter 1 for AFP extration: the distance of end-to-end
//------------------------------------------------------------------
double AFPCHAIN::End2End(int p1b, int p1e, int p2b, int p2e)
{
	double 	dis1 = distance(&(pro1->caCod[3 * p1b]), &(pro1->caCod[3 * p1e]));
	double 	dis2 = distance(&(pro2->caCod[3 * p2b]), &(pro2->caCod[3 * p2e]));
	double	min = dis1 - dis2;
	return (fabs(min));
}

//------------------------------------------------------------------
//filter 2 for AFP extration: the context 
//------------------------------------------------------------------
int AFPCHAIN::TermFilter(int p1b, int p1e, int p2b, int p2e)
{
	int	d1 = (p1b < p2b)?p1b:p2b;
	int	d2 = (pro1Len - p1e) < (pro2Len - p2e)?(pro1Len - p1e):(pro2Len - p2e);
	int	d3 = d1 + d2 + fragLen; //maximum alignment length from current AFP
	int	d4 = int(0.3 * minLen);
	if(d3 < d4)	return 1;	
	return 0;
}

//------------------------------------------------------------------
//Assign score to each AFP 
//------------------------------------------------------------------
double AFPCHAIN::ScoreAfp(AFP &afptmp)
{
	//longer AFP with low rmsd is better
	double	s, w;
	//s = (rmsdCut - afptmp.rmsd) * afptmp.len; //the same scroing strategy as that in the post-processing 
	w = afptmp.rmsd/ badRmsd;
	w = w * w;
	s = fragScore * (1.0 - w);
	return s;
}

//--------------------------------------------
//calculate the difference between transformations of two AFPs 
//refer: Gerrit Vriend & Chris Sander, Proteins, 1991, 11(1):52-58 
//--------------------------------------------
double	AFPCHAIN::DiffTrans(AFP &afp1, AFP &afp2)
{
	double	**p1 = NewMatrix <double> (3, 3);
	double	**p2 = NewMatrix <double> (3, 3);
	int	k, m;
	for(k = 0; k < 3; k ++)	{
		for(m = 0; m < 3; m ++)	{
			p1[k][m] = afp1.r[3 * k + m];
			p2[k][m] = afp2.r[3 * k + m];
		}
	}
	invmatrix(3, p2); //extract the inverse matrix
	double	**c = mulmatrix(3, 3, p1, 3, 3, p2);
	double ang = 0.5 * (tracematrix(3, c) - 1.0);
	DelMatrix <double> (p1, 3);
	DelMatrix <double> (p2, 3);
	DelMatrix <double> (c, 3);
	return ang;
}

//------------------------------------------------------------------
//Key function:
//clustering the AFPs and merged those compatible
//If two AFPs(k,m) are overlaped
//ie. in the same diagonal, i(k)-i(m) = j(k)-j(m)
//and the RMSD of the combined AFP of AFP k and AFP m
//then merge AFP k and m.
//repeat this process till all diagonals are considered
//------------------------------------------------------------------
void AFPCHAIN::MergeAfp(void)
{
	int	i, j, n, k;

	int	*invalid = NewArray <int> (afpNum);
	for(k = 0; k < afpNum; k ++)	invalid[k] = 0;
	int	*list = NewArray <int> (afpNum);

	//collect the AFPs in each diagonal
	j = 0;
	for(i = 0; i < pro1Len; i ++)	{
		n = 0;
		for(k = 0; k < afpNum; k ++)	{
			if(afpSet[k].IsDiagonal(i, j))	list[n ++] = k;	
		}
		MergeAfp(n, list, invalid); //merge the AFPs in a diagonal (i, 0)
	}
	i = 0;
	for(j = 0; j < pro2Len; j ++)	{
		n = 0;
		for(k = 0; k < afpNum; k ++)	{
			if(afpSet[k].IsDiagonal(i,j))	list[n ++] = k;
		}
		MergeAfp(n, list, invalid); //merge the AFPs in a diagonal (0, j)
	}
	//consider pro1Len + pro2Len - 1 diagonals in total

	//clean the merged AFPs
	AFPARRAY	tmpset = afpSet;
	afpSet.erase(afpSet.begin(), afpSet.end()); //clear the afpSet
	for(i = 0; i < afpNum; i ++)	{
		if(invalid[i])	continue;
		afpSet.insert(afpSet.end(),tmpset[i]);
	}
	afpNum = afpSet.size();

	if(debug)	printf("after merge, AFP number %d\n", afpNum);

	delete[] invalid;
	delete[] list;
}

//------------------------------------------------------------------
//merge the AFPs a given diagonal (listed in list0)
//------------------------------------------------------------------
void AFPCHAIN::MergeAfp(int len, int *list0, int *invalid)
{
	if(len <= 1)	return;

	int	i, j, k, f1, f2, i1, i2, j1, j2, frag;
	double	rmsd, score, r[9], t[3];
	int	*ipos = NewArray <int> (len);
	int	*sort = NewArray <int> (len);
	int	*list = NewArray <int> (len);

	//sort the AFPs by the starting point
	for(i = 0; i < len; i ++)	{
		k = list0[i];
		ipos[i] = afpSet[k].i;
		sort[i] = i;
	}
	STATIS::qksort <int> (len, ipos, sort);

	for(i = 0; i < len; i ++)	list[i] = list0[sort[i]];

	//merge the AFPs from left to right, in a progressive way 
	for(i = 0; i < len; i ++)	{
		f1 = list[i];
		i1 = afpSet[f1].i;
		j1 = afpSet[f1].j;
		if(invalid[f1])	continue;
		for(j = i + 1; j < len; j ++)	{
			f2 = list[j];
			i2 = afpSet[f2].i;
			j2 = afpSet[f2].j;
			if(i2 > (i1 + afpSet[f1].len))	break;	//not overlapped
			if((i2 + afpSet[f2].len) > (i1 + afpSet[f1].len))	{
				frag = i2 + afpSet[f2].len - i1;
				rmsd = kearsay(frag, &(pro1->caCod[3 * i1]), &(pro2->caCod[3 * j1]), r, t);
				if(rmsd < rmsdCut)	{
//printf("merge %d %d len %d rmsd %f\n", f1, f2, frag, rmsd);
					afpSet[f1].Set(i1, j1, frag);
					afpSet[f1].Set(rmsd, r, t);
					score = ScoreAfp(afpSet[f1]); //assign score to a AFP
					afpSet[f1].Set(score);
					invalid[f2] = 1;
				} //merge the f2 to f1, and invalidate the f2
			}
		}
	}

	delete[] ipos;
	delete[] sort;
	delete[] list;
}

//------------------------------------------------------------------
//Sort the AFPs in increase of their diagonals(i,j)
//------------------------------------------------------------------
void AFPCHAIN::SortAfp(void)
{
	//index the AFP for easy extraction of compabitle AFPs 
	int	i, j, k, a, b;
	for(i = 0; i < pro1Len; i ++)	{
		for(j = 0; j < pro2Len; j ++)	{
			afpIndex[i][j] = afpAftIndex[i][j] = afpBefIndex[i][j] = -1;
		}
	}
	int	b0 = 0;
	for(a = 0; a < afpNum; a ++)	{
		if(a == afpNum - 1 || afpSet[a].i != afpSet[a + 1].i)	{ 
			i = afpSet[a].i;
			for(b = b0; b <= a; b ++)	{
				j = afpSet[b].j;
				afpIndex[i][j] = afpBefIndex[i][j] = afpAftIndex[i][j] = b;
				if(afpSet[b].i != i)	{
					printf("Warning: wrong afp index %d %d\n", i, afpSet[b].i);
					exit(1);
				}
			}
			for(k = 1; k < pro2Len; k ++)	{
				if(afpBefIndex[i][k] == -1) afpBefIndex[i][k] = afpBefIndex[i][k - 1];
			}
			for(k = pro2Len - 2; k >= 0; k --)	{
				if(afpAftIndex[i][k] == -1) afpAftIndex[i][k] = afpAftIndex[i][k + 1];
			}
			b0 = a + 1;
		}
	}

}

//--------------------------------------------------------------------
//Key function: calculate the connectivity of AFP pairs
//no compatibility criteria is excuted
//note: afp1 is previous to afp2 in terms of the position
//--------------------------------------------------------------------
//this module must be optimized
int AFPCHAIN::AfpPairConn(int afp1, int afp2, double *conn, double *dvar)
{
        int     m = afpSet[afp2] / afpSet[afp1];
        int     g = afpSet[afp2] % afpSet[afp1];


	double	gp = misScore * m;	//on average, penalty for a mismatch is misScore, no modification on score 
	if(g > 0)	{
		gp += gapExtend * g;
	}
	if(gp < maxPenalty)	gp = maxPenalty; //penalty cut-off
	//note: use < (smaller) instead of >, because maxPenalty is a negative number

	double	d;
	d = CalAfpDis(afp1, afp2);
	//note: the 'dis' value is numerically equivalent to the 'rms' with exceptions

	int	ch = 0;
	double	tp = 0.0;	
	if(d >= disCut)	{
		tp = torsionPenalty;
		ch = 1;
	} //use the variation of the distances between AFPs
	else  if(d > disCut - disSmooth)	{
		double	wt = sqrt((d - disCut + disSmooth) / disSmooth); 
		//using sqrt: penalty increase with dis more quicker than linear function
		tp = torsionPenalty * wt; 
	}

	*dvar = d;
	*conn = tp + gp;

	return ch;
}

//-----------------------------------------------------------------------
//return the rmsd of the residues from the segments that form the given AFP list
//this value can be a measurement (1) for the connectivity of the AFPs
//-----------------------------------------------------------------------
double AFPCHAIN::CalAfpRmsd(int afpn, int *list)
{
	//assign the length of residues involved in the list of AFPs 

	focusResn = Afp2Res(afpn, list, focusRes1, focusRes2);
	double	rmsd = pro1->CalCaRmsd(pro2, focusResn, focusRes1, focusRes2);

	return rmsd;
}

//-----------------------------------------------------------------------
//return the root mean square of the distance matrix between the residues 
//from the segments that form the given AFP list
//this value can be a measurement (2) for the connectivity of the AFPs
//and its calculation is quicker than the measurement (1), rmsd
//currently only deal with AFP pair
//
//          |-d1--|
//         |--d2---|
//        |---d3----|
//-----------------------------------------------------------------------
//this module is optimized
double AFPCHAIN::CalAfpDis(int afp1, int afp2)
{
	int	i, j, ai, bi, aj, bj;
	double	d;
	double	rms = 0;
	for(i = 0; i < fragLen; i ++)	{
		ai = afpSet[afp1].i + i; 
		bi = afpSet[afp1].j + i; 
		for(j = 0; j < fragLen; j ++)	{
			aj = afpSet[afp2].i + j;
			bj = afpSet[afp2].j + j;	
			d = disTable1[aj][ai] - disTable2[bj][bi];
			rms += d * d; 
			if(rms > afpDisCut)	{ return (disCut); }
		}
	}
	return (sqrt(rms / fragLenSq));
}

//get the afp list and residue list for each block
//----------------------------------------------------------------------
void AFPCHAIN::BlockInfo(void)
{
	int	i, j, k, a, n;
	for(i = 0; i < blockNum; i ++)	{
		n = 0;
		for(j = 0; j < blockSize[i]; j ++)	{	
			//the index in afpChainList, not in the whole afp set
			a = afpChainList[block2Afp[i] + j];
			for(k = 0; k < afpSet[a].len; k ++)	{
				blockResList[i][0][n] = afpSet[a].i + k;
				blockResList[i][1][n] = afpSet[a].j + k;
				n ++;
			}
		}
		blockResSize[i] = n;
	}
}

//---------------------------------------------------------------------
//get the list of residues in the two proteins given a list of AFPs
//---------------------------------------------------------------------
int AFPCHAIN::Afp2Res(int afpn, int *list, int *res1, int *res2)
{
	int	i, j, n, a;
	n = 0;
	for(i = 0; i < afpn; i ++)	{
		a = list[i];
		//printf("afp %d %d\n", i, a);
		for(j = 0; j < afpSet[a].len; j ++)	{
			if(n >= minLen)	{
				printf("Warning: two many residues\n");
				exit(1);
			}
			res1[n] = afpSet[a].i + j;
			res2[n] = afpSet[a].j + j;
			n ++;
		}	
	}
	return n;
}

//-----------------------------------------------------------------------------------
//output the transformed CA coordinations of a given list of residues 
//----------------------------------------------------------------------------------
void AFPCHAIN::DebugPdb(int n, int *reslist, PROT *prot, char chain, double *tran_cod)
{
	int	i, j;
	for(i = 0; i < n; i ++)	{
		j = reslist[i];
       	       	printf("ATOM%7d %-5s%3s%2c%4d%4c%8.3f%8.3f%8.3f\n", 
			j + 1, "CA", pro1->res[j], chain, j + 1, ' ', 
			tran_cod[3*i], tran_cod[3*i+1], tran_cod[3*i+2]);
	}
}

int AFPCHAIN::QuickFilter(double probcut)
{
	int	maxGap_bk = maxGap;
	int	misCut_bk = misCut;
	int	maxTra_bk = maxTra;
	maxGap = 10;
	misCut = fragLen;
	maxTra = 0;
	DoChainAfp();	

	SIGEVA 	sig;
	probability = sig.calSigOldRigid(pro1Len, pro2Len, alignScore);

	maxGap = maxGap_bk;
	misCut = misCut_bk;
	maxTra = maxTra_bk;

	if(probability > probcut)	{
		printf("Align %s %d with %s %d\n", pro1Name, pro1Len, pro2Name, pro2Len); 
		printf("Excluded in QuickFilter with probability cut %.2f\n", probcut);	
		return 0;
	}
	else	return 1;
}

//run rigid AFP chaining process 
//--------------------------------------------------------------------
void AFPCHAIN::RChainAfp(void)
{
	maxTra = 0;
	ChainAfp();
}

//run AFP chaining allowing up to maxTra flexible regions
//--------------------------------------------------------------------
void AFPCHAIN::ChainAfp(void)
{
	if(showtime)	{
		printf("total AFP %d\n", afpNum);
		time_t nowtime = time(NULL);
		double	diff = difftime(nowtime, bgtime);
		printf("Afp extraction: time %.2fs\n", diff); 
	}

	//run AFP chaining
	
	DoChainAfp( );
	
	if(afpChainLen < 1)	{
		shortAln = 1;
		return;
	} //very short alignment

	if(showtime)	{
		time_t nowtime = time(NULL);
		double	diff = difftime(nowtime, bgtime);
		printf("Afp chaining: time %.2fs\n", diff); 
	}

	
	//PostProcess of chaining result
	
	blockNumIni = blockNum;

	//split blocks (introduce twists) with high RMSD 
	SplitBlock();
	blockNumSpt = blockNum;

	//redo: merge blocks with similar transformations & remove small blocks
	//if(blockNum >= 2)	ClustBlock();
	
	DeleteBlock();
	MergeBlock();
	blockNumClu = blockNum;

	OptimizeAln();

	BlockInfo();

	UpdateScore();

	GetAlign();

	TwistPdb();

	SIGEVA 	sig;
	probability = sig.calSigAll(maxTra, sparse, pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength, blockNum - 1);
	normAlignScore = sig.calNS(pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength, blockNum - 1);
	//if(maxTra == 0)	probability = sig.calSigRigid(pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength);
	//else	probability = sig.calSigFlexi(pro1Len, pro2Len, alignScore, totalRmsdOpt, optLength, blockNum - 1);

	if(showtime)	{
		time_t nowtime = time(NULL);
		double	diff = difftime(nowtime, bgtime);
		printf("Alignment optimization: time %.2fs\n", diff); 
	}
}

//--------------------------------------------------------------------
// Key function: chain (assembly) the AFPs
// a AFP (k) is defined as (i, j, k), with i and j are staring points
// AFP extension (eg. AFP(k-1) -> AFP(k) ) requiring
// AFP(k-1) < AFP(k)(refer AFP.h definition), 
// ie i(k-1) < i(k) and j(k-1) < j(k)
// in the figure, only (2) AFP can extend to AFP(k)
// Key features: a coordination transformation is allowed in the AFP extension
//               gap penalties are also considered
//
// 			      protein1	
//          ---------------------------
//          |        \                |
//          |         \(1)            |
//          |     \    \              |
//          |      \(2) \             |
//       p  |       \                 |
//       r  |   \                     |
//       o  |    \(3)  \(i,j, k)      |
//       t  |     \     \             |
//       e  |            \            |
//       i  |                         |
//       n  |                         |
//       2  ---------------------------
//           schametic of AFP chainning
//--------------------------------------------------------------------
void AFPCHAIN::DoChainAfp(void)
{
	if(afpNum <= 0)	return; //fix a bug reported by Zhanwen

	if(twi == NULL)	twi = NewArray <int> (afpNum);    
	//transformation, calculated at DoChainAfp, be used in List extraction

	//forward: calculate the score matrix
	int	i, j, j0, t, n;
	double	stmp, conn, dvar;
	double	*sco = NewArray <double> (afpNum); //the score ending at an AFP
	int	*pre = NewArray <int> (afpNum);    //the previous AFP
	double	maxsco = 0;
	int	maxafp = 0;
	int	*list = NewArray <int> (afpNum);

	disTable1 = pro1 -> DisTable(maxGap + 2 * fragLen + 1); 
	disTable2 = pro2 -> DisTable(maxGap + 2 * fragLen + 1); 

	for(i = 0; i < afpNum; i ++)	{
		sco[i] = afpSet[i].score; //start from itself
		pre[i] = -1;
		twi[i] = 0;
		if(afpSet[i].i < fragLen || afpSet[i].j < fragLen)	n = 0;
		else	n = CompatibleAfps(i, list); //get a compatible list
		//printf("afp %d, compatible %d\n", i, n);
		for(j0 = 0; j0 < n; j0 ++)	{
			j = list[j0];
			t = AfpPairConn(j, i, &conn, &dvar); //note: j, i
			if(twi[j] + t > maxTra)	continue;
				//two many transformation are disfavored
			stmp = sco[j] + afpSet[i].score + conn;
			if(stmp > sco[i])	{ //considered all previous compatible AFPs
				sco[i] = stmp;
			       	twi[i] = twi[j] + t; 
				pre[i] = j;
			}
		}
		if(maxsco < sco[i])	{
			maxsco = sco[i];
			maxafp = i;
		}
	}

	int	currafp = maxafp;
	if(debug)	printf("maximum score %f, %d\n", maxsco, twi[currafp]);

	//trace-back from maxafp (maxsco)
	alignScore = maxsco;
	alignScoreUpdate = alignScore;

	afpChainTwiNum = 0;
	TraceBack(pre, currafp, twi[currafp]);

	DelArray <int> (pre);
	DelArray <double> (sco);
	DelArray <int> (list);
}

//--------------------------------------------------------------------
//derive the compabitle AFP lists for AFP-chaining
//this is important for speeding up the process
//for a given AFP(i1,j1), there are three regions that could be the starting 
//point for the compabitle AFPs of AFP(i1,j1)
//              a1        a2   a3
//            i1-G    i1-f-c i1-f+1 i1
//              |          |   |   |
//           ----------------------------
//           | B               |   |
//b1  j1-G  -|  ---------------|   |
//           |  |          |   |   |
//           |  |     C    | 3 |   |
//           |  |          |   |   |
//b2 j1-f-c -|  |--------------|   |
//           |  |     2    | 1 |   |
//b3 j1-f+1 -|------------------   |
//           |                   A |
//       j1 -|---------------------\
//           |                      \ (AFP(i1,j1))
//           -----------------------------
//
//f: the length of AFPs (we use segments of same length)
//G: g + f, where g is the maximum allowed gaps 
//c: the maximum allowed cross-over in AFP-connection,
//   here we use c = f, and j1-f-c = j1-2f
//incompatible region A: its AFPs overlap with given AFP(i1,j1)
//incompatible region B: the gaps between its AFP with given AFP is larger than g
//incompatible region C: its AFPs connect with given AFP but cross a given threshold. 
//compatible region 1: [i1-f-c,i1-f+1>,[j1-f-c,j1-f+1> or [a2,a3],[b2,b3]
//compatible region 2: [i1-G,i1-f-c],[j1-f-c,j1-f] or [a1,a2],[b2,b3]
//combine 1 and 2    : [i1-G,i1-f],[j1-f-c,j1-f]   or [a1,a3],[b2,b3]
//compatible region 3: [i1-f-c,i1-f],[j1-G,j1-f-c] or [a2,a3],[b1,b2]
//c->misCut
//f->fragLen
//G->fragLen+maxGap->maxGapFrag
//--------------------------------------------------------------------
//this module is optimized
int AFPCHAIN::CompatibleAfps(int afp, int *list)
{
       	int     i, j, i1, j1, f, G, c, a1, a2, a3, b1, b2, b3, s1, s2;

        f = fragLen;
        G = maxGapFrag;
        c = misCut;

        i1 = afpSet[afp].i;
        j1 = afpSet[afp].j;
       	a3 = i1 - f;
       	a2 = a3 - c;
        a1 = i1 - G;
	a2 = a2>0?a2:0;
	a1 = a1>0?a1:0;

        b3 = j1 - f;
        b2 = b3 - c;
        b1 = j1 - G;
	b2 = (b2 > 0)?b2:0;
	b1 = (b1 > 0)?b1:0;

	int	n = 0;
	//compatible region 1-2, [a1,a3][b2,b3]
	for(i = a1; i <= a3; i ++)	{//i <= a3 instead of i < a3
		s1 = afpAftIndex[i][b2]; //note afpAftIndex, not afpIndex
		if(s1 < 0)	continue;//no AFP for the given i with j > b2 
	       	s2 = afpBefIndex[i][b3]; //afps is sorted by j given a i,it's sparse matrix	
		if(s2 < 0)	continue;//no AFP for the given i with j < b3
		for(j = s1; j <= s2; j ++)	{ //j <= s2 instead of j < s2
			if(twi[j] <= maxTra)	{	
				list[n ++] = j;
			}
		}
	}

	//compatible region 3  [a2,a3][b1,b2]
	for(i = a2; i <= a3; i ++)	{
		s1 = afpAftIndex[i][b1];
		if(s1 < 0)	continue;
	       	s2 = afpBefIndex[i][b2]; //afps is sorted by j given a i	
		if(s2 < 0)	continue;
		//note j < s2, as the cases of j == s2 is alread considered in previous region
		for(j = s1; j < s2; j ++)	{
			if(twi[j] <= maxTra)	{	
				list[n ++] = j;
			}
		}
	}

	return n;
}

//--------------------------------------------------------------------
//derive the optimal chainning of AFPs by trace-back
//--------------------------------------------------------------------
void AFPCHAIN::TraceBack(int *pre, int currafp0, int twist)
{
	//trace-back from currafp (maxsco)
	int	*afpchain = NewArray <int> (minLen);
	int	*afptwibin = NewArray <int> (minLen);
	double	*afptwilist = NewArray <double> (minLen);
	int	currafp = currafp0;
	int	s = 0;
	afptwibin[s] = 0;
	afpchain[s ++] = currafp;
	int	prevafp, t;
	int	alnlen = afpSet[afpchain[s]].len;
	double	conn, dvar;
	while((prevafp = pre[currafp]) != -1)	{
		t = AfpPairConn(prevafp, currafp, &conn, &dvar);
		afptwibin[s - 1] = t;  
		afptwilist[s - 1] = dvar;  
		//note s - 1: the transformation of prevafp-currafp is recorded in currafp
		currafp = prevafp;
		alnlen += afpSet[currafp].len; 
		afpchain[s ++] = currafp;
	}
	afpChainLen = s;
	afptwibin[s - 1] = 0; //first afp without transformation
	if(debug)	printf("including %d AFPs, %d residues\n", afpChainLen, alnlen);

	//record the optimal alignment in afpChainList (afpChainLen)
	if(afpChainList == NULL)	{
		afpChainList = NewArray <int> (s);
		afpChainTwiBin = NewArray <int> (s);
		afpChainTwiList = NewArray <double> (s);
	}
	int	i;
	for(i = 0; i < s; i ++)	{
		afpChainList[i] = afpchain[s - 1 - i];
		afpChainTwiBin[i] = afptwibin[s - 1 - i];
		afpChainTwiList[i] = afptwilist[s - 1 - i];
		afpChainTwiNum += afptwibin[s - 1 - i];
	}
	if(afpChainTwiNum != twist)	{
		printf("Warning: the twists number is't consistent %d %d\n", afpChainTwiNum, twist);
	}

	double	checkscore = afpSet[afpChainList[0]].score;
	for(i = 1; i < afpChainLen; i ++)	{
		t = AfpPairConn(afpChainList[i - 1], afpChainList[i], &conn, &dvar);
		checkscore = checkscore + afpSet[afpChainList[i]].score + conn;
	}
	if(fabs(checkscore - alignScore) > 1e-4)	{
		printf("Warning: fail in alignment score checking %.4f %.4f\n", alignScore, checkscore);
	}

	double	rmsd = CalAfpRmsd(afpChainLen, afpChainList);
	chainRmsd = rmsd;

	int	b1 = 0;
	int	bk = 0;
	int	a, b;
	chainLen = 0;
	block2Afp[0] = 0;
	for(i = 0; i < afpChainLen; i ++)	{
		a = afpChainList[i];
		chainLen += afpSet[a].len;
		if(i > 0)	{
			b = afpChainList[i - 1];
			misLen += afpSet[a] / afpSet[b]; 
			gapLen += afpSet[a] % afpSet[b];
		}
		if(afpChainTwiBin[i])	{
			rmsd = CalAfpRmsd(i - b1, afpChainList + b1);
			blockRmsd[bk] = rmsd;
			blockSize[bk] = i - b1;
			b1 = i;
			block2Afp[++bk] = i; //next-block
		}
	}
	rmsd = CalAfpRmsd(i - b1, afpChainList + b1);
	blockSize[bk] = i - b1;
	blockRmsd[bk] = rmsd;
	blockNum = ++bk;

	DelArray <int> (afpchain);
	DelArray <int> (afptwibin);
	DelArray <double> (afptwilist);
}

//--------------------------------------------------------------------
//remove the artifical small rigid-body superimpose in the middle
//clust the similar superimpositions (caused by the small flexible
//region, which is detected as a seperate rigid superimposing region by adding
//two twists before and after it(artifically!) 
//one possible solution: allowing long enough loops in the chaining process,
//which however increase the calculation complexity
//--------------------------------------------------------------------
void AFPCHAIN::DeleteBlock(void)
{
	//remove those blocks (both in terminals and in the middle) with only a AFP
	//but still keep those small blocks spaning large regions
	if(blockNum <= 1)	return;
	int	blockNumOld = blockNum;
	int	i, j, b1, b2, e1, e2, len;
	e1 = e2 = 0;
	for(i = 0; i < blockNum; i ++) {
		b1 = e1;
		b2 = e2;
		if(i < blockNum - 1)	{
			e1 = afpSet[afpChainList[block2Afp[i + 1]]].i;	
			e2 = afpSet[afpChainList[block2Afp[i + 1]]].j;	
		}
		else	{
			e1 = pro1Len;
			e2 = pro2Len;
		}
		if(blockSize[i] > 1)	continue;
		len = (e1 - b1) < (e2 - b2)?(e1 - b1):(e2 - b2);
		//if(i == blockNum - 1)	blockNum --;
		if(len < 2 * fragLen)	{	
			for(j = i; j < blockNum - 1; j ++)	{
				blockRmsd[j] = blockRmsd[j + 1];
				blockSize[j] = blockSize[j + 1];
				block2Afp[j] = block2Afp[j + 1];
			}
			blockNum --;
			i --;
		} //delete a block
	}
	if(blockNumOld > blockNum)	
		if(debug)	printf("Delete %d small blocks\n", blockNumOld - blockNum);
}

//--------------------------------------------------------------------
//Merge consecutive blocks with similar transformation
//--------------------------------------------------------------------
void AFPCHAIN::MergeBlock(void)
{
	//clustering the neighbor blocks if their transformations are similar
	int	i, j, b1, b2, minb1, minb2;
	double	minrmsd;
	int	merge = 0;
	int	blockNumOld = blockNum;
	double	**rmsdlist = NULL;
	if(blockNum > 1)	{
		rmsdlist = NewMatrix <double> (blockNumOld, blockNumOld);
		for(b1 = 0; b1 < blockNum - 1; b1 ++)	{
			for(b2 = b1 + 1; b2 < blockNum; b2 ++)	{
				rmsdlist[b1][b2] = CombineRmsd(b1, b2);
			}
		}
	}
	minb1 = 0;
	while(blockNum > 1)	{
		minrmsd = 1000;
		for(i = 0; i < blockNum - 1; i ++)	{
			j = i + 1; //only consider neighbor blocks
			if(minrmsd > rmsdlist[i][j])	{
				minrmsd = rmsdlist[i][j];
				minb1 = i;
			}	
		}	
		minb2 = minb1 + 1; //merge those most similar blocks
		//maxrmsd = (blockRmsd[minb1] > blockRmsd[minb2])?blockRmsd[minb1]:blockRmsd[minb2];
		if(minrmsd < badRmsd)	{
			if(debug)	printf("merge block %d (rmsd %.3f) and %d (rmsd %.3f), total rmsd %.3f\n", 
				minb1, blockRmsd[minb1], minb2, blockRmsd[minb2], minrmsd); 
			blockSize[minb1] += blockSize[minb2];
			blockRmsd[minb1] = minrmsd;
			for(i = minb2; i < blockNum - 1; i ++)	{
				block2Afp[i] = block2Afp[i + 1];
				blockSize[i] = blockSize[i + 1];
				blockRmsd[i] = blockRmsd[i + 1];
			} //update block information	
			afpChainTwiNum --;
			blockNum --;
			for(b1 = 0; b1 < blockNum - 1; b1 ++)	{
				for(b2 = b1 + 1; b2 < blockNum; b2 ++) {
					if(b1 == minb1 || b2 == minb1)	{
						rmsdlist[b1][b2] = CombineRmsd(b1, b2);
					}
					else if(b2 < minb1)	continue;	
					else if(b1 < minb1)	{
						rmsdlist[b1][b2] = rmsdlist[b1][b2 + 1];
					}
					else	{
						rmsdlist[b1][b2] = rmsdlist[b1 + 1][b2 + 1];
					}
				}
			} //update the rmsdlist
			merge ++;
		} //merge two blocks	
		else if(minrmsd >= 100)	break;
		else	{
			rmsdlist[minb1][minb2] += 100;
		} //not merge, modify the rmsdlist so that this combination is not considered in next iteration
	}
	if(rmsdlist != NULL)	DelMatrix <double> (rmsdlist, blockNumOld);
	if(merge)	{ 
		if(debug)	printf("Merge %d blocks, remaining %d blocks\n", merge, blockNum);	
	}
}

//to update the chaining score after block delete and merge processed
//the blockScore value is important for significance evaluation
//---------------------------------------------------------------------
void AFPCHAIN::UpdateScore(void)
{
	int	i, j, bknow, bkold, g1, g2;
	double	t, conn, dvar;
	alignScoreUpdate = 0;
	bkold = 0;
	for(i = 0; i < blockNum; i ++)	{
		blockScore[i] = 0;
		blockGap[i] = 0;
		for(j = 0; j < blockSize[i]; j ++)	{
			bknow = afpChainList[block2Afp[i] + j];
			if(j == 0)	{	
				blockScore[i] = afpSet[bknow].score;
			}
			else	{
				t = AfpPairConn(bkold, bknow, &conn, &dvar); //note: j, i
				blockScore[i] += afpSet[bknow].score + conn;
				g1 = afpSet[bknow].i - afpSet[bkold].i - afpSet[bkold].len;
				g2 = afpSet[bknow].j - afpSet[bkold].j - afpSet[bkold].len;
				blockGap[i] += (g1 > g2)?g1:g2;
			}
			bkold = bknow;
		}
		alignScoreUpdate += blockScore[i]; 
	}
	if(blockNum >= 2)	{
		alignScoreUpdate += double(blockNum - 1) * torsionPenalty;
	}
}

//return the rmsd of two blocks
//---------------------------------------------------------------------
double AFPCHAIN::CombineRmsd(int b1, int b2)
{
	int	i;
	int	afpn = 0;
	int	*list = NewArray <int> (blockSize[b1] + blockSize[b2]);
	for(i = block2Afp[b1]; i < block2Afp[b1] + blockSize[b1]; i ++)	{
		list[afpn ++] = afpChainList[i];
	}
	for(i = block2Afp[b2]; i < block2Afp[b2] + blockSize[b2]; i ++)	{
		list[afpn ++] = afpChainList[i];
	}
	double	rmsd = CalAfpRmsd(afpn, list);
	DelArray <int> (list);
	return rmsd;
}

//--------------------------------------------------------------------
//optimize the alignment by dynamic programming
//--------------------------------------------------------------------
void AFPCHAIN::OptimizeAln(void)
{
	int	i, a, k, p1, p2, bk, b1, b2, e1, e2, a1, a2;
	int	iniLen;
	int	**iniSet = NewMatrix <int> (2, minLen);
	int	maxi = 100;
	if(optAln == NULL)	{
		optAln = NewArray3 <int> (maxTra + 1, 2, minLen);
		optLen = NewArray <int> (maxTra + 1);
		optRmsd = NewArray <double> (maxTra + 1);
	}

	//optimize each alignment defined by a block
	b1 = b2 = e1 = e2 = optLength = 0;
	for(bk = 0; bk < blockNum; bk ++)	{
		//initial aligned position
		iniLen = 0;
		if(bk > 0)	{
			b1 = e1;
			b2 = e2; 
		}
		if(bk < blockNum - 1)	{
			a1 = afpChainList[block2Afp[bk] + blockSize[bk] - 1]; //the last AFP in current block
			a2 = afpChainList[block2Afp[bk + 1]];  //the first AFP in next block
			e1 = (afpSet[a1].i + fragLen +  afpSet[a2].i) / 2;
			e2 = (afpSet[a1].j + fragLen +  afpSet[a2].j) / 2;
		} //use the middle point of the current and next AFPs. old (starting point of next AFP)
		else	{
			e1 = pro1Len;
			e2 = pro2Len;
		}
		for(i = block2Afp[bk]; i < block2Afp[bk] + blockSize[bk]; i ++)	{
			a = afpChainList[i];
			p1 = afpSet[a].i;
			p2 = afpSet[a].j;
			for(k = 0; k < afpSet[a].len; k ++)	{
				iniSet[0][iniLen] = p1 + k - b1; //note -b1
				iniSet[1][iniLen] = p2 + k - b2; //note -b2
				iniLen ++;
			}
		}

		//optimize the align by dynamic programming & constraint the optimization region
		if(debug)	printf("optimize block %d (%d afp), region %d-%d(len %d), %d-%d(len %d)\n", 
				bk, blockSize[bk], b1, e1, e1-b1, b2, e2, e2-b2);
		SALNOPT	*opt = new SALNOPT(e1-b1, &pro1->caCod[3 * b1], e2-b2, &pro2->caCod[3 * b2], iniLen, iniSet, maxi);
		optRmsd[bk] = opt->OptimizeResult(&optLen[bk], optAln[bk]);
		if(debug)	printf(" optimized len=%d, rmsd %f\n", optLen[bk], optRmsd[bk]);

		for(i = 0; i < optLen[bk]; i ++)	{
			optAln[bk][0][i] += b1; //restore the position
			optAln[bk][1][i] += b2; //restore the position
		}
		optLength += optLen[bk];

		delete opt;
	}
	DelMatrix <int> (iniSet, 2);
	if(debug)	printf("complete AlignOpt\n");
}

//--------------------------------------------------------------------
//in some special cases, there is no maginificent twists in the
//final chaining result; however, their rmsd (original and after
//optimizing) are very large. Therefore, a post-process is taken
//to split the blocks further at the ralative bad connections (
//with relative high distance variation)  
//to be tested: 
//  split or not according to optimized or initial chaining???
//--------------------------------------------------------------------
void AFPCHAIN::SplitBlock(void)
{
	int	i, a, bk, cut;
	double	maxs, maxt;
	int	blockNum0 = blockNum;

	bk = 0;
	while(blockNum < maxTra + 1)	{
		maxs = 0;
		for(i = 0; i < blockNum; i ++)   {
			if(blockRmsd[i] > maxs && blockSize[i] > 2) { //according to the optimized alignment
				maxs = blockRmsd[i];
				bk = i;
			} //!(Note: optRmsd, not blockRmsd, according to the optimized alignment
		}
		if(maxs < badRmsd)	break;
		maxt = 0;
		cut = 0;
		for(i = 1; i < blockSize[bk]; i ++)	{
			a = i + block2Afp[bk];
			if(afpChainTwiList[a] > maxt)	{
				maxt = afpChainTwiList[a];
				cut = i;
			}
		}
		if(debug)	printf("block %d original size %d rmsd %.3f maxt %.2f cut at %d\n", bk, blockSize[bk], maxs, maxt, cut);
		for(i = blockNum - 1; i > bk; i --)	{
			block2Afp[i + 1] = block2Afp[i];
			blockSize[i + 1] = blockSize[i];
			blockRmsd[i + 1] = blockRmsd[i];
		} //update block information
		block2Afp[bk + 1] = cut + block2Afp[bk];
		blockSize[bk + 1] = blockSize[bk] - cut;
		blockSize[bk] = cut; 
		if(debug)	printf("  split into %d and %d sizes\n", blockSize[bk], blockSize[bk + 1]);
		blockRmsd[bk + 1] = CalAfpRmsd(blockSize[bk + 1], afpChainList + block2Afp[bk + 1]);
		blockRmsd[bk] = CalAfpRmsd(blockSize[bk], afpChainList + block2Afp[bk]);
		//splict a block at the biggest position
		blockNum ++;
	}
	if(blockNum - blockNum0 > 0)	{
		if(debug)	printf("Split %d times:\n", blockNum - blockNum0);
		for(i = 0; i < blockNum; i ++)	{
			if(debug)	printf("  block %d size %d from %d rmsd %.3f\n", i, blockSize[i], block2Afp[i], blockRmsd[i]);
		}
	}
}

//--------------------------------------------------------------------
//display the AFPs in the optimal chainning
//--------------------------------------------------------------------
void AFPCHAIN::ShowAfpChainText(void)
{
	if(shortAln)    return;

	int	i, k, a, len, pi, pj, t, mis, gap, prev;
	double	dis, conn;
	char	info[1000];
	char	outfile[100];
	sprintf(outfile, "%s.chain.txt", output0);
	ofstream out(outfile);
	if(!out) { printf("open outfile %s error in ShowAfpChainText\n", outfile); exit(1); }

	sprintf(info, "Flexible alignment for %s and %s\n", pro1Name, pro2Name);
	out<<info;
	time_t  lt = time(NULL);
	double	diff = difftime(lt, bgtime);
	sprintf(info, "Program: FATCAT (by Y.Y), CPU = %.0fs, on %s", diff, ctime(&lt)); 
	out<<info;

	int	ifprtPar = 1;
	if(ifprtPar)	{
		out<<"\n/***Parameters****/\n";
		out<<"   Fragment Length "<<fragLen<<endl;
		out<<"   RMSD cut for AFP "<<rmsdCut<<endl;
		out<<"   The cutoff for a twist introduced or not "<<disCut<<endl; 
		out<<"   Penalty per structural dissimilar residue "<<misScore<<endl;
		out<<"   Penalty for introducing a twist "<<torsionPenalty<<endl;
		out<<"   Maximum length of gaps "<<maxGap<<endl;
       		out<<"   Maximum length of passing regions "<<misCut<<endl;	
		out<<"   Maximum twists allowed "<<maxTra<<endl;
		out<<"   Gap create "<<gapCreate<<" extension "<<gapExtend<<endl;
	}
	out<<"\n/**** Result summary ****/\n";
	sprintf(info, "ChainLen %d mis %d gap %d RMSD %f, with %d transformation introduced (original %d)\n", 
			chainLen, misLen, gapLen, chainRmsd, blockNum - 1, blockNumIni - 1);
	out<<info;
	cout<<info;

	sprintf(info, "Divided into %d blocks alignment-score %f\n", blockNum, alignScore);
	cout<<info;
	out<<info;
	for(i = 0; i < blockNum; i ++)	{
		sprintf(info, "   Block %d: %d (AFP) %d/%.3f (position/rmsd); after optimization %d/%.3f (position/rmsd)\n", 
				i + 1,  blockSize[i], blockSize[i] * fragLen, blockRmsd[i], optLen[i], optRmsd[i]);
		out<<info;
		cout<<info;
	}

	out<<"\nAFP chainning (positions are consistent with PDB)\n"; 

	for(i = 0; i < afpChainLen; i ++)	{
		a = afpChainList[i];
		len = afpSet[a].len;
		pi = afpSet[a].i;
		pj = afpSet[a].j;
		if(i > 0)	{
			prev = afpChainList[i - 1];
			t = AfpPairConn(prev, a, &conn, &dis);
			mis = afpSet[a] / afpSet[prev];
			gap = afpSet[a] % afpSet[prev];
		}
		else	{
			dis = conn = 0;
			t = mis = gap = 0;
		}
		sprintf(info, "AFP %d(%d) len %d rmsd %.3f score %.2f connect dis %.2f mis %d gap %d pnl %.2f twist? %d\n", 
			i, a, afpSet[a].len, afpSet[a].rmsd, afpSet[a].score, dis, mis, gap, conn, afpChainTwiBin[i]);
		out<<info;
		for(k = 0; k < len; k ++)	{
			sprintf(info, "%5s", pro1->index[pi + k]);
			out<<info;
		}	
		out<<"\n";
		for(k = 0; k < len; k ++)	{
			sprintf(info, "%5s", pro2->index[pj + k]);
			out<<info;
		}	
		out<<"\n";
	}
	if(debug)	printf("complete ShowChainText\n");
}

//------------------------------------------------------------------
//Show the AFPs in the optimal chain in postscript format 
//------------------------------------------------------------------
void AFPCHAIN::ShowAfpChainPs(void)
{
	if(shortAln)	return;

	int	fgnum = afpChainLen;
	int	ptnum = 2 * fgnum;
        double  **point = NewMatrix <double> (ptnum, 2); //refer myArrayTemp.h
        double  **color = NewMatrix <double> (ptnum, 3);
	double	*linewidth = NewArray <double> (ptnum);
	int	i, a, k;
	int	bk = 0;
	int	ck = 0;
	for(i = 0; i < afpChainLen; i ++)	{
		a = afpChainList[i];
		point[2 * i][0] = afpSet[a].i;
		point[2 * i][1] = afpSet[a].j;
		point[2 * i + 1][0] = afpSet[a].i + afpSet[a].len;
		point[2 * i + 1][1] = afpSet[a].j + afpSet[a].len;
		linewidth[2 * i] = 3;
		bk += afpChainTwiBin[i];
		ck += afpChainTwiBin[i]; //note: cause the color scheme could be different from ShowAfp()
		if(ck > colorn)	ck = 0;
		for(k = 0; k < 3; k ++)	{	
			color[2 * i][k] = colorlist[3 * ck + k];
		}
	}

	char	outfile[100];
	sprintf(outfile, "%s.chain.ps", output0);
	PsShow	*psShow = new PsShow(outfile, "alignment graph from FATCAT");	//refer PsShow.C and PsShow.h
	psShow->PageArea(200, 200, 300, 200);
	psShow->setupsidedown();
	psShow->drawAxis(point, ptnum, pro1Name, pro2Name);
	psShow->setTransform(); //set the scale & translate
	psShow->prtColorStepLine(point, ptnum, color, linewidth, 0);

	DelArray <double> (linewidth);
	DelMatrix <double> (point, ptnum);
	DelMatrix <double> (color, ptnum);
	delete psShow;

	if(debug)	printf("complete ShowChainPs\n");
}

//------------------------------------------------------------------
//Show all possible AFPs and the final AFP chain in postscript format 
//------------------------------------------------------------------
void AFPCHAIN::ShowAfp(int ifcolor)
{
	if(shortAln)	return;
	if(alnsymb == NULL)	GetAlign();

	char	outfile[100];
	if(ifcolor)	sprintf(outfile, "%s.afp.color.ps", output0);
	else		sprintf(outfile, "%s.afp.ps", output0);
	PsShow	*psShow = new PsShow(outfile, "%!PS-Adobe-3.0 EPSF-3.0\n%%BoundingBox: 0 0 500 500", "alignment graph from FATCAT");	//refer PsShow.C and PsShow.h
	psShow->PageArea(50, 400, 50, 400);
	psShow->setupsidedown();
	psShow->setMarkFontSize(12);
	psShow->setTextFontSize(12);

	int	i, f, k;

	//extract the AFP coordinations
	int	ptnum = 2 * afpNum;
        double	**point = NewMatrix <double> (ptnum, 2); //refer myArrayTemp.h
        double	*gray = NewArray <double> (ptnum);
	double	*linewidth = NewArray <double> (ptnum);
        double	**color = NewMatrix <double> (ptnum, 3);
	int	min1 = optAln[0][0][0] - 30;
	int	min2 = optAln[0][1][0] - 30;
	int	max1 = optAln[blockNum - 1][0][optLen[blockNum - 1] - 1] + 30;
	int	max2 = optAln[blockNum - 1][1][optLen[blockNum - 1] - 1] + 30;
	int	g = 0;
	for(f = 0; f < afpNum; f ++)	{
		if(afpSet[f].i < min1 || afpSet[f].j < min2)	continue;
		if(afpSet[f].i + afpSet[f].len > max1 || afpSet[f].j + afpSet[f].len > max2)	continue;
		//to focus the AFP-graph
		point[2 * g][0] = afpSet[f].i + 1;
		point[2 * g][1] = afpSet[f].j + 1;
		point[2 * g + 1][0] = afpSet[f].i + afpSet[f].len + 1;
		point[2 * g + 1][1] = afpSet[f].j + afpSet[f].len + 1;
		linewidth[2 * g] = 0.3;
		gray[2 * g] = 0.8;
		g ++;
	} //residue count from 1 instead of 0 in the graph
	int	real_ptnum = 2 * g;

	//draw the background AFPs in black-and-white
	int	ifcodtransform = 0;
	psShow->drawAxis(point, real_ptnum, pro1Name, pro2Name);
	psShow->setTransform(); //set the scale & translate
	psShow->prtGrayStepLine(point, 2 * g, gray, linewidth, ifcodtransform); //draw step-lines

	//extract the AFP-chain coordinations
	int	a, b, ta, tb, block, colorblock;
	a = alnbeg1; 
	b = alnbeg2; 
	block = 0;
	ta = optAln[block][0][optLen[block] - 1];
	tb = optAln[block][1][optLen[block] - 1];
	g = 0;
	for(i = 0; i < alnLength; i ++)	{
		point[2 * g][0] = point[2 * g + 1][0] = a;
		point[2 * g][1] = point[2 * g + 1][1] = b;
		if(alnseq1[i] != '-')	{
			point[2 * g + 1][0] ++;
			a ++;
		}
		if(alnseq2[i] != '-')	{
			point[2 * g + 1][1] ++;
			b ++;
		}
		if(alnsymb[i] != ' ')	linewidth[2 * g] = 1;
		else	linewidth[2 * g] = 0.3;
		gray[2 * g] = 0.0;
		if(block < blockNum - 1)	{
			if((optAln[block + 1][0][0] == a) || (optAln[block + 1][1][0] == b))	{
				block ++;
				ta = optAln[block][0][optLen[block] - 1];
				tb = optAln[block][1][optLen[block] - 1];
				colorblock = block;
			}
			else	if((a>(optAln[block+1][0][0]-ta)/2+ta) || (b>(optAln[block+1][1][0]-tb)/2+tb))	{
					colorblock = block + 1;
			}
			else	colorblock = block;
		}
		else	{ colorblock = block; }
		for(k = 0; k < 3; k ++)	{
			if(colorblock > colorn - 1)	color[2 * g][k] = colorlist[3 * (colorn - 1) + k];
			else	color[2 * g][k] = colorlist[3 * colorblock + k];
		}
		g ++;
	}
	real_ptnum = 2 * g;

	//draw the AFP-chain in color or black-and-white
	if(ifcolor)	{
		psShow->prtColorStepLine(point, real_ptnum, color, linewidth, ifcodtransform);
	}
	else	{
		psShow->prtGrayStepLine(point, real_ptnum, gray, linewidth, ifcodtransform);
	}

	//draw the arrows in the twist position
	double	**arrowpoint = NewMatrix <double> (blockNum - 1, 2);
	for(i = 1; i < blockNum; i ++)	{
		arrowpoint[i - 1][0] = optAln[i][0][0];
		arrowpoint[i - 1][1] = optAln[i][1][0];
	}
	psShow->prtArrow(arrowpoint, blockNum - 1, ifcodtransform, "downleftarrow");
	DelMatrix <double> (arrowpoint, blockNum - 1);

	DelArray <double> (linewidth);
	DelArray <double> (gray);
	DelMatrix <double> (point, ptnum);
	DelMatrix <double> (color, ptnum);

	delete psShow;
}

//--------------------------------------------------------
//calculate the total rmsd of the blocks
//output a merged pdb file for both proteins
//protein 1, in chain A
//protein 2 is twisted accroding to the twists detected, in chain B
//--------------------------------------------------------
void AFPCHAIN::TwistPdb(void)
{
	if(shortAln)	return;

	int	i, bk, b2, e2;
	int	*bound = NewArray <int> (blockNum);
	e2 = 0;
	//superimposing according to the initial AFP-chaining
	if(iniTwistPdb == NULL)	{
       		iniTwistPdb = new PROT;
		iniTwistPdb->AssignPdb(pro2); //refer Prot.C
	}
	b2 = 0;
	int	focusAfpn = 0;
	int	focusResn = 0;
	for(bk = 0; bk < blockNum; bk ++)	{
		TransPdb(blockResSize[bk], blockResList[bk][0], blockResList[bk][1]); 
		//transform pro2 according to comparison of pro1 and pro2 at give residues
		if(bk > 0)	{ b2 = e2; }
		if(bk < blockNum - 1)	{ //bend at the middle of two consecutive AFPs
			e2 = afpSet[afpChainList[block2Afp[bk] + blockSize[bk] - 1]].j;
			e2 = (afpSet[afpChainList[block2Afp[bk + 1]]].j - e2)/ 2 + e2;
		}
		else	{ e2 = pro2Len; }
		ModifyCod(iniTwistPdb, pro2, b2, e2);
		bound[bk] = e2;
		for(i = 0; i < blockSize[bk]; i ++)	{	
			focusAfpList[focusAfpn ++] = afpChainList[block2Afp[bk] + i];
		}
	}

	focusResn = Afp2Res(focusAfpn, focusAfpList, focusRes1, focusRes2);
	totalLenIni = focusResn;
	if(debug)	printf("calrmsdini for %d residues\n", focusResn);
	totalRmsdIni = pro1->CalCaRmsd(iniTwistPdb, focusResn, focusRes1, focusRes2);

	//superimposing according to the optimized alignment 
	if(optTwistPdb == NULL)	{
	       	optTwistPdb = new PROT;
		optTwistPdb->AssignPdb(pro2);
	}
	b2 = 0;
	focusResn = 0;
	for(bk = 0; bk < blockNum; bk ++)	{
		TransPdb(optLen[bk], optAln[bk][0], optAln[bk][1]);
		//transform pro2 according to comparison of pro1 and pro2 at give residues
		if(bk > 0)	{ b2 = e2; }
		if(bk < blockNum - 1)	{ //bend at the middle of two consecutive blocks
			e2 = optAln[bk][1][optLen[bk] - 1]; 
			e2 = (optAln[bk + 1][1][0] - e2)/ 2 + e2;
		}
		else	{ e2 = pro2Len; }
		ModifyCod(optTwistPdb, pro2, b2, e2);
		bound[bk] = e2;
		for(i = 0; i < optLen[bk]; i ++)	{
			focusRes1[focusResn] = optAln[bk][0][i];
			focusRes2[focusResn] = optAln[bk][1][i];
			focusResn ++;
		}
	}
	totalLenOpt = focusResn;
	if(debug)	printf("calrmsdopt for %d residues\n", focusResn);
	totalRmsdOpt = pro1->CalCaRmsd(optTwistPdb, focusResn, focusRes1, focusRes2);
	if(totalLenOpt != optLength)	{
		printf("Warning: final alignment length is different %d %d\n", totalLenOpt, optLength);
	}
	if(debug) printf("final alignment length %d, rmsd %.3f\n", focusResn, totalRmsdOpt);

	DelArray <int> (bound);
}

//--------------------------------------------------------
void AFPCHAIN::WriteTwistPdb(void )
{
	if(shortAln)	return;

	char	pdbfile[100], scriptfile[100], info[1000];

	//output twisted superimposing according to the initial AFP-chaining
	sprintf(pdbfile, "%s.ini.twist.pdb", output0);
	sprintf(info, "REMARK  protein %s chain A with twisted protein %s in chain B", pro1Name, pro2Name); 
	sprintf(info, "%s\nREMARK  initial chaining result %d blocks %d AFPs %d residues %.3f rmsd",
			info, blockNum, focusAfpn, totalLenIni, totalRmsdIni);
	WritePdb(info, pdbfile, pro1, iniTwistPdb);
	sprintf(scriptfile, "%s.ini.twist.script", output0);
	WriteScript(pro1, iniTwistPdb, blockNum, 0, pdbfile, scriptfile);

	//output twisted superimposing according to the optimized alignment 
	sprintf(pdbfile, "%s.opt.twist.pdb", output0);
	sprintf(info, "REMARK  protein %s chain A with twisted protein %s in chain B", pro1Name, pro2Name); 
	sprintf(info, "%s\nREMARK  result after optimizing %d blocks %d residues %.3f rmsd",
			info, blockNum, totalLenOpt, totalRmsdOpt);
	WritePdb(info, pdbfile, pro1, optTwistPdb);
	sprintf(scriptfile, "%s.opt.twist.script", output0);
	WriteScript(pro1, optTwistPdb, blockNum, 1, pdbfile, scriptfile);
}

//---------------------------------------------------------------
void AFPCHAIN::ModifyCod(PROT *p1, PROT *p2, int r1, int r2)
{
	int	i, k;
	for(i = r1; i < r2; i ++) {
		for(k = 0; k < 3; k ++)	{
			p1->caCod[3 * i + k] = p2->caCod[3 * i + k];
		}
	} //modify caCod
	int	a2 = p2->totalAtm;
	if(r2 < p2->length)	{ a2 = p2->res2Atm[r2]; }
	for(i = p2->res2Atm[r1]; i < a2; i ++) {
		for(k = 0; k < 3; k ++)	{
			p1->atmCod[3 * i + k] = p2->atmCod[3 * i + k];
		}
	} //modify atmCod
}
//--------------------------------------------------------
//output the superimposed pdbs in seperate files 
//and write the corresponding rasmol scripts
//according to different superimposing of the blocks
//--------------------------------------------------------
void AFPCHAIN::SuperPdb(void)
{
	int	bk;
	char	pdbfile[100], scriptfile[100], info[1000];
	for(bk = 0; bk < blockNum; bk ++)	{
		//superimposing according to the initial AFP-chaining
		sprintf(pdbfile, "%s.ini.%d.pdb", output0, bk + 1);
		sprintf(scriptfile, "%s.ini.%d.script", output0, bk + 1);
		sprintf(info, "REMARK  according to block %d, including %d AFPs (%d residues), with rmsd %.3f", 
				bk, blockSize[bk], blockSize[bk] * fragLen, blockRmsd[bk]);
		TransPdb(blockResSize[bk], blockResList[bk][0], blockResList[bk][1]);
		WritePdb(info, pdbfile);
		WriteScript(pro1, pro2, bk, 0, pdbfile, scriptfile);

		//superimposing according to the optimized alignment
		sprintf(pdbfile, "%s.opt.%d.pdb", output0, bk + 1);
		sprintf(scriptfile, "%s.opt.%d.script", output0, bk + 1);
		sprintf(info, "REMARK  optimized alignment (%d residues with rmsd %.3f) from block %d (%d residues with rmsd %.3f)", 
				optLen[bk], optRmsd[bk], bk, blockSize[bk] * fragLen, blockRmsd[bk]);
		TransPdb(optLen[bk], optAln[bk][0], optAln[bk][1]);
		WritePdb(info, pdbfile);
		WriteScript(pro1, pro2, bk, 1, pdbfile, scriptfile);
	}
}

//--------------------------------------------------------
//write the pdb according to the superimposing of the given position pairs
//--------------------------------------------------------
void AFPCHAIN::TransPdb(int n, int *res1, int *res2)
{
	double	*cod1 = pro1->Cod4Res(n, res1);
	double	*cod2 = pro2->Cod4Res(n, res2);
	double	r[9], t[3], rmsd;
	rmsd = kearsay(n, cod1, cod2, r, t);
	DelArray <double> (cod1);
	DelArray <double> (cod2);

	//transformed the coordinations of whole protein 2 accordingly
	//printf("rotation %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]);
	//printf("translate %.3f %.3f %.3f\n", t[0], t[1], t[2]);
	pro2->TransCod(r, t); //refer Prot.C
}

//--------------------------------------------------------
//write the superimposed pdb
//--------------------------------------------------------
void AFPCHAIN::WritePdb(char *info, char *pdbfile)
{
	WritePdb(info, pdbfile, pro1, pro2);
}

//write two proteins into one file
//--------------------------------------------------------
void AFPCHAIN::WritePdb(char *info, char *pdbfile, PROT *p1, PROT *p2)
{
	ofstream pdb(pdbfile);
	if(!pdb) { printf("open pdbfile %s error in WritePdb\n", pdbfile); exit(1); }
	pdb<<"REMARK  the superimposed protein file"<<endl;
	pdb<<info<<endl;
	pdb<<"REMARK  "<<pro1Name<<" chain A"<<endl;
	pdb<<"REMARK  "<<pro2Name<<" chain B"<<endl;
	p1->PrintPdb('A', pdb);
	pdb<<"TER"<<endl;
	p2->PrintPdb('B', pdb);
	pdb<<"END"<<endl;
	pdb.close();
	
}

//--------------------------------------------------------
//write the rasmol script for superimposing a block 
//--------------------------------------------------------
void AFPCHAIN::WriteScript(PROT *p1, PROT *p2, int bk, int ifopt, char *pdbfile, char *scriptfile)
{
	char	region[5000] = (" ");
	int	i, b1, b2, e1, e2;	
	for(i = 0; i < blockNum; i ++)	{
		if(bk != blockNum && i != bk)	continue; //if bk == blockNum, all blocks
		if(ifopt)	{
			b1 = optAln[i][0][0];
			b2 = optAln[i][1][0];
			e1 = optAln[i][0][optLen[i] - 1];
			e2 = optAln[i][1][optLen[i] - 1];
		}
		else	{
			b1 = blockResList[i][0][0];
			b2 = blockResList[i][1][0];
			e1 = blockResList[i][0][blockResSize[i] - 1];
			e2 = blockResList[i][1][blockResSize[i] - 1];
		}
		sprintf(region, "%s\nselect %s-%sA\ncolor [181,181,181]", region, p1->index[b1], p1->index[e1]);//gray
		sprintf(region, "%s\nselect %s-%sB\ncolor [%d,%d,%d]\n", region, p2->index[b2], p2->index[e2], 
				colorlist[3 * i], colorlist[3 * i + 1], colorlist[3 * i + 2]); //different color
	}
	ofstream script(scriptfile);
	script<<"#!rasmol -script\n# File: "<<scriptfile<<"\n# Creator: FATCAT\n"<<endl;
	script<<"zap\nload pdb \""<<pdbfile<<"\""<<endl;
	script<<"background white\nset ambient 40\nset specular on\n";
	script<<"wireframe off\n\nselect all\ncartoon\n"<<endl;
	script<<"select *A\ncolor [219,219,219]\nselect *B\ncolor [232,232,232]\n"<<endl;
	                  //light gray         //light blue 
	script<<region<<endl;
	script.close();
}

//--------------------------------------------------------
//report the alignment result in special format for list-model
//--------------------------------------------------------
void AFPCHAIN::Report(void)
{
	int	bk, f1, f2, bi, bj, ei, ej, obi, obj, oei, oej;
	char	str1[1000], str2[1000];
	str1[0] = str2[0] = 0;
	int	afpLen = 0;
	for(bk = 0; bk < blockNum; bk ++)	{
		afpLen += blockSize[bk] * fragLen;
	}
	sprintf(str1, "%-20s%4d%4d%6d/%.2f", pro1Name, pro1Len, blockNum - 1, totalLenIni, totalRmsdIni);
	sprintf(str2, "%-20s%4d%4d%6d/%.2f", pro2Name, pro2Len, blockNum - 1, totalLenOpt, totalRmsdOpt);
	for(bk = 0; bk < blockNum; bk ++)	{
		f1 = afpChainList[block2Afp[bk]];
		f2 = afpChainList[block2Afp[bk] + blockSize[bk] - 1];
		bi = afpSet[f1].i;
		bj = afpSet[f1].j;
		ei = afpSet[f2].i + fragLen - 1;
		ej = afpSet[f2].j + fragLen - 1;
		obi = optAln[bk][0][0];
		obj = optAln[bk][1][0];
		oei = optAln[bk][0][optLen[bk] - 1];
		oej = optAln[bk][1][optLen[bk] - 1];
		sprintf(str1, "%s  [%3s-%3s]/[%3s-%3s]/%.1f", 
			str1, pro1->index[bi], pro1->index[ei], pro2->index[bj], pro2->index[ej], blockRmsd[bk]);
		sprintf(str2, "%s  [%3s-%3s]/[%3s-%3s]/%.1f", 
			str2, pro1->index[obi], pro1->index[oei], pro2->index[obj], pro2->index[oej], optRmsd[bk]);
	}
	char	mark[] = ("----------------------------------------------------------------------------\n");
	printf("%s%s\n%s\n", mark, str1, str2);
}

//output short report, for statistical purpose
void AFPCHAIN::ShtReportDes(void)
{
	printf("# pro1  pro2  len1 len2 twist afp rmsd-int equvilent rmsd-equ  chain-Rmsd  score alignLen gaps block-len..\n");
}

//output short report, for statistical purpose
void AFPCHAIN::ShtReport(void)
{
	printf("%-20s %-20s %4d %4d %2d %4d %5.2f %4d %5.2f    %5.2f %8.2f %4d %4d  ", 
		pro1Name, pro2Name, pro1Len, pro2Len, blockNum - 1, 
		totalLenIni, totalRmsdIni, optLength, totalRmsdOpt, chainRmsd, alignScore, alnLength, alnLength-optLength);
	int	i, j, tmp;
	int	*len = NewArray <int> (blockNum);
	for(i = 0; i < blockNum; i ++)	len[i] = optLen[i];
	for(i = 0; i < blockNum; i ++)	{
		for(j = i + 1; j < blockNum; j ++)	{
			if(len[j] > len[i])	{
				tmp = len[j];
				len[j] = len[i];
				len[i] = tmp;
			}
		}
	}
	for(i = 0; i < blockNum; i ++)	{
		printf(" %4d", optLen[i]);
	}
	printf("\n");
	DelArray <int> (len);
}

char* AFPCHAIN::ShtReportStr(void)
{
	char	*str = new char[300];
	
	ALIGN0  *aln = new ALIGN0(alnseq1, alnseq2);
        double  identity = aln->GetIdentity();
	double  similarity = aln->GetSimilarity();
        delete aln;


       	sprintf(str, "%-20s %-20s %4d %4d %2d %4d %5.2f %4d %5.2f    %5.2f %8.2f %4d %4d %.3e %.1f %.2f%% %.2f%% \n",
                   pro1Name, pro2Name, pro1Len, pro2Len, blockNum - 1, totalLenIni, totalRmsdIni, totalLenOpt, totalRmsdOpt, 
		chainRmsd, alignScore, alnLength, gapLen, probability, normAlignScore, identity * 100, similarity * 100);

	return str;
}

//--------------------------------------------------------
//extract the alignment output
// eg
//    STWNTWACTWHLKQP--WSTILILA
//    111111111111     22222222
//    SQNNTYACSWKLKSWNNNSTILILG
// those position pairs labeled by 1 and 2 are equivilent positions, belongs to
// two blocks 1 and 2. The residues between labeled residues are non-equivalent,
// with '-' filling in their resulting gaps 
//--------------------------------------------------------
void AFPCHAIN::GetAlign(void)
{
	int	i, j, k, p1, p2, p1b, p2b, lmax;
	int	len = 0;
	char	tmp[2];
	p1b = p2b = 0;
	if(alnsymb == NULL)	{
		alnseq1 = NewArray <char> (pro1Len + pro2Len + 1);
		alnseq2 = NewArray <char> (pro1Len + pro2Len + 1);
		alnsymb = NewArray <char> (pro1Len + pro2Len + 1);
	}
	for(i = 0; i < blockNum; i ++)	{
		for(j = 0; j < optLen[i]; j ++)	{
			p1 = optAln[i][0][j];
			p2 = optAln[i][1][j];
			if(len > 0)	{
				lmax = (p1 - p1b - 1)>(p2 - p2b - 1)?(p1 - p1b - 1):(p2 - p2b - 1);
				for(k = 0; k < lmax; k ++)	{
					if(k >= (p1 - p1b - 1))	alnseq1[len] = '-';
					else	alnseq1[len] = pro1->saa[p1b + 1 + k];
					if(k >= (p2 - p2b - 1))	alnseq2[len] = '-';
					else	alnseq2[len] = pro2->saa[p2b + 1 + k];
					alnsymb[len ++] = ' ';
				}
			}	
			else	{
				alnbeg1 = p1; //the first position of sequence in alignment 
				alnbeg2 = p2;
			}
			alnseq1[len] = pro1->saa[p1]; 
			alnseq2[len] = pro2->saa[p2]; 
			sprintf(tmp, "%d", i + 1);
			alnsymb[len ++] = tmp[0];
			p1b = p1;
			p2b = p2;
		}
	}	
	alnLength = len;
	alnseq1[len] = alnseq2[len] = alnsymb[len] = '\0';
	gapLen = alnLength - optLength;
}

//--------------------------------------------------------
//output the alignment format result
//mod = 1, print to stdout
//mod = 2, print to a file
//--------------------------------------------------------
void AFPCHAIN::Display(int mod)
{
	char	filename[100];
	ofstream	out;

	if(mod == 2)	{
		sprintf(filename, "%s.aln", output0);
		out.open(filename);
		if(!out) { printf("open file %s error in Display\n", filename); exit(1); }
	}

	char	info[1000];
	sprintf(info, "Align %s %d with %s %d\n", pro1Name, pro1Len, pro2Name, pro2Len); 
	if(mod == 2)	out<<info;
	else	cout<<info;
	if(shortAln)    {
		if(mod == 2)	{
			out<<"Short match\n\n";
			out.close();
		}
		else	{
			cout<<"Short match\n\n";
		}
		return;
	}

	ALIGN0	*aln = new ALIGN0(alnseq1, alnseq2);
	double	identity = aln->GetIdentity();
	double	similarity = aln->GetSimilarity();
	delete aln;

	if(alnsymb == NULL)	GetAlign();
	sprintf(info, "Twists %d ini-len %d ini-rmsd %.2f opt-equ %d opt-rmsd %.2f chain-rmsd %.2f Score %.2f align-len %d gaps %d (%.2f%%)\n", 
		blockNum - 1, totalLenIni, totalRmsdIni, optLength, totalRmsdOpt, chainRmsd, alignScore, 
		alnLength, gapLen, 100.0 * double(gapLen)/double(alnLength));
	if(mod == 2)	out<<info;
	else		cout<<info;
	sprintf(info, "P-value %.2e Afp-num %d Identity %.2f%% Similarity %.2f%%\n", probability, afpNum, identity * 100, similarity * 100);
	if(mod == 2)	out<<info;
	else		cout<<info;

	int	i;
	double	gap;
	for(i = 0; i < blockNum; i ++)	{
		gap = double(blockGap[i]) / double(blockGap[i] + fragLen * blockSize[i]); 
		sprintf(info, "Block %2d afp %2d score %5.2f rmsd %5.2f gap %d (%.2f%%)\n", 
				i, blockSize[i], blockScore[i], blockRmsd[i], blockGap[i], gap);
		if(mod == 2)	out<<info;
		else		cout<<info;
	}

	int	linelen = 70;
	char	*a = NewArray <char> (alnLength + 1);
	char	*b = NewArray <char> (alnLength + 1);
	char	*c = NewArray <char> (alnLength + 1);

	int	t = 0;
	int	ap = alnbeg1;
	int	bp = alnbeg2;
	int	k, len;
	while((alnLength - t) > 0)	{
		if(alnLength - t > linelen)	len = linelen;
		else	len = alnLength - t;
		strcpy(a, alnseq1 + t);
		strcpy(b, alnseq2 + t);
		strcpy(c, alnsymb + t);
		a[len] = b[len] = c[len] = '\0';
		sprintf(info, "\n%14s", " ");
		out<<info;
		for(k = 10; k <= len; k += 10)
			out<<"    .    :";
		if(k <= len + 5) out<<"    .";
	     	sprintf(info, "\nChain 1:%5s %s\n%14s%s\nChain 2:%5s %s\n",
	        	pro1->index[ap], a, " ", c, pro2->index[bp], b);
		if(mod == 2)	out<<info;
		else		cout<<info;
		for(k = 0; k < len; k ++)	{	
			if(a[k] != '-')	ap ++;
			if(b[k] != '-')	bp ++;
		}
		t += len;

	}
	if(mod == 2)	{
		out<<"\nNote: positions are from PDB; the numbers between alignments are block index\n\n";
		out.close();
	}
	else	cout<<"\nNote: positions are from PDB; the numbers between alignments are block index\n\n";

	DelArray <char> (a);
	DelArray <char> (b);
	DelArray <char> (c);
}

//return alignment score
double 	AFPCHAIN::GetAlnSco(void)
{
	return alignScore;
}

//return alignment length (after optimization)
int 	AFPCHAIN::GetOptLen(void)
{
	return optLength;
}
//return the rmsd of alignment (after optimization)
double 	AFPCHAIN::GetOptRms(void)
{
	return totalRmsdOpt;
}
//return alignment length (before optimization)
int 	AFPCHAIN::GetIniLen(void)
{
	return afpChainLen * fragLen;
}
//return the rmsd of alignment (before optimization)
double 	AFPCHAIN::GetIniRms(void)
{
	return totalRmsdIni;
}
//return the blockNum
int 	AFPCHAIN::GetBlkNum(void)
{
	return blockNum; 
}

double  AFPCHAIN::CalNS(void)
{
	return normAlignScore;
}

double 	AFPCHAIN::GetProb(void)
{
	return probability;
}

int  AFPCHAIN::GetGapLen(void)
{
	return gapLen;
}
