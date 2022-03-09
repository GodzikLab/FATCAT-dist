//-------------------------------------------------------------------------
//The main class of flexible structural alignment
//Main variables: pro1(protein 1), pro2(protein 2) and afpSet (AFP list)
//Main functions: Constructors() (assignning two proteins and their AFPs)
//                ChainAfp()     (chainning AFPs for a complete alignment) 
//By Y.Y, 10/24/02
//Latest Update: 10/24/02
//-------------------------------------------------------------------------

#ifndef AFPCHAIN_H
#define AFPCHAIN_H

#include "AFP.h"
#include "Prot.h" 
#include "PsShow.h"
#include "matrix.h"

using namespace MATRIX;

class AFPCHAIN	{
	private:
		//protein information
		char	pro1Name[100];  //name of protein 1
		char	pro2Name[100];  //name of protein 2
		int	pro1Len;        //length of protein 1
		int	pro2Len;        //length of protein 2
		int	minLen;         //the minimum of pro1Len and pro2Len
		PROT	*pro1;          //protein 1 class
		PROT	*pro2;	        //protein 2 class
		PROT    *iniTwistPdb;   //twisted protein 2 according to initial alignment
		PROT	*optTwistPdb;   //twisted protein 2 according to final alignment
		char	output0[100];   //the prefix of output files

		//AFP information
		int	afpNum;	        //the number of AFP (aligned fragment pair)
		AFPARRAY afpSet;        //the list of AFP
		int	**afpIndex;     //the index of AFP for position pair(i,j)
		int	**afpBefIndex;  //the AFPs in the up-left-diagonal of the current AFP
		int	**afpAftIndex;  //the AFPs in the down-right-diagonal of the current AFP

		//AFP chaining parameters
		int	fragLen;	//the threshould for AFP length, ie. fragment length
		int	fragLenSq;
		double	rmsdCut;	//the rmsd threshould for a AFP denifition
		int	maxGap;		//the maximum gap allowed in AFP extension
		int	maxGapFrag;	//i.e. maxGap + fragLen
		int	maxTra;		//the maximum transformations allowed
		double	twiCut;		//the cut of twisting
		double	shfCut;		//the cut of shifting
		double	disFilter;	//the cut for end-to-end distance difference
		double	disCut;		//the cut for transformation or not
		double	afpDisCut;	//the cut for rmsd calculation between two AFPs
		double	afpDisCut0;	//the cut for rmsd calculation between two AFPs
		double	disSmooth;      //the cut for transformation or not smoothing
		int	misCut;		//the cut for allowed crossing region	
		double	badRmsd;        //the rmsd threshold for spliting blocks	
		double	resScore;	//the score associated with a good aligned position pair
		double	fragScore;	//the score associated with a good aligned fragment pair
		double	gapCreate;	//parameters in chainning
		double	gapExtend;	//the penalty for extending a gap
		double	misScore;	//scores associated with structural mis-match positions		
		double	torsionPenalty;	//scores associated with twist of segments
		double	maxPenalty;	//the maximum penalty for structural-dissimilar and gaps
		int	sparse;		//if use sparse sampling of fragments

		//AFP chaining results
		int	*afpChainList;	//the AFPs list in the optimal alignments
		int	afpChainLen;    //the number of AFPs in the optimal alignments
		int	*afpChainTwiBin;//the list of binary transformations in the AFP chain, 1: yes, 0: no
		double	*afpChainTwiList;//the list of transformations in the AFP chain, 1: yes, 0: no
		int	afpChainTwiNum;	//the total number of coordination transformations introduced
		int	afpChainTwiNum0;//the original total number of coordination transformations introduced
		double	chainRmsd;	//the rmsd of the whole AFP chain
		int	chainLen;	//the length of AFP chain
		int	misLen;		//the length of mis-matched regions between AFPs
		int	gapLen;		//the length of gaps between AFPs
		int	*twi;		//the twist introduced at j afp	

		int	blockNum;	//the final block number
		int	blockNumIni;	//block number before block clustering and split
		int	blockNumClu;	//block number after clustering blocks
		int	blockNumSpt;	//block number after spliting blocks
		double	*blockRmsd;	//the RMSD of each block
		int	*block2Afp;	//the index of afp for each block
		int	*blockSize;	//the number of AFPs involved in a block	
		double	*blockScore;	//the score associated with each block
		int	*blockGap;	//the gaps in each block
		int	*blockResSize;	//the number of residues involved in a block
		int	***blockResList;//the list of AFP for each block

		int	***optAln;	//the equivalent positions after optimization of each block
		int	*optLen;	//the number of equivalent positions for each block
		double	*optRmsd;	//the rmsd of superimposing after optimization of each block
		int	optLength;	//the number of equivalent positions for all blocks`

		double	alignScore;	//the original afp chaining score
		double	alignScoreUpdate; //the afp chaining score after clustering and spliting blocks
		double	normAlignScore; //(ns - mu)/beta
		double	probability;    //the probability of a hit with alignScore
		int	totalLenIni;	//the original alignment length
		double	totalRmsdIni;   //the original alignment rmsd
		int	totalLenOpt;    //the final alignment length after alignment optimization
		double	totalRmsdOpt;   //the final alignment rmsd after alignment optimization

		//introduced only for calculate RMSD between two sets of residues from two proteins
		int	focusResn;	//the size of the set
		int	*focusRes1;	//the residues from protein 1
		int	*focusRes2;	//the residues from protein 2
		int	focusAfpn;	//the AFP number
		int	*focusAfpList;	//the AFP list

		//introduced for calculate the distance between two AFPs
		double	**disTable1;
		double	**disTable2;

		//for alignment display
		char	*alnseq1;       //strings record the alignments: protein 1
		char	*alnseq2;       //strings record the alignments: protein 2
		char	*alnsymb;       //symbols to label block-index of each aligned position pairs
		int	alnLength;	//alignment length with gaps & non-equivalent positions
		int	alnbeg1;	//the start position of protein 1 in alignment
		int	alnbeg2;	//the start position of protein 2 in alignment

		int	shortAln;	//the minimum meanful alignment length

		int	showtime;
		time_t	bgtime;		//the starting time for this job 

		void	Initial(void);	//initial the class
		void	SetParameters(void);	//setup the parameters
		void	Calloc(void);	//calloc the variants
		void	ExtractAfp(void); 	//assign AFPs simply based on equal length fragments
		void	MergeAfp(void); 	//merge the AFPs
		void	MergeAfp(int len, int *list, int *invalid); //merge the AFPs in a diagonal
		double	ScoreAfp(AFP &afp); 	//assign the score of a AFP
		int	AfpPairConn(int afp1, int afp2, double *conn, double *dvar);
		double	CalAfpRmsd(int num, int *list);
		double	CalAfpTrans(int num, int *list, double *twi, double *shf);
		double	CalAfpDis(int afp1, int afp2);
		//calculate the transformation of two AFPs
		int	Afp2Res(int afpn, int *list, int *res1, int *res2);
		void	BlockInfo();
		void	SortAfp(void);
		double	End2End(int p1b, int p1e, int p2b, int p2e);
		double	DiffTrans(AFP &afp1, AFP &afp2);
		int	TermFilter(int p1b, int p1e, int p2b, int p2e);
		void	DebugPdb(int n, int *reslist, PROT *prot, char chain, double *tran_cod);
		void	WriteScript(PROT *p1, PROT *p2, int bk, int ifopt, char *pdbfile, char *scriptfile);
		void	TransPdb(int n, int *res1, int *res2);
		void	ModifyCod(PROT *p1, PROT *p2, int r1, int r2);
		void	WritePdb(char *info, char *pdbfile);
		void	WritePdb(char *info, char *pdbfile, PROT *p1, PROT *p2);
		void	DeleteBlock(void);
		void	MergeBlock(void);
		void	UpdateScore(void);
		double	CombineRmsd(int b1, int b2);
		void 	DoChainAfp(void); //main-package for chaining the AFPs for a complete alignment
		int 	CompatibleAfps(int afp, int *list);
		void	TraceBack(int *tra, int currafp0, int twist);
		void	OptimizeAln(void);//main-package for optimize the structural alignment
		void	SplitBlock(void);
		void	GetAlign(void);
	public:
		AFPCHAIN(char *file1, char *file2, char *output);
		AFPCHAIN(char *file1, char *file2, int sps, char *output);
		AFPCHAIN(char *file1, int s1, int l1, char *file2, int s2, int l2, int sps, char *output);
		~AFPCHAIN();
		void	ShowAfp(int ifcolor); //show AFP in ps format
		void 	ChainAfp(void);	 //chain the AFPs
		void 	RChainAfp(void); //chain the AFPs without flexibility(rigid structural alignment)
		void	ShowAfpChainText(void);
		void	ShowAfpChainPs(void);
		void	SuperPdb(void);
		void	TwistPdb(void);
		void	WriteTwistPdb(void);
		void	Report(void);
		void	ShtReportDes(void);
		void	ShtReport(void);
		char*	ShtReportStr(void);
		void	Display(int mod); //output the alignments to stdout (mod=1) or file(mod=2)
		void	ChangeDisCut(double value);
		void	ChangeBadRmsd(double value);
		void	ChangeMaxTra(int value);
		void	ChangeTorsionPenalty(double value);
		void	ChangeMaxPenalty(double value);
		double	GetAlnSco(void);
		int	GetIniLen(void);
		double	GetIniRms(void);
		int	GetOptLen(void);
		double	GetOptRms(void);
		int	GetBlkNum(void);
		double	GetProb(void);
		int	GetGapLen(void);
		void	ShowTimeSet(void);
		int	QuickFilter(double probcut);
		double	CalNS(void);
};

#endif
