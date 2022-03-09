//-------------------------------------------------------------------------
//The main class of flexible structural alignment
//Main variables: pro1(protein 1), pro2(protein 2) and afpSet (AFP list)
//Main functions: Constructors() (assignning two proteins and their AFPs)
//                ChainAfp()     (chainning AFPs for a complete alignment) 
//By Y.Y, 10/24/02
//Latest Update: 10/24/02
//-------------------------------------------------------------------------

#ifndef POSTALIGN_H
#define POSTALIGN_H

#include "Prot.h" 
#include "PsShow.h"

class POSTALIGN	{
	private:
		//protein information
		char	pro1Name[100];  //name of protein 1
		char	pro2Name[100];  //name of protein 2
		int	pro1Len;        //length of protein 1
		int	pro2Len;        //length of protein 2
		int	minLen;         //the minimum of pro1Len and pro2Len
		int	maxTra;
		PROT	*pro1;          //protein 1 class
		PROT	*pro2;	        //protein 2 class
		char	*proseq1;	//the sequence of protein 1 from structure
		char	*proseq2;	//the sequence of protein 2 from structure
		PROT	*optTwistPdb;   //twisted protein 2 according to final alignment
		char	*alnseq1;	//the sequence of protein 1 from alignment
		char	*alnseq2;	//the sequence of protein 2 from alignment
		int	alnbeg1;
		int	alnbeg2;
		char	chain1;
		char	chain2;

		int	blockNum;	//the final block number
		int	***optAln;	//the equivalent positions after optimization of each block
		int	*optLen;	//the number of equivalent positions for each block
		double	*optRmsd;	//the rmsd of superimposing after optimization of each block
		int	optLength;	//the number of equivalent positions for all blocks`

		int	focusResn;	//the size of the set
		int	*focusRes1;	//the residues from protein 1
		int	*focusRes2;	//the residues from protein 2

		int	totalLenOpt;
		double	totalRmsdOpt;

		int	hotSpotNum;
		int	*hotSpotIdx;
		int	rotateIdx;
		double	*rotateList;
		double	*translateList;
		double	*rmsdList;

		void	Initial(void);	//initial the class
		void	Calloc(void);	//calloc the variants
		void	WriteScript(PROT *p1, PROT *p2, int bk, int ifopt, char *pdbfile, char *scriptfile);
		void	TransPdb(int n, int *res1, int *res2);
		void	ModifyCod(PROT *p1, PROT *p2, int r1, int r2);
		void	WritePdb(char *info, char *pdbfile);
		void	WritePdb(char *info, char *pdbfile, PROT *p1, PROT *p2);
	public:
		POSTALIGN(char *file1, char *file2, char *name1, char *name2, char *alnfile);
		POSTALIGN(char *file1, char *file2, char *name1, char *name2, char *alnfile, char *format);
		~POSTALIGN();
		void	GetInput(char *file1, char *file2, char *name1, char *name2, char *alnfile, char *format);
		void	ReadPdb(char *file1, char *file2);
		void	ReadAln(char *alnfile, char *name1, char *name2);
		void	ReadFastaAln(char *alnfile, char *name1, char *name2);
		void	ReadClustalAln(char *alnfile, char *name1, char *name2);
		void	GetAlignPos(char *blockseq);
		void	UniName(char *name);
		char*	Aln2Seq(char *alnseq);
		void	WriteAln(char *alnfile, char *name1, char *name2, char *outfile);
		void	WriteAln(char *alnfile, char *name1, char *name2, ofstream &stream);
		void	TwistPdb(void);
		void	WriteTwistPdb(char *pdbfile, char chain1, char chain2, char *scriptfile);
		void	WriteTwistPdb(char *pdbfile, char chain1, char chain2);
		void	ShowAfp(char *outfile, int ifcolor);
		void	GetHotSpot(int num, char **hotspot);
		void	WriteMatrix(char *outfile);
};

#endif
