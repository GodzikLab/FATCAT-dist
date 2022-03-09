//this class is for alignment by dynamic programing
//using trace-back strategy
//local-model, affine-gap penalty
//accepts two types of input
//1. two sequences, the program calculate the similarity matrix
//2. a similarity matrix (any type, eg. defined by structural similarity)
#ifndef ALIGN0_H
#define ALIGN0_H
#include "basic.h"

class ALIGN0
{
	private:
		int	M; //length of protein 1
		int	N; //length of protein 2
		double	g; //gap-create
		double	h; //gap-extend
		double	m; //g + h
		double	**sij;
		char	**trace; //trace-record
		char	**etrace; //trace-record
		char	**dtrace; //trace-record
		int	B1; //beginning position of protein 1 in alignment
		int	B2; //beginning position of protein 2 in alignment
		int	E1; //end position of protein 1 in alignment
		int	E2; //end position of protein 2 in alignment
		double	alignScore;
		double	identity;
		double	similarity;
		int	*sapp;
		int	*sapp0;
		int	last;

		char	*seq1;
		char	*seq2;
		char	*aln1;
		char	*aln2;
		char	*mark;

		void	Initial(int M0, int N0, double g0, double h0);
		void	DoAlign(void);
		void	Trace(char mod, int i, int j);
		void	Trace(double score, int i, int j);
		void	DEL(int k);
		void	INS(int k);
		void	REP(void);
		double	CheckScore(void);
		void	CheckAlign(void);
	public:
		ALIGN0(double **sij, int M0, int N0, double g0, double h0); //input matrix
		ALIGN0(char *seq1, char *seq2, double g0, double h0); //input proteins
		ALIGN0(char *seq1, char *seq2); //input alignments
		~ALIGN0(void);
		void	Display(void);
		void	Display(char *A, char *B);
		int*	MatchSeqByAlign(void);
		void	FastaAlign(void);	
		void	MapAlign(void);
		void	PlainAlign(void);
		int	GetAlignPos(int **alignList);
		double	GetIdentity(void);
		double	GetSimilarity(void);
};

#endif
