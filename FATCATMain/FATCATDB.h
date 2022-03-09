#ifndef _FATCATDB_H
#define _FATCATDB_H

#include "basic.h"
#include "PsShow.h"

class FATCATDB
{
	private:
		char	inputFile[100];
		int	itemNum;
		int	itemNum0;
		double	*itemvalue;
		char	**code1;
		char	**code2;
		int	*len1;
		int	*len2;
		int	*twi;
		int	*iniLen;
		int	*optLen;
		double	*iniRms;
		double	*optRms;
		double	*chainRms;
		double	*score;
		double	*probability;
		int	*alnLen;
		double	*gap;
		int	*bkNum;
		double	**bkScore;
		int	*sortorder;
		double	*sortscore;
		int	sort;
		int	cut;

		int	step;
		double	*scalelist;
		double	*distrib;

		void	Calloc(int num0);

	public:
		FATCATDB(int sort_in, int cut_in);
		~FATCATDB(void);
		void	ReadReport(char *alnfile);
		void	ReadAlignList(char *listfile, char *dir);
		void	ReadAlign(char *alnfile);
		void	SortByProb(void);
		void	SortByProbAdOther(void);
		void	SortAlign(void);
		void	WriteAlign(char *outfile);
		void	WriteAlign(char *name1, char *name2, ofstream &out);
		void	WriteAlign(int num, char **name1, char **name2, ofstream &out);
		void	WriteReport(double cut, char *outfile);
		void	WriteHtml(double cut, char *cgi, char *dir, char *outfile);
		void	WriteProb(char *outfile);
		void	ExtractAlign(char *pro1, char *pro2, char *input, char *outfile);
		void	Distribute(int step_in);
		int	GetTarget(char *target);
		void	ExtTarget(double *value);
		void	WriteDistribute(char *target, char *filename);
		void	PsDistribute(char *target, char *filename);
		int	GetItemNum(void);
};

#endif
