//-----------------------------------------------
//The class for extraction of aligned fragment pairs (AFP)
//a AFP (k) is associated with: 
//     len: the length of AFP (resonable value is 8)
//     i(k): the AFP(k)-starting residue position in protein 1
//     j(k): the AFP(k)-starting residue position in protein 2
//     score: the score of AFP(k), it might be RMSD, or others 
//-----------------------------------------------

#ifndef AFPSET_H
#define AFPSET_H

#include "basic.h"

class AFPCHAIN;

class AFP 
{
	private:
		int	i;      //AFP's staring point in protein 1
		int	j;      //AFP's staring point in protein 2
		int	len; 	//AFP length
		double	rmsd;   //AFP RMSD
		double	score;  //AFP Score
		double	t[3];	//the translation
		double	r[9];	//the rotation matrix
		void	Initial()
		{
			i = j = len = 0;
			rmsd = score = 0.0;
			int	k;
			for(k = 0; k < 3; k ++)	t[k] = 0;
			for(k = 0; k < 9; k ++)	r[k] = 0;	
		}
	public:
		//Constructors
		AFP(void)	
		{
			Initial();
		}

		AFP(int i1, int j1, int len1)
		{
			Initial();
			Set(i1, j1, len1);
		}

		AFP(int i1, int j1, int len1, double rmsd1, double *r1, double *t1, double score1) 
		{
			Initial();
			Set(i1, j1, len1, rmsd1, r1, t1, score1);
		}
		~AFP(void)
		{
		}

		//assign AFP's staring position in both proteins and the length
		void	Set(int i1, int j1, int len1)
		{
			i = i1;
			j = j1;
			len = len1;
		}
		//assign rmsd and related translation information
		void 	Set(double rmsd1, double *r1, double *t1)
		{
			rmsd = rmsd1;
			int	k;
			for(k = 0; k < 3; k ++)	t[k] = t1[k];	
			for(k = 0; k < 9; k ++)	r[k] = r1[k];	
		}
		//assign score
		void 	Set(double score1)	{
			score = score1;
		}

		void	Set(int i1, int j1, int len1, double rmsd1, double *r1, double *t1, double score1) 
		{
			Set(i1, j1, len1);
			Set(rmsd1, r1, t1);
			Set(score1);
		}
		
		int	FPosi(void)	{ return i; }
		int	FPosj(void)	{ return j; }
		int	FLen(void)	{ return len; }
		double	FRmsd(void)	{ return rmsd; }
		double	FScore(void)	{ return score; }

		//to check if this AFP is in the diagonal (i,j)
		int	IsDiagonal(int i1, int j1)	{
			if(i - i1 == j - j1)	return 1;
			return 0;
		}

		//to check if this AFP is in the same diagonal of given afp 
		int	IsDiagonal(AFP &afp)	{
			if(i - afp.i == j - afp.j)	return 1;
			return 0;
		}

		//assign a given afp to this
		void	operator = (const AFP &afp) 
		{
			i = afp.i;
			j = afp.j;
			len = afp.len;
			rmsd = afp.rmsd;
			score = afp.score;
			for(int k = 0; k < 3; k ++)	t[k] = afp.t[k];	
			for(int k = 0; k < 9; k ++)	r[k] = afp.r[k];	
		}

		//--------------------------------------------
		//return the gaps between this and afp
		//requiring this AFP > given afp
		//--------------------------------------------
		int	operator % (AFP &afp)
		{
			int	g = (i - afp.i) - (j - afp.j);
			if(g < 0)	g = -g;
			return g;
		}

		//--------------------------------------------
		//return the mis-matched between this and afp
		//requiring this AFP > given afp
		//--------------------------------------------
		int	operator / (AFP &afp)
		{
			int	l1 = i - afp.i - afp.len;
			int	l2 = j - afp.j - afp.len;
			return (l1 > l2?l2:l1);
		}

		//--------------------------------------------
		//merge this AFP and given afp
		//requiring the two AFPs are compatible
		//--------------------------------------------
		void	operator += (AFP &afp)
		{
			if(i > afp.i)	{
				i = afp.i;
				j = afp.j;
				len = afp.len + i - afp.i;
				
			}
			else	{
				len = len + afp.i - i;
			}
		}

		//--------------------------------------------
		//to check if this AFP and given afp are compatible
		//--------------------------------------------
		int	operator > (AFP &afp)
		{
			if(i > afp.i + afp.len - 1 && j > afp.j + afp.len - 1)	return 1; //this > afp
			return 0;
		}

		int	operator < (AFP &afp)
		{
			if(afp.i > i + len - 1 && afp.j > j + len - 1)	return 1; 	  //afp > this
			return 0;
		}

		//--------------------------------------------
		//calculate the shift (translate) between two AFPs 
		//--------------------------------------------
		double	operator >> (AFP &afp)
		{
			double	dis = distance(t, afp.t);
			return dis;
		}

		friend class AFPCHAIN;
};

typedef vector <AFP> AFPARRAY;

#endif
