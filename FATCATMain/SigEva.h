#ifndef _SIGEVA_H
#define	_SIGEVA_H
#include "basic.h"

class SIGEVA
{
	private:
		double	mu;
		double	beta;
		double	mu_a;
		double	beta_a;
		double	mu_b;
		double	beta_b;
	public:
		SIGEVA(void);
		void	GetPara(int set, int len);
		double	calSigOldRigid(int len1, int len2, double score);
		double	calSigOldFlexi(int len1, int len2, double score, int r);
		double	calSigRigid(int len1, int len2, double score, double rmsd, int r);
		double	calSigFlexi(int len1, int len2, double score, double rmsd, int optLen, int r);
		double	calSigAll(int twist, int sparse, int len1, int len2, double score, double rmsd, int optLen, int r);
		double	modScore(double score, int r);
		double	normScore(double score, double rmsd, int optLen, int r);
		double	calSig(int len1, int len2, int bknum, double *bkscore);
		int	aveLen(int len1, int len2);
		void	calMu(int len);
		void	calBeta(int len);
		double	calNS(int len1, int len2, double score, double rmsd, int optlen, int r);
};

#endif
