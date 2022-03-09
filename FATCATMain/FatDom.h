#ifndef FATDOM_H
#define FATDOM_H

#include "basic.h"
#include "PsShow.h"
#include "Prot.h"
#define	maxpair	5000
#define	maxopt	100

//a class for FATCAT-based domain dissection

class FATDOM
{
	private:
		PROT	*prot;
		char	*proname;
		int	prolen;
		char	*proseq;
		int	pairn;
		int	tpairn;
		int	*optlen;
		int	**optpos;
		int	*count;
		double	*freq;
		double	*smooth;

		void	Initial(void);
		void	Calloc(void);

	public:
		FATDOM(char *reference, char *pdbfile, char *alnfile, double pvaluecut);
		~FATDOM(void);
		void	GetTwistFreq(void);
		void	ReadAlign(char *reference, char *alignfile, double pcut);
		int	ReadEquPos(int linen, char **lines, int *pos, int reverse);
		void	UniName(char *name);
		void	DrawTwistCurve(char *curvefile);
		void	PrintTwistFreq(char *curvetxt);
};

#endif
