#include <cstdlib> 

#include "SigEva.h"
#define CDF(t) (exp(-exp(-t)))
#define SF(t) (1 - CDF(t))
/* CDF: cumultive density function, SF: survial function for EVD distribution, Gumbel maximum function*/

/* EVD (extreme value distribution)
 * Gumbel maximum distribution (longer right tail) (Type II)
 * pdf (probability density function) f(x) = 1/b exp(-(x-u)/b)) * exp(-exp(-(x-u)/b)))
 * cdf (cumulative distribution function)  F(x) = exp(-exp(-x))
 * sf (survial function)	      S(x) = 1 - F(x) = 1 - exp(-exp(-x))
 * Gumbel minimum distribution (longer left tail) (Type I)
 * pdf (probability density function) f(x) = 1/b exp((x-u)/b)) * exp(-exp((x-u)/b)))
 * cdf (cumulative distribution function)  F(x) = 1 - exp(-exp(x))
 * sf (survial function)	      S(x) = 1 - F(x) = exp(-exp(x))
 * the FATCAT follows Gumbel maximum distribution
 * The survival function is the probability that the variate takes a value greater than x.
*/

SIGEVA::SIGEVA(void)
{
	GetPara(3, 0);
}

void SIGEVA::GetPara(int set, int len)
{
	if(set == 1)	{
		mu_a = 0.2461;
		mu_b = 17.1530;
		beta_a = 0.1284;
		beta_b = 1.3756;
	} //fitting, based on the old benchmark (the fold.list is redundant), old parameters
	else if(set == 2)	{
		mu_a = 1.1137;
		mu_b = 6.5574;
		beta_a = 0.6448;
		beta_b = -12.2793;
	} //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, "NS3"
	else if(set == 3)	{
		mu_a = 0.8440;
		mu_b = 30.2160;
		beta_a = 0.3525;
		beta_b = 18.6652;
	} //score * sqrt(optLen / (rmsd * (twist + 1))), the new best parameters for flexible FATCAT, "FNS8"
	else if(set == 4)	{
		mu_a = 0.4708;
		mu_b = 38.9863;
		beta_a = 0.2511;
		beta_b = 12.9228;
	} //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, sparse-sampling=2, modified scoring
	else if(set == 5)	{
		mu_a = 1.3794; 
		mu_b = -14.4778; 
		beta_a = 0.7465; 
		beta_b = -22.9452; 
	} //score * sqrt(optLen / RMSD), the new best parameters for flexible FATCAT, sparse-sampling=2, modified scoring 3
	else if(set == 6)	{
		mu_a = 0.6036; 
		mu_b = 35.4783; 
		beta_a = 0.3136;
		beta_b = 13.3922; 
	} //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, sparse-sampling=1
	else if(set == 7)	{
		mu_a = 0.7183; 
		mu_b = 27.9647; 
		beta_a = 0.4688; 
		beta_b = -1.3293; 
	} //score * sqrt(optLen / RMSD), the new best parameters for flexible FATCAT, sparse-sampling=1
	else if(set == 8)	{
		mu_a = 0.4813;
		mu_b = 34.2051;
		beta_a = 0.2618;
		beta_b = 12.4581; 
	} //score * sqrt(optLen / RMSD), the new best parameters for rigid FATCAT, sparse-sampling=3, modified scoring
	else if(set == 9)	{
		mu_a = 0.6672; 
		mu_b = 26.5767; 
		beta_a = 0.4373; 
		beta_b = -1.4017;
	} //score * sqrt(optLen / RMSD), the new best parameters for flexible FATCAT, sparse-sampling=3, modified scoring 2
	else	{
		printf("no corresponding parameters\n");
		exit(1);
	}

	mu = beta = .0;
	calMu(len);
	calBeta(len);
}

void SIGEVA::calMu(int len)
{
	mu = mu_a * len + mu_b;
	//printf("mu = %.4f\n", mu);
}

void SIGEVA::calBeta(int len)
{
	beta = beta_a * len + beta_b;
	//printf("beta = %.4f\n", beta);
}

int SIGEVA::aveLen(int len1, int len2)
{
	//int	len = (int(sqrt(len1 * len2)));
	//int	len = len1 < len2?len1:len2; //use the minimum length, bad discriment
	int	len = int(0.5 * (len1 + len2));
	return len;
}

//calculate the probability for a given score (without twists)
double  SIGEVA::calSigOldRigid(int len1, int len2, double score)
{
	int	len = aveLen(len1, len2);
	GetPara(1, len);
	double	t = (score - mu) / beta;
	double	sf = SF(t);
	return sf;
}

//calculate the probability for a given chaining score with r twists
//in this case, the EVD is Gumbel maximum distribution
double  SIGEVA::calSigOldFlexi(int len1, int len2, double score, int r)
{
	int	len = aveLen(len1, len2);
	GetPara(1, len);
	double	mods;
	if(r >= 1)	mods = modScore(score, r); //adjust flexible score to rigid score
	else	mods = score; //no twist

	double	t = (mods - mu) / beta;
	double	sf = SF(t);
	return sf;
}

//calculate the probability for a given score in rigid-FATCAT mod 
//in this case, the EVD is Gumbel maximum distribution
double  SIGEVA::calSigRigid(int len1, int len2, double score, double rmsd, int optLen)
{
	int	len = aveLen(len1, len2);
	GetPara(2, len); //NS3
	double	mods = normScore(score, rmsd, optLen, 0);

	double	t = (mods - mu) / beta;
	double	sf = SF(t);
	return sf;
}

//calculate the probability for a given score with r twists in flexible-FATCAT mod
//in this case, the EVD is Gumbel maximum distribution
double  SIGEVA::calSigFlexi(int len1, int len2, double score, double rmsd, int optLen, int r)
{
	int	len = aveLen(len1, len2);
	GetPara(3, len); //FNS8
	double	mods = normScore(score, rmsd, optLen, r);

	double	t = (mods - mu) / beta;
	double	sf = SF(t);
	return sf;
}

double  SIGEVA::calSigAll(int twist, int sparse, int len1, int len2, double score, double rmsd, int optLen, int r)
{
	int	len = aveLen(len1, len2);
	if(sparse == 2)	{ //use sparse sampling = 2
		if(twist == 0)	GetPara(4, len); //rigid-FATCAT
		else	GetPara(5, len);	//flexible-FATCAT
	}
	else if(sparse == 3)	{ //sparse sampling = 3
		if(twist == 0)	GetPara(8, len); //rigid-FATCAT
		else	GetPara(8, len); //flexible-FATCAT
	}
	else if(sparse == 1)	{ //sparse sampling = 1
		if(twist == 0)	GetPara(6, len); //rigid-FATCAT
		else	GetPara(7, len); //flexible-FATCAT
	}
	else	{ //no sparse sampling or the corresponding sparse parameters are not fitted
		if(twist == 0)	GetPara(2, len); //rigid-FATCAT
		else	GetPara(3, len);	//flexible-FATCAT
	}

	double	mods = normScore(score, rmsd, optLen, r);

	double	t = (mods - mu) / beta;
	double	sf = SF(t);
	return sf;
}

double SIGEVA::calNS(int len1, int len2, double score, double rmsd, int optlen, int r)
{
	int	len = aveLen(len1, len2);
	GetPara(3, len);
	double	ns0 = normScore(score, rmsd, optlen, r);
	double	ns1 = (ns0 - mu) / beta;
	return ns1;
}

//the chaining score is normalized by rmsd, twist and optimal alignment length
double	SIGEVA::normScore(double score, double rmsd, int optLen, int r)
{
	//double	score1 = modScore(score, r);
	double	score1 = score;
	if(r > 0)	score1 /= sqrt(double(r + 1)); 
	 //it is tested that flexible score is more linear relevent to 1/r2 than 1/r
	if(rmsd < 0.5)	score1 *= sqrt(double(optLen) / 0.5);
	else	score1 *= sqrt(double(optLen) / rmsd);
	return score1;
}

//adjust the chaining score of fatcat-flexible with twists to scores without twists
//old evaluation method (adjust the score of flexible alignment to rigid alignment)
double SIGEVA::modScore(double score, int r)
{
	if(r <= 0)	return score;
	int	twist = 5;
	double	shift[] = {0.23, 3.47, 8.74, 11.15, 16.30};  //this set of parameter works with parameter 1
	//the shift values are fitted by random-set: scores with twists and without twists have coefficient better than 5
	if(r > twist)	return (score - shift[4]);	
	else	return (score - shift[r - 1]);
}

//calculate the probability for a list of given scores (with twists)
//ref: PNAS, S Karlin, SF Altschul, 90:5873-5877, 1993
//the tail probability of sum of the r highest normalized scores Tr = S1 + S2..+Sr
//behaves as 
//prob(Tr>x) = e^(-x)x^(r-1)/(r!(r-1)!))
//note: this method is not used in the calculation
double  SIGEVA::calSig(int len1, int len2, int bknum, double *bkscore)
{
	int	i, j, k;
	double	r1, r2, spower, sum, cdf;

	int	len = aveLen(len1, len2);
	calMu(len);
	calBeta(len);

	double	scoresort[10], ns[10];
	double	tmp;
	for(i = 0; i < bknum; i ++)	scoresort[i] = bkscore[i];
	for(i = 0; i < bknum - 1; i ++)	{
		for(j = i + 1; j < bknum; j ++)	{
			if(scoresort[i] < scoresort[j])	{
				tmp = scoresort[i];
				scoresort[i] = scoresort[j];
				scoresort[j] = tmp;	
			}
		}
	} //sort scores (desreasing order)	

	for(i = 0; i < bknum; i ++)	{
		ns[i] = (scoresort[i] - mu) / beta;
		//printf("block %d, raw-score %f, normalize-score %f\n", i, scoresort[i], ns[i]);
	} //normalize the scores

	double	prob0;
       	prob0 = calSigOldRigid(len1, len2, scoresort[0]); //un-normalized score
	
	//printf("block 1, score %f prob %f\n", scoresort[0], prob0);

	double	prob = 1.0;
	r1 = r2 = 1;
	sum = .0;
	for(i = 1; i <= bknum; i ++)	{
		r1 *= i;
		if(i > 1)	r2 *= (i - 1);
		sum += ns[i - 1];
		if(sum < 0)	break;
		spower = 1.0;
		for(k = 1; k <= i - 1; k ++)	{
			spower *= sum;
		}
		cdf = exp(-sum) * spower/(r1 * r2);
		if(cdf > prob)	break; 
		//printf("sum block %d, score %f, r1 %.0f, r2 %.0f, spower %f, cdf %f\n", i, sum, r1, r2, spower, cdf);
		prob = cdf;
	}
	//printf("Final probability %e\n", prob);
	//getchar();
	return prob;
}
