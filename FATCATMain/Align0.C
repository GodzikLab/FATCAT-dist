#include "Align0.h"
#include "Amino.h"

using namespace ARRAY;
using namespace AMINO;

//-----------------------------------------------------------------------------
//given matrix
//-----------------------------------------------------------------------------
ALIGN0::ALIGN0(double **sij0, int M0, int N0, double g0, double h0)
{
	Initial(M0, N0, g0, h0);
	int	i, j;
	for(i = 0; i < M; i ++)	{
		for(j = 0; j < N; j ++)	sij[i][j] = sij0[i][j];
	}
	DoAlign();
}

//-----------------------------------------------------------------------------
//given protein sequences
//-----------------------------------------------------------------------------
ALIGN0::ALIGN0(char *seq10, char *seq20, double g0, double h0)
{
	int	M0 = strlen(seq10);
	int	N0 = strlen(seq20);

	Initial(M0, N0, g0, h0);
	strcpy(seq1, seq10);
	strcpy(seq2, seq20);

	int	i, j;
	for(i = 0; i < M0; i ++)	{
		for(j = 0; j < N0; j ++)	{
			sij[i][j] = aaScore(seq1[i], seq2[j]); //refer Amino.C
		}
	}

	DoAlign();
}

//-----------------------------------------------------------------------------
//given alignment sequences
//-----------------------------------------------------------------------------
ALIGN0::ALIGN0(char *seq10, char *seq20)
{
	int	M0 = strlen(seq10);
	int	N0 = strlen(seq20);
	if(M0 != N0)	{
		printf("input are not alignments\n");
		exit(1);
	}

	Initial(M0, N0, 0, 0);
	strcpy(seq1, seq10);
	strcpy(seq2, seq20);

	identity = 0.0;
	similarity = 0.0;

	int	i;
	for(i = 0; i < M0; i ++)	{
		if(seq1[i] == seq2[i])	identity += 1.0;	
		else if(seq1[i] == '-' || seq1[i] == '*' || seq1[i] == '.' ||
		   seq2[i] == '-' || seq2[i] == '*' || seq2[i] == '.' )	
			continue;
		else if(aaScore(seq1[i], seq2[i]) > 0)	similarity += 1.0;
	}

	similarity = double(identity + similarity) / double(M0);
	identity /= double(M0);

}

//-----------------------------------------------------------------------------
void ALIGN0::Initial(int M0, int N0, double g0, double h0)
{
	M = M0;
	N = N0;
	g = g0; //gap-create
	h = h0; //gap-extend
	m = g + h; //gap-create + gap-extend
	trace = NewMatrix <char> (M + 1, N + 1);
	etrace = NewMatrix <char> (M + 1, N + 1);
	dtrace = NewMatrix <char> (M + 1, N + 1);
	B1 = B2 = E1 = E2 = 0;
	alignScore = 0;
	last = 0;
	sapp = NewArray <int> (M + N);
	sapp0 = sapp;
	sij = NewMatrix <double> (M, N);
	seq1 = NewArray <char> (M + 1);
	seq2 = NewArray <char> (N + 1);
	aln1 = aln2 = mark = NULL;
	identity = similarity = 0;
}

//-----------------------------------------------------------------------------
ALIGN0::~ALIGN0(void)
{
	if(trace != NULL)	DelMatrix <char> (trace, M + 1);
	if(etrace != NULL)	DelMatrix <char> (etrace, M + 1);
	if(dtrace != NULL)	DelMatrix <char> (dtrace, M + 1);
	if(sapp0 != NULL)	DelArray <int> (sapp0);
	if(sij != NULL)		DelMatrix <double> (sij, M);
	if(seq1 != NULL)	DelArray <char> (seq1);
	if(seq2 != NULL)	DelArray <char> (seq2);
	if(aln1 != NULL)	DelArray <char> (aln1);
	if(aln2 != NULL)	DelArray <char> (aln2);
	if(mark != NULL)	DelArray <char> (mark);
}

//-----------------------------------------------------------------------------
//trace-back strategy
//affine-gap penalty
//local-model
//-----------------------------------------------------------------------------
void ALIGN0::DoAlign(void)
{
	int	i, j;
	double	s, e, c, d, wa;
	double	*CC = NewArray <double> (N + 1); //note N + 1
	double	*DD = NewArray <double> (N + 1);
	double	maxs = -100;
	char	trace_e, trace_d;

	//forward-phase
	CC[0] = 0;
	for(j = 1; j <= N; j ++)	{
		CC[j] = 0;
		DD[j] = -g; 
	} //local-alignment, no terminal penalty
	for(i = 1; i <= M; i ++)	{
		CC[0] = c = s = 0; 
		e = -g;
		for(j = 1; j <= N; j ++)	{
			trace_e = 'e';
			if ((c =   c   - m) > (e =   e   - h)) {
				e = c;  trace_e = 'E';
			}//insertion
			trace_d = 'd';
			if ((c = CC[j] - m) > (d = DD[j] - h)) {
				d = c;  trace_d = 'D';
			}//deletion
			//ie   CC[j]==CC[i-1][j]   DD[j]==DD[i-1][j]
			wa = sij[i - 1][j - 1]; //note i - 1, j - 1
			c = s + wa; //s==CC[i-1][j-1]
			trace[i][j] = 's';
			if (e > c) {
				c = e;
				trace[i][j] = trace_e;
			}
			if (d > c) {
				c = d;
				trace[i][j] = trace_d;
			}
			etrace[i][j] = trace_e;
			dtrace[i][j] = trace_d;
			s = CC[j]; //important for next replace
			CC[j] = c; //CC[i][j]
			DD[j] = d; //DD[i][j]
			if(c < 0)	{
				CC[j] = 0;
				DD[j] = -g;
				c = 0;
				e = -g; 
				trace[i][j] = '0';
			} //local-N
			if(c > maxs)	{
				E1 = i;
				E2 = j;
				maxs = c;
			} //local-C
		}
	}
	alignScore = maxs;
	//printf("alignment score %f\n", alignScore);

	DelArray <double> (CC);
	DelArray <double> (DD);

	//trace-back
	if(trace[E1][E2] != 's')	{
		printf("Not end with substitution\n");
		exit(1);
	}
	//Trace(maxs, E1, E2);
	Trace('s', E1, E2);
	//printf("B1 %d B2 %d, E1 %d E2 %d\n", B1, B2, E1, E2);

	//check-alignment
	CheckAlign();
}

//-----------------------------------------------------------------------------
//trace-back, recorded in sapp, wrong method!
//-----------------------------------------------------------------------------
void ALIGN0::Trace(char mod, int i, int j)
{
	if(mod == '0' || i <= 0 || j <= 0)	{
		B1 = i + 1;
		B2 = j + 1;
	}
	if(mod == 's')	{
		Trace(trace[i - 1][j - 1], i - 1, j - 1);
		REP();
	}
	else if(mod == 'D')	{
		Trace(trace[i - 1][j], i - 1, j);
		DEL(1);
	}
	else if(mod == 'd')	{
		Trace(dtrace[i - 1][j], i - 1, j);
		DEL(1);
	}
	else if(mod == 'E')	{
		Trace(trace[i][j - 1], i, j - 1);
		INS(1);
	}
	else if(mod == 'e')	{
		Trace(etrace[i][j - 1], i, j - 1);
		INS(1);
	}
}

//-----------------------------------------------------------------------------
//record the alignment in sapp
//deletion, sapp < 0, sequences in i, gaps in j
//-----------------------------------------------------------------------------
void ALIGN0::DEL(int k)
{
        if(last < 0)    last = sapp[-1] -= (k);
        else            last = *sapp++ = -(k);
}

//Insertion, sapp > 0, gaps in i, sequences in j
//-----------------------------------------------------------------------------
void ALIGN0::INS(int k)
{
	//if(last < 0) { sapp[-1] = (k); *sapp++ = last; }
        //else    last = *sapp++ = (k);
	if(last > 0) 	last = sapp[-1] += k; 
        else    	last = *sapp++ = (k);
}

//-----------------------------------------------------------------------------
void ALIGN0::REP(void)
{
	last = *sapp++ = 0;
}

//-----------------------------------------------------------------------------
void ALIGN0::FastaAlign(void)
{
	if(aln1 == NULL)	Display();
	printf(">%s\n%s\n>%s\n%s\n", "seq1", aln1, "seq2", aln2);
}

//-----------------------------------------------------------------------------
void ALIGN0::MapAlign(void)
{
	if(aln1 == NULL)	Display();
	int	*index = MatchSeqByAlign();
	int	i;
	for(i = 0; i < M; i ++)	{
		if(index[i] == -1)	continue;
		printf("%d %d\n", i + 1, index[i] + 1);
	}
	delete[] index;
}

//-----------------------------------------------------------------------------
void ALIGN0::PlainAlign(void)
{
	if(aln1 == NULL)	Display();
	char	a[100], b[100], c[100];
	int	n = 0;
	int	m, k;
	int	i = B1;
	int	j = B2;
	while(n < int(strlen(mark)))	{
		if((strlen(mark) - n) < 70)	m = strlen(mark) - n;
		else	m = 70;
		a[0] = b[0] = c[0] = 0;
		strncpy(a, aln1, m); 
		strncpy(b, aln2, m); 
		strncpy(c, mark, m); 
		a[m] = b[m] = c[m] = 0;
		printf("%-10s %3d %s\n", "Chain 1:", i, a);
		printf("%-10s %3s %s\n", " ", " ", c);
		printf("%-10s %3d %s\n\n", "Chain 2:", j, b);
		for(k = 0; k < m; k ++)	{
			if(a[k] != '-')	i ++;
			if(b[k] != '-')	j ++;
		}
		n += m;
	}
}

//-----------------------------------------------------------------------------
void ALIGN0::Display(void)
{
	Display(seq1, seq2);
}

//Note: A, B recorded from 0, while B1, B2, E1, E2 from 1
//-----------------------------------------------------------------------------
void ALIGN0::Display(char *A, char *B)
{
	int	i = B1;
	int	j = B2;
	int	s = 0;
	char	*a = NewArray <char> (M + N + 1);
	char	*b = NewArray <char> (M + N + 1);
	char	*c = NewArray <char> (M + N + 1);
	int	n1, n2, n3, op, k;
	n1 = n2 = n3 = 0;

	a[0] = b[0] = c[0] = 0;
  	while (i <=  E1 && j <= E2) {
		op = sapp0[s ++];
		if(op == 0)	{
			a[n1] = A[i - 1];
			b[n2] = B[j - 1];
			c[n3 ++] = (a[n1 ++] == b[n2 ++])?'|':' ';
			i ++;
			j ++;
		}
		else if(op > 0)	{
			for(k = 0; k < op; k ++)	{
	            		a[n1 ++] = '-';
				b[n2 ++] = B[j - 1];
				c[n3 ++] = ' ';
				j ++;
			}
		}
          	else	{
			for(k = 0; k < -op; k ++)	{
	            		a[n1 ++] = A[i - 1];
				b[n2 ++] = '-';
				c[n3 ++] = ' ';
				i ++;
			}
		}
        }
	a[n1] = b[n2] = c[n3] = '\0';
	int	len = strlen(a);
	if(aln1 == NULL)	aln1 = NewArray <char> (len + 1);
	if(aln2 == NULL)	aln2 = NewArray <char> (len + 1);
	if(mark == NULL)	mark = NewArray <char> (len + 1);
	strcpy(aln1, a);
	strcpy(aln2, b);
	strcpy(mark, c);
	DelArray <char> (a);
	DelArray <char> (b);
	DelArray <char> (c);
}

//-----------------------------------------------------------------------------
void ALIGN0::CheckAlign(void)
{
	if(sapp[0] != 0)	{
		printf("warn: not a local-alignment result, first operation %d\n", sapp[0]);
	}
	double	sco = CheckScore();
	if(fabs(sco - alignScore) > 1e-3)	{
		printf("warn: alignment scores are different %f(check) %f(align)\n", sco, alignScore); 	
	}	
}

//-----------------------------------------------------------------------------
// checkscore - return the score of the alignment stored in sapp
//-----------------------------------------------------------------------------
double ALIGN0::CheckScore(void)
{ 
  	int   	i, j, op, s;
  	double 	sco;

  	sco = 0;
	op = 0;
	s = 0;

	i = B1;
	j = B2;
  	while (i <= E1 && j <= E2) {
		op = sapp0[s ++];
		if (op == 0) 	{
			sco += sij[i - 1][j - 1];
			//printf("%d-%d %f\n", i - 1, j - 1, sij[i - 1][j - 1]);
			i ++;
			j ++;
		} 
		else if (op > 0) {
			sco -= g+op*h;
			j = j+op;
		} 
		else {
			sco -= g-op*h;
			i = i-op;
		}
	}
  	return(sco);
}

//record the aligned pairs in alignList[][0], alignList[][1];
//return the number of aligned pairs
//-----------------------------------------------------------------------------
int ALIGN0::GetAlignPos(int **alignList)
{
	int	i = B1;
	int	j = B2;
	int	s = 0;
	int	a = 0;
	int	op;
	while(i <= E1 && j <= E2)	{
		op = sapp0[s ++];
		if (op == 0) 	{
			alignList[0][a] = i - 1; //i - 1
			alignList[1][a] = j - 1;
			a ++;
			i ++;
			j ++;
		} 
		else if (op > 0) {
			j += op;
		} 
		else {
			i -= op;
		}
	}
	return a;
}

//map the aligned positions from seq 1 to seq 2
//return the index list, 
//index = -1 for position i means i in protein 1 is not aligned to any residues from pro2
//-----------------------------------------------------------------------------
int* ALIGN0::MatchSeqByAlign(void)
{
	int     *index = NewArray <int> (M);
	for(int k = 0; k < M; k ++)      index[k] = -1;
	int	i = B1;
	int	j = B2;
	int	s = 0;
	int	op;
	while(i <= E1 && j <= E2)	{
		op = sapp0[s ++];
		if (op == 0) 	{
			index[i - 1] = j - 1;
			i ++;
			j ++;
		} 
		else if (op > 0) {
			j += op;
		} 
		else {
			i -= op;
		}
	}
	return index;
}

double ALIGN0::GetIdentity(void)
{
	return identity;
}

double ALIGN0::GetSimilarity(void)
{
	return similarity;
}
