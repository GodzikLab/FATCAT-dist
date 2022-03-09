#include "FatDom.h"

using namespace ARRAY;

//----------------------------------------------------------
FATDOM::FATDOM(char *reference, char *pdbfile, char *alnfile, double pcut)
{
	Initial();

	prot = new PROT(pdbfile);	
	prolen = prot->GetLength();

	Calloc();

	strcpy(proname, reference);
	prot->GetSAA(proseq);

	ReadAlign(reference, alnfile, pcut);
	GetTwistFreq( );

	printf("protein %s: %daa, %d alignments with p-value better than %f, %d alignmnents with twists\n", proname, prolen, pairn, pcut, tpairn);
}

//----------------------------------------------------------
void FATDOM::Initial(void)
{
	prolen = 0;
	proname = NULL;
	prot = NULL;
	pairn = 0;
	tpairn = 0;
	optlen = NULL;
	optpos = NULL;
	count = NULL;
	freq = NULL;
	smooth = NULL;
}

//----------------------------------------------------------
void FATDOM::Calloc(void)
{
	proname = NewArray <char> (100);
	proseq = NewArray <char> (prolen + 1);
	optlen = NewArray <int> (maxpair); //the block information for each alignment
	optpos = NewMatrix <int> (maxpair, maxopt); 
	count = NewArray <int> (prolen);  //the count of twist at each position
	freq = NewArray <double> (prolen);  // the freq of twist at each position
	smooth = NewArray <double> (prolen); //the smoothed frequency of twist at each position
}

//----------------------------------------------------------
FATDOM::~FATDOM(void)
{
	if(prot != NULL)	delete prot;
	if(proseq != NULL)	delete[] proseq;

	if(freq != NULL)	DelArray <double> (freq);
	if(smooth != NULL)	DelArray <double> (smooth);

	if(optlen != NULL)	DelArray <int> (optlen);
	if(optpos != NULL)	DelMatrix <int> (optpos, maxpair);
}

//----------------------------------------------------------
void FATDOM::GetTwistFreq(void)
{
	int	i, j, p;

	//freq the twists
	tpairn = 0;
	for(i = 0; i < pairn; i ++)	{
		//printf("pair %d, optlen %d\n", i, optlen[i]);
		if(optlen[i] <= 1)	continue; //no twist
		for(j = 0; j < optlen[i] - 1; j ++){ //note - 1, terminal is not considered as twist
			p = int((optpos[i][2 * j + 1] + optpos[i][2 * (j + 1)]) / 2);
				//mid-point of the end of current block and begin of next block
			//printf("pair %d  block %d (%d - %d), p1 %d p2 %d\n", i, j, optpos[i][2 * j] + 1, optpos[i][2 * j + 1], p1, p2);
			count[p] += 1;
		}
		tpairn ++; //pair of alignments with twist
	}

	//normalize the twist number by the number of pairs
	for(i = 0; i < prolen; i ++)	{
		freq[i] = double(count[i]) / double(tpairn);
	}

	//smooth the freqing
	int	b1, b2;
	int	window = 3;
	double	add;
	for(i = 0; i < prolen; i ++)	{
		b1 = (i - window) > 0?(i - window):0;
		b2 = (i + window) < (prolen - 1) ? (i + window):(prolen - 1);
		add = 0;
		for(j = b1; j <= b2; j ++)	{
			add += freq[j];	
		}
		smooth[i] = add / double(b2 - b1 + 1);
		//printf("pos %d, twist %.3f, range b1 %d, b2 %d, add %d, ave %f\n", i, freq[i], b1, b2, add, smooth[i]);
	}
}

//get the alignment of reference with other proteins
//the alignment is kept for domain dissection only if 
//   (1) the length of the alignment
//   (2) the similarity significance of the alignment < pcut
//----------------------------------------------------------
void FATDOM::ReadAlign(char *reference, char *alignfile, double pcut)
{
	UniName(reference);
	//printf("reference %s\n", reference);

	ifstream aln(alignfile);
	if(!aln)	{	
		printf("open file %s error\n", alignfile);
		exit(1);
	}
	int	lab = 0;
	int	linen = 0;
	double	prob = 0;
	char	**lines = NewMatrix <char> (1000, 200);
	int	reverse = 0;
	char	str[1000], tmp1[1000], tmp2[1000];
	while(!aln.eof())	{
		str[0] = 0;
		aln.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			sscanf(str, "%*s%s%*s%*s%s", tmp1, tmp2);
			UniName(tmp1);
			UniName(tmp2);
			lab = 1;
			reverse = 0;
			if(!strcmp(tmp1, reference))	{ reverse = 0; }
			else if(!strcmp(tmp2, reference))	{ reverse = 1; }
			else	lab = 0;
			linen = 0;
		}
		else if(!strncmp(str, "Identity", 8))	{
			sscanf(str, "%*s%*s%*s%*s%*s%*s%*s%lf", &prob);
		} //old-fatcat-output-format
		else if(!strncmp(str, "P-value", 7))	{
			sscanf(str, "%*s%lf", &prob);
		} //new-fatcat-output-format
		if(lab)	{
			strcpy(lines[linen ++], str);
		}
		if(lab && !strncmp(str, "Note", 4) && prob <= pcut)	{
			optlen[pairn] = ReadEquPos(linen, lines, optpos[pairn], reverse);
			pairn ++;
		}
	}
	aln.close();


	DelMatrix <char> (lines, 1000);
}

//----------------------------------------------------------
int FATDOM::ReadEquPos(int linen, char **lines, int *pos, int reverse)
{
	int	block = 0;
	char	str[201], str2[201], str3[201], tmp[10];
	int	n = 0;
	int	i;
	int	maxaln = 10000;
	char	*alnseq = NewArray <char> (maxaln);
	char	*blockseq = NewArray <char> (maxaln);
	int	alnbeg1 = 0;
	for(i = 0; i < linen; i ++)	{
		strcpy(str, lines[i]);
		if(!strncmp(str, "Hinges", 6) || !strncmp(str, "Twists", 6))	sscanf(str, "%*s%d", &block);
		else if(!strncmp(str, "Chain 1", 7))	{
			strcpy(str2, lines[++i]);
			strcpy(str3, lines[++i]);
			if(alnbeg1 == 0)	{
				if(reverse)	sscanf(str3, "%*s%*s%s", tmp); //B-A
				else	sscanf(str, "%*s%*s%s", tmp); //A-B
				if(reverse)	strcpy(alnseq, str3 + 14);
				else	strcpy(alnseq, str + 14);
				strcpy(blockseq, str2 + 14);
				alnbeg1 = 1;
			}
			else	{
				if(reverse)	strcat(alnseq, str3 + 14);
				else	strcat(alnseq, str + 14);
				strcat(blockseq, str2 + 14);
			}
		}
	}

	char	*purialnseq = NewArray <char> (maxaln);
	int	purilen = 0;

	char	blockmark = '1';
	int	bg, len, end;
	bg = len = end = 0;
	for(i = 0; i < int(strlen(alnseq)); i ++)	{
		//printf("pos %d %c (total %d)\n", i, alnseq1[i], strlen(alnseq1));
		if(alnseq[i] != '-')	{
			len ++;	
			if(blockseq[i] != ' ')	end = len;
			purialnseq[purilen ++] = alnseq[i];
		}
		if(blockseq[i] != ' ' && blockseq[i] != blockmark)	{ //a new block
			pos[n ++] = bg;  //records the beginning of a block
			pos[n ++] = end; //records the ending of a block
			bg = len; 
			blockmark = blockseq[i];
		}
	}
	pos[n ++] = bg;
	pos[n ++] = end; //note: up to this point, pos records according to the index of sequence from the alignment

	if((block + 1) != n / 2)	{ printf("read alignment incorrectly\n"); exit(1); } 

	for(i = 0; i < prolen; i ++)	{
		if(!strncmp(proseq + i, purialnseq, purilen))	break;
	}
	bg = i;
	if(bg == prolen)	{ printf("alnseq and pdb-seq do not match\n"); exit(1); }
	if(bg != 0)	{
		for(i = 0; i < n; i ++)	pos[i] += bg;
	} //note: rejust pos: pos records the index of block positions according to the sequence from pdb file, index is from 0

	DelArray <char> (alnseq);
	DelArray <char> (purialnseq);
	DelArray <char> (blockseq);

	return (block + 1);
}

//----------------------------------------------------------
void FATDOM::UniName(char *name)
{
	int	len = strlen(name);
	if(!strcmp(name + len - 4, ".pdb"))	{
		name[len - 4] = '\0';
	}
}

//draw the curve
//----------------------------------------------------------
void FATDOM::DrawTwistCurve(char *curvefile)
{
	if(tpairn <= 0)	{ printf("no twists are detected for %s\n", proname); return; }
	int	ptnum = prolen;
	double	**point = NewMatrix <double> (ptnum, 2);
	double	**point2 = NewMatrix <double> (ptnum, 2);
	int	i, idx;
	for(i = 0; i < prolen; i ++)	{
		idx = prot->GetIndex(i); //note: use the pdb position instead of the index
		point[i][0] = idx;
		point[i][1] = smooth[i];
		point2[i][0] = idx;
		point2[i][1] = freq[i];
	}
	PsShow *ps2 = new PsShow(curvefile, "%!PS-Adobe-3.0 EPSF-3.0\n%%BoundingBox: 0 0 550 300","Twist Frequency");
       	ps2->PageArea(50, 450, 80, 200);
	ps2->drawAxis(point2, ptnum, "position", "frequency");
	double	green[3] = {0, 1, 0};
	ps2 -> assignColor(green);
	ps2 -> prtContLine(point2, ptnum, 1);
	double	red[3] = {1, 0, 0};
	ps2 -> assignColor(red);
	ps2 -> prtContLine(point, ptnum, 1);
	char	title[1000];
	sprintf(title, "The schematic representation of twist frequency in each position of %s, raw data in green and smoothed curve in red", proname);
	ps2 -> prtTitle(title);
	delete ps2;
	DelMatrix <double> (point2, ptnum);
	DelMatrix <double> (point, ptnum);
}	

//----------------------------------------------------------
void FATDOM::PrintTwistFreq(char *curvetxt)
{
	if(tpairn <= 0)	{ printf("no twists are detected for %s\n", proname); return; }
	ofstream out2(curvetxt);
	char	str[1000];
	if(!out2)	{
		printf("open %s error\n", curvetxt);
		exit(1);
	}
	sprintf(str, "#%s: %daa, %d alignments with p-value better than the cutoff, %d alignmnents with twists\n", proname, prolen, pairn, tpairn);
	out2<<str;
	out2<<"# the twist frequence at each position for "<<proname<<endl;
	out2<<"# position-index position-in-pdb frequency smoothed-frequency count\n";
	int	i, j;
	for(i = 0; i < prolen; i ++)	{
		j = prot->GetIndex(i);
		sprintf(str, "%d\t%d\t%.3f\t%.3f\t%d\n", i + 1, j, freq[i], smooth[i], count[i]);
		out2<<str;
	}
	out2.close();
}

