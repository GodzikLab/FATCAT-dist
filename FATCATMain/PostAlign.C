//------------------------------------------------------------------
//By Y.Y, 6/11/03
//Latest update: 6/11/03
//construct the complex structure according to the alignment file
//------------------------------------------------------------------

#include "tempkit.h"
#include "geometry.h"
#include "PostAlign.h"
#include "Align0.h"

using namespace ARRAY;
using namespace GEOMETRY;

int colorn = 7;
int colorlist[] = {255,0,0,  205,0,205,  255,255,0, 0,255,0,  0,0,255,  0,255,255, 139,101,8};
                  //red      magenta     yellow     green     blue      cyan       darkgolden

//------------------------------------------------------------------
//Constructors: assign two proteins and extract their AFPs
//------------------------------------------------------------------
POSTALIGN::POSTALIGN(char *file1, char *file2, char *name1, char *name2, char *alnfile)
{
	GetInput(file1, file2, name1, name2, alnfile, "fatcat");
	//default: fatcat-alignment output
}

POSTALIGN::POSTALIGN(char *file1, char *file2, char *name1, char *name2, char *alnfile, char *format)
{
	GetInput(file1, file2, name1, name2, alnfile, format);
}

//------------------------------------------------------------------
void POSTALIGN::GetInput(char *file1, char *file2, char *name1, char *name2, char *alnfile, char *format)
{
	Initial();

	strcpy(pro1Name, name1);
	strcpy(pro2Name, name2);

	ReadPdb(file1, file2);

	Calloc();

	if(!strcmp(format, "fatcat"))	ReadAln(alnfile, name1, name2);
	else if(!strcmp(format, "fasta"))	ReadFastaAln(alnfile, name1, name2);
	else if(!strcmp(format, "clustal"))	ReadClustalAln(alnfile, name1, name2);
	else	{ printf("format %s is not supported currently\n", format); exit(1); }
}

//------------------------------------------------------------------
void POSTALIGN::ReadPdb(char *file1, char *file2)
{
	pro1 = new PROT(file1);
	pro1Len = pro1->GetLength();
	proseq1 = new char[pro1Len + 1];
	pro1->GetSAA(proseq1);

	pro2 = new PROT(file2);
	pro2Len = pro2->GetLength();
	proseq2 = new char[pro2Len + 1];
	pro2->GetSAA(proseq2);

	if(pro1Len <= 0 || pro2Len <= 0)	{
		printf("Wrong pdb read %d %d\n", pro1Len, pro2Len);
		exit(1);
	}

	if(pro1Len < pro2Len)	minLen = pro1Len;
	else	minLen = pro2Len;

}

//------------------------------------------------------------------
//Initialize the class
//------------------------------------------------------------------
void POSTALIGN::Initial(void)
{
	pro1 = NULL;
	pro2 = NULL;
	proseq1 = proseq2 = NULL;
	pro1Len = pro2Len = 0;
	minLen = 0;
	optTwistPdb = NULL;
	focusRes1 = NULL;
	focusRes2 = NULL;
	focusResn = 0;
	totalLenOpt = 0;
	totalRmsdOpt = 0;
	optAln = NULL;
	optLen = NULL;
	optRmsd = NULL;
	maxTra = 5;
	alnseq1 = NULL;
	alnseq2 = NULL;
	alnbeg1 = alnbeg2 = 0;
	hotSpotNum = 0;
	hotSpotIdx = NULL;
	rotateList = NULL;
	translateList = NULL;
	rmsdList = NULL;
	rotateIdx = 0;
	chain1 = 'A';
	chain2 = 'B';
}

//------------------------------------------------------------------
void POSTALIGN::Calloc(void)
{
	if(focusRes1 == NULL)	focusRes1 = NewArray <int> (minLen); 
	if(focusRes2 == NULL)	focusRes2 = NewArray <int> (minLen); 
      	if(optAln == NULL)      {
               	optAln = NewArray3 <int> (maxTra + 1, 2, minLen);
	        optLen = NewArray <int> (maxTra + 1);
	        optRmsd = NewArray <double> (maxTra + 1);
	}
	if(alnseq1 == NULL)	alnseq1 = NewArray <char> (pro1Len + pro2Len + 1);
	if(alnseq2 == NULL)	alnseq2 = NewArray <char> (pro1Len + pro2Len + 1);
}

POSTALIGN::~POSTALIGN(void)
{
	if(pro1 != NULL)	delete pro1;
	if(pro2 != NULL)	delete pro2;
	if(proseq1 != NULL)	delete[] proseq1;
	if(proseq2 != NULL)	delete[] proseq2;
	if(optTwistPdb != NULL)	delete optTwistPdb;
	if(focusRes1 != NULL)	delete[] focusRes1;
	if(focusRes2 != NULL)	delete[] focusRes2;
        if(optAln != NULL)      DelArray3 <int> (optAln, maxTra + 1, 2);
        if(optLen != NULL)      DelArray <int> (optLen); //calloc size maxTra + 1
        if(optRmsd != NULL)     DelArray <double> (optRmsd);
	if(alnseq1 != NULL)	DelArray <char> (alnseq1);
	if(alnseq2 != NULL)	DelArray <char> (alnseq2);
	if(hotSpotIdx != NULL)	DelArray <int> (hotSpotIdx);
	if(rotateList != NULL)	DelArray <double> (rotateList);
	if(translateList != NULL)	DelArray <double> (translateList);
	if(rmsdList != NULL)	DelArray <double> (rmsdList);
}

//------------------------------------------------------------------
void POSTALIGN::WriteAln(char *alnfile, char *name1, char *name2, char *outfile)
{
	ofstream out(outfile);
	if(!out)	{
		printf("open file %s error\n", outfile);
		exit(1);
	}
	WriteAln(alnfile, name1, name2, out);
	out.close();
}

//------------------------------------------------------------------
void POSTALIGN::WriteAln(char *alnfile, char *name1, char *name2, ofstream &out)
{
	char	str[1000];
	ifstream in(alnfile);
	if(!in)	{
		printf("open file %s error\n", alnfile);
		exit(1);
	}
	int	lab = 0;
	char	code1[200], code2[200];
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			sscanf(str, "%*s%s%*s%*s%s", code1, code2);
			UniName(code1);
			UniName(code2);
			if((!strcmp(code1, name1) && !strcmp(code2, name2)) ||
			   (!strcmp(code1, name2) && !strcmp(code2, name1)))	{
				lab = 1;
				out<<str<<endl;
			}
			else	lab = 0;
		}
		else if(lab)	{
			out<<str<<endl;
			if(!strncmp(str, "Note", 4))	break;
		}
	}
	in.close();
}

//hot spot is defined according to pdb 1, using the position in the original pdb instead of index
//------------------------------------------------------------------
void POSTALIGN::GetHotSpot(int num, char **hotspot)
{
	hotSpotNum = num;
	hotSpotIdx = new int[num]; 
	int	i;
	for(i = 0; i < num; i ++)	{
		hotSpotIdx[i] = pro1->ResIndex(hotspot[i]);
		if(hotSpotIdx[i] < 0)	{ printf("can not find hot spot %s\n", hotspot[i]); exit(1); }
	}
}

//------------------------------------------------------------------
void POSTALIGN::ReadAln(char *alnfile, char *name1, char *name2)
{
	UniName(name1);
	UniName(name2);

	char	str[1000];
	ifstream in(alnfile);
	if(!in)	{
		printf("open file %s error\n", alnfile);
		exit(1);
	}
	int	len1, len2, reverse;
	int	lab = 0;
	char	code1[200], code2[200], str2[201], str3[201], blockseq[20000];
	blockNum = reverse = 0;
	alnseq1[0] = alnseq2[0] = blockseq[0] = 0;
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			sscanf(str, "%*s%s%d%*s%s%d", code1, &len1, code2, &len2);
                        UniName(code1);
                        UniName(code2);
			if(!strcmp(code1, name1) && !strcmp(code2, name2))	{ 
				lab = 1;
				reverse = 0;
			}
			else if(!strcmp(code1, name2) && !strcmp(code2, name1))	{
                        	lab = 1;
				reverse = 1;
			}
			else	lab = 0;
		}
		else if(!strncmp(str, "Chain 1", 7) && lab)	{
			in.getline(str2, 200);
			in.getline(str3, 200);
			strcat(blockseq, str2 + 14);
			if(reverse == 0)	{
				strcat(alnseq1, str + 14);
				strcat(alnseq2, str3 + 14);
			}
			else	{
				strcat(alnseq2, str + 14);
				strcat(alnseq1, str3 + 14);
			}
		}
		else if(!strncmp(str, "Note", 4) && lab)	{
			break;
		}
	}
	in.close();

	if(lab == 0)	{ 
		printf("No alignment is found between the given proteins %s %s, in %s\n", name1, name2, alnfile); 
		exit(1); 
	}

	GetAlignPos(blockseq);
}

void POSTALIGN::GetAlignPos(char *blockseq)
{
	char	*alnseq1_puri = Aln2Seq(alnseq1);
	char	*alnseq2_puri = Aln2Seq(alnseq2);

	ALIGN0	*aln1 = new ALIGN0(alnseq1_puri, proseq1, 11, 1);
	ALIGN0	*aln2 = new ALIGN0(alnseq2_puri, proseq2, 11, 1);
	int	*idx1 = aln1->MatchSeqByAlign();
	int	*idx2 = aln2->MatchSeqByAlign();
	//printf("alnseq1 %s\nproseq1 %s\nalnseq2 %s\nproseq2 %s\n", alnseq1_puri, proseq1, alnseq2_puri, proseq2);

	int	i, j;
	char	blockmark = '0';
	int	totalLen;
	int	len1 = -1;
	int	len2 = -1;
	blockNum = -1;
	totalLen = 0;
	for(i = 0; i < int(strlen(alnseq1)); i ++)	{
		if(alnseq1[i] != '-')	{ len1 ++; }
		if(alnseq2[i] != '-')	{ len2 ++; }
		if(hotSpotNum > 0 && blockseq[i] != ' ')	{
			for(j = 0; j < hotSpotNum; j ++)	{
				if(idx1[len1] == hotSpotIdx[j])	{
					blockNum ++;
					optLen[blockNum] = 0;
					break;
				}
			}
		} //new-block forced by the hot-spot
		else if(blockseq[i] != ' ' && blockseq[i] != blockmark)	{ //new-block
			blockNum ++;
			optLen[blockNum] = 0;
			blockmark = blockseq[i];
		}
		if(blockseq[i] != ' ' && alnseq1[i] != '-' && alnseq2[i] != '-')	{
			//optAln[blockNum][0][optLen[blockNum]] = len1 + alnbeg1;
			//optAln[blockNum][1][optLen[blockNum]] = len2 + alnbeg2;
			if(idx1[len1] < 0) { printf("can not find residue %d %c in pdb 1\n", len1, alnseq1[i]); exit(1); }
			if(idx2[len2] < 0) { printf("can not find residue %d %c in pdb 2\n", len2, alnseq2[i]); exit(1); }
			optAln[blockNum][0][optLen[blockNum]] = idx1[len1];
			optAln[blockNum][1][optLen[blockNum]] = idx2[len2];
			//modified by Y.Y, Nov, 20, 2003
			optLen[blockNum] ++;
			totalLen ++;
		}
	}
	blockNum ++;
	if(blockNum >= 1)	{
		rotateList = new double[9 * blockNum];
		translateList = new double[3 * blockNum];
		rmsdList = new double[blockNum];
	}
	
	//printf("blockNum %d, len %d\n", blockNum, totalLen);

	if(totalLen <= 0)	{
		printf("Wrong alignment format\n");
		exit(1);
	}

	delete[] idx1;
	delete[] idx2;
	delete aln1;
	delete aln2;
	delete[] alnseq1_puri;
	delete[] alnseq2_puri;
}

//--------------------------------------------------------
char* POSTALIGN::Aln2Seq(char *alnseq)
{
	int	len = strlen(alnseq);
	char	*seq = new char[len + 1];
	seq[0] = 0;	
	int	i, j;
	j = 0;
	for(i = 0; i < len; i ++)	{
		if(alnseq[i] != '-')	seq[j ++] = alnseq[i];
	}
	seq[j] = 0;
	return seq;
}

//----------------------------------------------------------
void POSTALIGN::UniName(char *name)
{
        int     len = strlen(name);
        if(!strcmp(name + len - 4, ".pdb"))     {
                name[len - 4] = '\0';
        }
}

//--------------------------------------------------------
//calculate the total rmsd of the blocks
//output a merged pdb file for both proteins
//protein 1, in chain A
//protein 2 is twisted accroding to the twists detected, in chain B
//--------------------------------------------------------
void POSTALIGN::TwistPdb(void)
{
	int	i, bk, b2, e2;
	i = bk = b2 = e2 = 0;
	int	*bound = NewArray <int> (blockNum);

	//superimposing according to the optimized alignment 
	if(optTwistPdb == NULL)	{
	       	optTwistPdb = new PROT;
		optTwistPdb->AssignPdb(pro2);
	}
	b2 = 0;
	focusResn = 0;
	for(bk = 0; bk < blockNum; bk ++)	{
		TransPdb(optLen[bk], optAln[bk][0], optAln[bk][1]);
		//transform pro2 according to comparison of pro1 and pro2 at give residues
		if(bk > 0)	{ b2 = e2; }
		if(bk < blockNum - 1)	{ //bend at the middle of two consecutive blocks
			e2 = optAln[bk][1][optLen[bk] - 1]; 
			e2 = (optAln[bk + 1][1][0] - e2)/ 2 + e2;
		}
		else	{ e2 = pro2Len; }
		ModifyCod(optTwistPdb, pro2, b2, e2);
		bound[bk] = e2;
		for(i = 0; i < optLen[bk]; i ++)	{
			focusRes1[focusResn] = optAln[bk][0][i];
			focusRes2[focusResn] = optAln[bk][1][i];
			focusResn ++;
		}
	}
	totalLenOpt = focusResn;
	totalRmsdOpt = pro1->CalCaRmsd(optTwistPdb, focusResn, focusRes1, focusRes2);
	printf("RMSD = %.2f\n", totalRmsdOpt);

	DelArray <int> (bound);
}

//--------------------------------------------------------
void POSTALIGN::WriteTwistPdb(char *pdbfile, char c1, char c2, char *scriptfile)
{
	WriteTwistPdb(pdbfile, c1, c2);

	WriteScript(pro1, optTwistPdb, blockNum, 1, pdbfile, scriptfile);
}

//--------------------------------------------------------
void POSTALIGN::WriteTwistPdb(char *pdbfile, char c1, char c2)
{
	chain1 = c1;
	chain2 = c2;

	TwistPdb();

	char	info[1000];

	//output tiwsted superimposing according to the optimized alignment 
	sprintf(info, "REMARK  protein %s chain A with twisted protein %s in chain B", pro1Name, pro2Name); 
	sprintf(info, "%s\nREMARK  result after optimizing %d blocks %d residues %.3f rmsd",
			info, blockNum, focusResn, totalRmsdOpt);
	WritePdb(info, pdbfile, pro1, optTwistPdb);
}

//---------------------------------------------------------------
void POSTALIGN::ModifyCod(PROT *p1, PROT *p2, int r1, int r2)
{
	int	i, k;
	for(i = r1; i < r2; i ++) {
		for(k = 0; k < 3; k ++)	{
			p1->caCod[3 * i + k] = p2->caCod[3 * i + k];
		}
	} //modify caCod
	int	a2 = p2->totalAtm;
	if(r2 < p2->GetLength())	{ a2 = p2->res2Atm[r2]; }
	for(i = p2->res2Atm[r1]; i < a2; i ++) {
		for(k = 0; k < 3; k ++)	{
			p1->atmCod[3 * i + k] = p2->atmCod[3 * i + k];
		}
	} //modify atmCod
}

//--------------------------------------------------------
//write the pdb according to the superimposing of the given position pairs
//--------------------------------------------------------
void POSTALIGN::TransPdb(int n, int *res1, int *res2)
{
	double	*cod1 = pro1->Cod4Res(n, res1);
	double	*cod2 = pro2->Cod4Res(n, res2);
	double	r[9], t[3], e[3];
	double	rmsd = kearsay(n, cod1, cod2, r, t);
	rot2euler(r, e);
	DelArray <double> (cod1);
	DelArray <double> (cod2);

	//record the transformation matrix for each block
	if(rotateIdx >= blockNum)	{ printf("rotate number error\n"); exit(1); }
	rmsdList[rotateIdx] = rmsd;
	int	i;
	for(i = 0; i < 3; i ++)	{
		translateList[3 * rotateIdx + i] = t[i]; 
	}
	for(i = 0; i < 9; i ++)	{
		rotateList[9 * rotateIdx + i] = r[i];
	}
	rotateIdx ++;

	//transformed the coordinations of whole protein 2 accordingly
	pro2->TransCod(r, t); //refer Prot.C
}

//--------------------------------------------------------
//write the superimposed pdb
//--------------------------------------------------------
void POSTALIGN::WritePdb(char *info, char *pdbfile)
{
	WritePdb(info, pdbfile, pro1, pro2);
}

//write two proteins into one file
//--------------------------------------------------------
void POSTALIGN::WritePdb(char *info, char *pdbfile, PROT *p1, PROT *p2)
{
	ofstream pdb(pdbfile);
	if(!pdb) { printf("open pdbfile %serror\n", pdbfile); exit(1); }
	pdb<<"REMARK  the superimposed protein file"<<endl;
	pdb<<info<<endl;
	pdb<<"REMARK  "<<pro1Name<<" chain "<<chain1<<endl;
	pdb<<"REMARK  "<<pro2Name<<" chain "<<chain2<<endl;
	p1->PrintPdb(chain1, pdb);
	pdb<<"TER"<<endl;
	p2->PrintPdb(chain2, pdb);
	pdb<<"END"<<endl;
	pdb.close();
	
}

//--------------------------------------------------------
//write the rasmol script for superimposing a block 
//--------------------------------------------------------
void POSTALIGN::WriteScript(PROT *p1, PROT *p2, int bk, int ifopt, char *pdbfile, char *scriptfile)
{
	char	region[5000] = (" ");
	int	i, b1, b2, e1, e2;	
	b1 = b2 = e1 = e2 = 0;
	for(i = 0; i < blockNum; i ++)	{
		if(bk != blockNum && i != bk)	continue; //if bk == blockNum, all blocks
		if(ifopt)	{
			b1 = optAln[i][0][0];
			b2 = optAln[i][1][0];
			e1 = optAln[i][0][optLen[i] - 1];
			e2 = optAln[i][1][optLen[i] - 1];
		}
		sprintf(region, "%s\nselect %s-%sA\ncolor [181,181,181]", region, p1->index[b1], p1->index[e1]);//gray
		sprintf(region, "%s\nselect %s-%sB\ncolor [%d,%d,%d]\n", region, p2->index[b2], p2->index[e2], 
				colorlist[3 * i], colorlist[3 * i + 1], colorlist[3 * i + 2]); //different color
	}
	ofstream script(scriptfile);
	script<<"#!rasmol -script\n# File: "<<scriptfile<<"\n# Creator: FATCAT\n"<<endl;
	script<<"zap\nload pdb \""<<pdbfile<<"\""<<endl;
	script<<"background white\nset ambient 40\nset specular on\n";
	script<<"wireframe off\n\nselect all\ncartoon\n"<<endl;
	script<<"select *"<<chain1<<"\ncolor [219,219,219]\nselect *"<<chain2<<"\ncolor [232,232,232]\n"<<endl;
	                  //light gray         //light blue 
	script<<region<<endl;
	script.close();
}

//------------------------------------------------------------------
//Show AFP chain in postscript format 
//------------------------------------------------------------------
void POSTALIGN::ShowAfp(char *outfile, int ifcolor)
{
	PsShow	*psShow = new PsShow(outfile, "%!PS-Adobe-3.0 EPSF-3.0\n%%BoundingBox: 0 0 500 500", "alignment graph from FATCAT");	//refer PsShow.C and PsShow.h
	psShow->PageArea(50, 400, 50, 400);
	psShow->setupsidedown();
	psShow->setMarkFontSize(12);
	psShow->setTextFontSize(12);


	//extract the AFP-chain coordinations
	int	alnLength = strlen(alnseq1);
	int	ptnum = 2 * alnLength;
        double  **point = NewMatrix <double> (ptnum, 2); //refer myArrayTemp.h
        double  *gray = NewArray <double> (ptnum);
        double  *linewidth = NewArray <double> (ptnum);
        double  **color = NewMatrix <double> (ptnum, 3);

	int	i, k, real_ptnum, notgap, notgapp, a, b, ap, bp, ta, tb, block, g;
	a = ap = alnbeg1;
	b = bp  = alnbeg2;
	ta = tb = -1;
	block = notgapp = 0;
	g = 0;
	for(i = 0; i < alnLength; i ++)	{
		if(alnseq1[i] == '-' || alnseq2[i] == '-')	notgap = 0;	
		else	notgap = 1;
		if((i == alnLength - 1 || alnseq1[i + 1] == '-' || alnseq2[i + 1] == '-') && notgap == 1)	{
			if(g > ptnum)	{
				printf("please increase ptnum %d\n", ptnum);
			}
			if(ta != -1)	{
				linewidth[2 * g] = 0.5;
				gray[2 * g] = 0.0;
				point[2 * g][0] = ta + 1;
				point[2 * g][1] = tb + 1;
				point[2 * g + 1][0] = ap + 1;
				point[2 * g + 1][1] = bp + 1;
				for(k = 0; k < 3; k ++)	{
					if(block > colorn - 1)	color[2 * g][k] = colorlist[3 * (colorn - 1) + k];
					else	color[2 * g][k] = colorlist[3 * block + k];
				}
				g ++;
			} //connection between AFP
			linewidth[2 * g] = 2;
			gray[2 * g] = 0.0;
			point[2 * g][0] = ap + 1;
			point[2 * g][1] = bp + 1;
			point[2 * g + 1][0] = a + 1;
			point[2 * g + 1][1] = b + 1;
			for(k = 0; k < 3; k ++)	{
				if(block > colorn - 1)	color[2 * g][k] = colorlist[3 * (colorn - 1) + k];
				else	color[2 * g][k] = colorlist[3 * block + k];
			}
			g ++;
			//AFPs
			notgap = 0;
			ta = a;
			tb = b;
		}
		if(notgapp == 0 && notgap == 1)	{
			ap = a;
			bp = b;
		}
		if(alnseq1[i] != '-')	a ++;
		if(alnseq2[i] != '-')	b ++;
		for(k = 0; k < blockNum; k ++)	{
			if(optAln[k][0][0] == a)	{
				block = k;
			}
		}
		notgapp = notgap;
	}
	real_ptnum = 2 * g;

	//draw the AFP-chain in color or black-and-white
	int     ifcodtransform = 0;

        psShow->drawAxis(point, real_ptnum, pro1Name, pro2Name);
      	psShow->setTransform(); //set the scale & translate

	if(ifcolor)	{
		psShow->prtColorStepLine(point, real_ptnum, color, linewidth, ifcodtransform);
	}
	else	{
		psShow->prtGrayStepLine(point, real_ptnum, gray, linewidth, ifcodtransform);
	}

	//draw the arrows in the twist position
	double	**arrowpoint = NewMatrix <double> (blockNum - 1, 2);
	for(i = 1; i < blockNum; i ++)	{
		arrowpoint[i - 1][0] = optAln[i][0][0];
		arrowpoint[i - 1][1] = optAln[i][1][0];
	}
	psShow->prtArrow(arrowpoint, blockNum - 1, ifcodtransform, "downleftarrow");
	DelMatrix <double> (arrowpoint, blockNum - 1);

	DelArray <double> (linewidth);
	DelArray <double> (gray);
	DelMatrix <double> (point, ptnum);
	DelMatrix <double> (color, ptnum);

	delete psShow;
}

void POSTALIGN::WriteMatrix(char *outfile)
{
	if(optTwistPdb == NULL) { TwistPdb(); } //YY Jul15,06

	printf("blockNum %d\n", blockNum);
	ofstream out(outfile);
	if(!out)	{ printf("open %s error\n", outfile); exit(1); }
	if(blockNum == 1)	{
		out<<"No twist found between "<<pro1Name<<" and "<<pro2Name<<endl;
		//out.close();
		//return;
	}
	out<<"#The transformation matrix for "<<blockNum<<" blocks"<<endl;
	int	i, k;
	char	str[1000];
	for(i = 0; i < blockNum; i ++)	{
		sprintf(str, "#block %d rmsd %.2f\nrotate ", i + 1, rmsdList[i]);
		out<<str;
		for(k = 0; k < 9; k ++)	{
			if(k == 8)	sprintf(str, "%.3f\n", rotateList[9 * i + k]);
			else	sprintf(str, "%.3f ", rotateList[9 * i + k]);
			out<<str;
		}
		sprintf(str, "translate %.3f %.3f %.3f\n", translateList[3*i], translateList[3*i+1], translateList[3*i+2]); 
		out<<str;
		/*
		sprintf(str, "%.3f %.3f %.3f %.3f %.3f %.3f # %.3f\n", 
			translateList[3*i], translateList[3*i+1], translateList[3*i+2],
			rotateList[3*i], rotateList[3*i+1], rotateList[3*i+2], rmsdList[i]);
		out<<str;
		*/
	}
	out.close();
}

void POSTALIGN::ReadFastaAln(char *file, char *name1, char *name2)
{
	//read alignment file in fasta format
	ifstream in(file);
	if(!in)	{ printf("open %s error\n", file); exit(1); }
	int	lab = 0;
	alnseq1[0] = alnseq2[0] = 0;
	char	str[1001], tmp[1000];
	while(!in.eof())	{
		in.getline(str, 1000);
		if(str[0] == '>')	{ 
			sscanf(str + 1, "%s", tmp);
			if(!strcmp(tmp, name1))	{	
				lab = 1;
			}
			else if(!strcmp(tmp, name2))	{
				lab = 2;
			}
			else	lab = 0;
		}
		else if(strlen(str) >= 1)	{
			if(lab == 1)	{
				strcat(alnseq1, str);
			}
			else if(lab == 2)	{
				strcat(alnseq2, str);
			}
		}
	}
	in.close();

	if(alnseq1[0] == 0)	{ printf("can not find %s sequence in the alignment file or the format is not correct\n", name1); exit(1); }
	if(alnseq2[0] == 0)	{ printf("can not find %s sequence in the alignment file or the format is not correct\n", name2); exit(1); }

	int	len = strlen(alnseq1);
	char	*blockseq = new char[len + 1];
	int	i;
	for(i = 0; i < len; i ++)	blockseq[i] = '1'; //in this case, no block information 

	GetAlignPos(blockseq);

	delete[] blockseq;
}

void POSTALIGN::ReadClustalAln(char *file, char *name1, char *name2)
{
	//read alignment file in clustal format
	ifstream in(file);
	if(!in)	{ printf("open %s error\n", file); exit(1); }
	alnseq1[0] = alnseq2[0] = 0;
	char	str[1001], tmp[1000], tmp2[1000];
	while(!in.eof())	{
		in.getline(str, 1000);
		if(strlen(str) > 3)	{
			sscanf(str, "%s", tmp);
			if(!strcmp(tmp, name1))	{	
				sscanf(str, "%*s%s", tmp2);
				strcat(alnseq1, tmp2);
			}
			else if(!strcmp(tmp, name2))	{
				sscanf(str, "%*s%s", tmp2);
				strcat(alnseq2, tmp2);
			}
		}
	}
	in.close();

	if(alnseq1[0] == 0)	{ printf("can not find %s sequence in the alignment file or the format is not correct\n", name1); exit(1); }
	if(alnseq2[0] == 0)	{ printf("can not find %s sequence in the alignment file or the format is not correct\n", name2); exit(1); }

	int	len = strlen(alnseq1);
	char	*blockseq = new char[len + 1];
	int	i;
	for(i = 0; i < len; i ++)	blockseq[i] = '1'; //in this case, no block information 

	GetAlignPos(blockseq);

	delete[] blockseq;
}
