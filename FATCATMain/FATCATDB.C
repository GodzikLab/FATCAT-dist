#include "FATCATDB.h"
#include "HtmlReport.h"
#include "SigEva.h"
#define MAX 70000  /*increased from 50000 to 70000*/

using namespace ARRAY;
using namespace STATIS;
/*
 a post-processing of the structure database searching by FATCAT
 1. collect the alignment-result 
 2. draw a score-distribution graph in postscript
 3. output potential candidates with high score
*/

FATCATDB::FATCATDB(int sort_in, int cut_in)
{
	inputFile[0] = 0;
	itemNum = 0;
	itemvalue = NULL;
	code1 = NULL; 
	code2 =  NULL;
	len1 =  NULL;
	len2 =  NULL;
	twi =  NULL;
	iniLen =  NULL;
	optLen =  NULL;
	iniRms =  NULL;
	optRms =  NULL;
	chainRms =  NULL;
	score =  NULL;
	probability = NULL;
	alnLen =  NULL;
	gap =  NULL;
	bkNum = NULL;
	bkScore = NULL;
	scalelist = NULL;
	distrib = NULL;
	sort = sort_in;
	sortorder = NULL;
	sortscore = NULL;
	cut = cut_in;

	Calloc(MAX);  //calloc the variants initially
}

void FATCATDB::Calloc(int num0)
{
	itemNum0 = num0; //calloc-number, not the real alignment number 
	itemNum = 0; //real alignment number
	code1 = NewMatrix <char> (num0, 50);
	code2 = NewMatrix <char> (num0, 50);
	len1 = NewArray <int> (num0);
	len2 = NewArray <int> (num0);
	twi = NewArray <int> (num0);
	iniLen = NewArray <int> (num0);
	optLen = NewArray <int> (num0);
	iniRms = NewArray <double> (num0);
	optRms = NewArray <double> (num0);
	chainRms = NewArray <double> (num0);
	score = NewArray <double> (num0);
	probability = NewArray <double> (num0);
	alnLen = NewArray <int> (num0);
	gap = NewArray <double> (num0);
	bkNum = NewArray <int> (num0);
	bkScore = NewMatrix <double> (num0, 10); //maximum 10 blocks
}

FATCATDB::~FATCATDB(void)
{
	if(code1 != NULL)	DelMatrix <char> (code1, itemNum0);
	if(code2 != NULL)	DelMatrix <char> (code2, itemNum0);
	if(len1 != NULL)	DelArray <int> (len1);
	if(len2 != NULL)	DelArray <int> (len2);
	if(twi != NULL)	DelArray <int> (twi);
	if(iniLen != NULL)	DelArray <int> (iniLen);
	if(optLen != NULL)	DelArray <int> (optLen);
	if(iniRms != NULL)	DelArray <double> (iniRms);
	if(optRms != NULL)	DelArray <double> (optRms);
	if(chainRms != NULL)	DelArray <double> (chainRms);
	if(score != NULL)	DelArray <double> (score);
	if(probability != NULL)	DelArray <double> (probability);
	if(alnLen != NULL)	DelArray <int> (alnLen);
	if(gap != NULL)	DelArray <double> (gap);
	if(bkNum != NULL)	DelArray <int> (bkNum);
	if(bkScore != NULL)	DelMatrix <double> (bkScore, itemNum0);
	if(distrib != NULL)	DelArray <double> (distrib);
	if(scalelist != NULL)	DelArray <double> (scalelist);
	if(itemvalue != NULL)	DelArray <double> (itemvalue);
}

//read database searching result from one line per one pair alignment result
void FATCATDB::ReadReport(char *alnfile)
{
	strcpy(inputFile, alnfile);
	ifstream in(alnfile);
	if(!in)	{
		printf("open file %s error\n", alnfile);
		exit(1);
	}
	char	str[1000];
	int	n = 0;
	while(!in.eof())	{
		in.getline(str, 200);
		if(str[0] == '#' || str[0] == '\n')	continue;
		sscanf(str, "%s%s%d%d%d%d%lf%d%lf%lf%lf%d%lf", 
			code1[n], code2[n], &len1[n], &len2[n], &twi[n], &iniLen[n], &iniRms[n], 
			&optLen[n], &optRms[n], &chainRms[n], &score[n], &alnLen[n], &gap[n]);
		if(len1[n] < cut || len2[n] < cut)	continue;
		//if(optLen[n] == 0)	gap[n] = 0;
		//else	gap[n] /= (gap[n] + optLen[n]); //percentage
		n ++;
	}
	itemNum = n;
	in.close();
}

//read alignments for a list of pairs
void FATCATDB::ReadAlignList(char *listfile, char *dir)
{
	ifstream in(listfile);
	if(!in)	{
		printf("open listfile %s error\n", listfile);
		exit(1);
	}
	char	alnfile[200], tmp1[100], tmp2[100], str[201];
	int	itemold = 0;
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(str[0] == '#' || strlen(str) < 2)	continue;
		sscanf(str, "%s%s", tmp1, tmp2); 
		sprintf(alnfile, "%s/%s.%s.aln", dir, tmp1, tmp2);
		ReadAlign(alnfile);
		if(itemNum > itemold)	{
			strcpy(code1[itemNum - 1], tmp1);
			strcpy(code2[itemNum - 1], tmp2);
		}
		itemold = itemNum;
	}
	in.close();
	printf("Read %d pairs\n", itemNum);
}

//read database searching result from detailed alignment output
void FATCATDB::ReadAlign(char *alnfile)
{
	strcpy(inputFile, alnfile);
	char	str[1000], tmp[100];
	ifstream in(alnfile);
	if(!in)	{
		printf("open file %s error\n", alnfile);
		exit(1);
	}
	int	i;
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			if(itemNum >= MAX)	{
				printf("Please, increase MAX %d\n", MAX);
				exit(1);
			}
			sscanf(str, "%*s%s%d%*s%s%d", code1[itemNum], &len1[itemNum], code2[itemNum], &len2[itemNum]);
			//printf("align %s %s\n", code1[itemNum], code2[itemNum]);
			bkNum[itemNum] = 0;
		}
		else if(!strncmp(str, "Hinges", 6) || !strncmp(str, "Twists", 6))	{
			sscanf(str, "%*s%d%*s%d%*s%lf%*s%d%*s%lf%*s%lf%*s%lf%*s%d%*s%s",
				&twi[itemNum], &iniLen[itemNum], &iniRms[itemNum], &optLen[itemNum], &optRms[itemNum], &chainRms[itemNum], &score[itemNum], &alnLen[itemNum], tmp); 
			for(i = 0; i < int(strlen(tmp)); i ++)	if(tmp[i] == '(')	tmp[i] = ' ';
			sscanf(tmp, "%lf", &gap[itemNum]);

		}
		else if(!strncmp(str, "P-value", 7))	{
			sscanf(str, "%*s%lf", &probability[itemNum]);
		}
		else if(!strncmp(str, "Block", 5))	{
			sscanf(str, "%*s%*s%*s%*s%*s%lf", &bkScore[itemNum][bkNum[itemNum] ++]);
			//printf("block %s, %f\n", str, bkScore[n][bkNum[n] - 1]);
			//getchar();
		}
		else if(!strncmp(str, "Note", 4))	{
			if(len1[itemNum] >= cut && len2[itemNum] >= cut)	itemNum ++;
		}
	}
	in.close();
}

//sort alignment by probability
void FATCATDB::SortByProb(void)
{
	if(sortorder != NULL)	return; 
	sortorder = NewArray <int> (itemNum);
	sortscore = NewArray <double> (itemNum);

	int	i;
	for(i = 0; i < itemNum; i ++)	sortorder[i] = i;

	if(itemvalue == NULL)	{
		printf("GetTarget must have been excuted\n");
		exit(1);
	}

	//GetTarget("fprobability", probability);

	for(i = 0; i < itemNum; i ++)	sortscore[i] = probability[i];

	qksort <double> (itemNum, sortscore, sortorder);
}

//sort alignment by probability
//then by rmsd & alignment-length-rmsd
void FATCATDB::SortByProbAdOther(void)
{
	if(sortorder != NULL)	return; 
	sortorder = NewArray <int> (itemNum);
	sortscore = NewArray <double> (itemNum);

	int	i;
	for(i = 0; i < itemNum; i ++)	sortorder[i] = i;

	if(itemvalue == NULL)	{
		printf("GetTarget must have been excuted\n");
		exit(1);
	}

	//GetTarget("fprobability", probability);

	for(i = 0; i < itemNum; i ++)	sortscore[i] = probability[i];

	qksort <double> (itemNum, sortscore, sortorder);

	double	*sortscore2 = NewArray <double> (itemNum);
	int	j, j1, add;
	for(i = 0; i < itemNum;)	{
		add = 0;
		for(j = i; j < itemNum; j ++)	{
			if(fabs(sortscore[j] - sortscore[i]) > 1e-16)	{
				break;
			}
			j1 = sortorder[j];
			if(optLen[j1] != 0)	sortscore2[add ++] = optRms[j1]/optLen[j1];
			else	sortscore2[add ++] = 0;
		}
		if(add > 1)	{
			qksort <double> (add, sortscore2, sortorder + i);
		}
		i = j;
	}
	delete[] sortscore2;
}

//sort alignment by probability, or whatever value given
void FATCATDB::SortAlign(void)
{
	if(sortorder != NULL)	return; 
	sortorder = NewArray <int> (itemNum);
	sortscore = NewArray <double> (itemNum);

	int	i;
	for(i = 0; i < itemNum; i ++)	sortorder[i] = i;

	if(itemvalue == NULL)	{
		printf("GetTarget must have been excuted\n");
		exit(1);
	}

	//GetTarget("fprobability", probability);

	for(i = 0; i < itemNum; i ++)	sortscore[i] = itemvalue[i];

	qksort <double> (itemNum, sortscore, sortorder);
}

//output the alignment result of the database-searching result (sorted by probability) 
void FATCATDB::WriteAlign(char *outfile)
{
	ofstream out(outfile);
	if(!out)	{
		printf("open file %s error\n", outfile);
		exit(1);
	}

	//SortAlign();
	SortByProbAdOther();
	int	i, n;
	int	batch = 500;
	if(!batch)	{
		for(i = 0; i < itemNum; i ++)	{
			n = sortorder[i];
			WriteAlign(code1[n], code2[n], out);	
			out<<"\n";
		}
	}
	else	{
		int 	num = 0;
		char	**name1 = NewMatrix <char> (batch, 100);
		char	**name2 = NewMatrix <char> (batch, 100);
		int	num2 = 0;
		printf("batch mod for item %d\n", itemNum);
		for(i = 0; i < itemNum; i ++)	{
			n = sortorder[i];
			strcpy(name1[num], code1[n]);
			strcpy(name2[num], code2[n]);
			num ++;
			if(num >= batch)	{
				WriteAlign(num, name1, name2, out);
				num2 += num;
				num = 0;
			}
		}
		printf("total read %d\n", num2);
		if(num > 0)	{
			WriteAlign(num, name1, name2, out);
			num2 += num;
		}
		printf("total read %d\n", num2);
		DelMatrix <char> (name1, batch);
		DelMatrix <char> (name2, batch);
	}
	out.close();
}

//------------------------------------------------------------------
void FATCATDB::WriteAlign(char *name1, char *name2, ofstream &out)
{
	char	str[1000];
	ifstream in(inputFile);
	if(!in)	{
		printf("open file %s error\n", inputFile);
		exit(1);
	}
	int	lab = 0;
	char	code1[200], code2[200];
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			sscanf(str, "%*s%s%*s%*s%s", code1, code2);
			if(!strcmp(code1, name1) && !strcmp(code2, name2))	{
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

//------------------------------------------------------------------
void FATCATDB::WriteAlign(int num, char **name1, char **name2, ofstream &out)
{
	if(num <= 0)	return;
	char	str[1000];
	char	**strrecord = NewMatrix <char> (num, 10000);
	ifstream in(inputFile);
	if(!in)	{
		printf("open file %s error\n", inputFile);
		exit(1);
	}
	int	i;
	int	lab = -1;
	char	code1[200], code2[200];
	int	readnum = 0;
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			sscanf(str, "%*s%s%*s%*s%s", code1, code2);
			lab = -1;
			for(i = 0; i < num; i ++)	{
				if(!strcmp(code1, name1[i]) && !strcmp(code2, name2[i]))	{
					lab = i;
					sprintf(strrecord[i], "%s\n", str);
					readnum ++;
					break;
				}
			}
		}
		else if(lab >= 0)	{
			sprintf(strrecord[lab], "%s%s\n", strrecord[lab], str); 
			if(!strncmp(str, "Note", 4) && readnum >= num)	break;
		}
	}
	in.close();
	if(readnum != num)	{
		printf("writealn error: %d vs %d\n", readnum, num);
	}

	for(i = 0; i < num; i ++)	{
		out<<strrecord[i]<<endl;
	}

	DelMatrix <char> (strrecord, num);
}

//------------------------------------------------------------------

//output a report of the database-searching result in one line per pair format
void FATCATDB::WriteReport(double cut, char *outfile)
{
	SortAlign();

	int	i, n;
	for(i = 0; i < itemNum; i ++)	{
		if(sortscore[i] > cut)	break;
	}
	int	showNum = itemNum;

	ofstream out(outfile);
	if(!out)	{
		printf("open outfile %s error\n", outfile);
		exit(1);
	}
	char	str[1000];
	out<<"#code1 code2 len1 len2 twist iniLen iniRms optLen optRms chainRms score alnLen gap probability\n";
	for(i = 0; i < showNum; i ++)	{
		n = sortorder[i];
		sprintf(str, "%-20s %-20s %4d %4d %2d %4d %5.2f %4d %5.2f    %5.2f %8.2f %4d %4.0f %.3e\n",
			code1[n], code2[n], len1[n], len2[n], twi[n], iniLen[n], iniRms[n], optLen[n], optRms[n],
			chainRms[n], score[n], alnLen[n], gap[n], itemvalue[n]);
		out<<str;
	}
	out.close();
}

//-----------------------------------------------------------------------------------------
void FATCATDB::WriteHtml(double cut, char *cgi, char *dir, char *outfile)
{
	SortAlign();

	ofstream out(outfile);
	if(!out)	{
		printf("open outfile %s error\n", outfile);
		exit(1);
	}
	char	str[1000], tmp[100], link[100];
	int	i, n;
	for(i = 0; i < itemNum; i ++)	{
		//printf("score %d %.2f\n", i, sortscore[i]);
		if(sortscore[i] > cut)	break;
	}
	int	showNum = i - 1;

	char	file[100];
	strcpy(file, inputFile);
	for(i = strlen(inputFile) - 1; i >= 0; i --)	{
		if(inputFile[i] == '/')	break;
	}	
	strcpy(file, inputFile + i + 1);
	sprintf(str, "FATCAT database searching report");
	WriteHtmlHead(str, out);
	CenterCaption(1, str, out);
	if(showNum <= 0)	{
		sprintf(str, "<p><center>No similar structures are found with probability < %.2e</center></p>\n", cut);
		out<<str;
	}
	else	{
		sprintf(str, "<p><center>The list of %d similar structures with probability < %.2e</center></p>\n", showNum, cut);
		out<<str;
		BegTable(out);
		sprintf(str, "<tr><td>pro1</td><td>pro2</td><td>len1</td><td>len2</td><td>score</td><td>P-value</td><td>twist</td>\n");
		out<<str;
		sprintf(str, "</td><td>opt-len</td><td>opt-rmsd</td><td>chain-rmsd</td><td>align-len</td><td>gap</td><td>alignment</td></tr>\n");
		out<<str;
		for(i = 0; i < showNum; i ++)	{
			n = sortorder[i];
			strcpy(tmp, code2[n]);
			if(!strncmp(tmp + strlen(tmp) - 4, ".pdb", 4))	{
				tmp[strlen(tmp) - 4] = '\0';
				sprintf(link, "http://scop.mrc-lmb.cam.ac.uk/scop/search.cgi?key=%s", tmp);
			}
			else	link[0] = 0;
			sprintf(str, "<tr><td>%s</td><td><a href=\"%s\">%s</td><td>%d</td><td>%d</td><td>%.2f</td><td>%.2e</td><td>%d</td><td>%d</td><td>%.2f</td><td>%.2f</td><td>%d</td><td>%.0f</td>\n",
				code1[n], link, code2[n], len1[n], len2[n], score[n], probability[n], twi[n], optLen[n], optRms[n],
				chainRms[n], alnLen[n], gap[n]);
			out<<str;
			sprintf(str, "<td><a href=\"%s\?-file=%s&-func=align&-n1=%s&-n2=%s\">view</a></td></tr>\n", cgi, dir, code1[n], code2[n]); 
			out<<str;
		}
		EndTable(out);
	}
	out<<"<p><p><center>\n";
	sprintf(str, "<form enctype=\"multipart/form-data\" action = \"%s\" method = \"get\">\n", cgi);
	out<<str;
	out<<"<input type=\"hidden\" checked name=\"-func\" value=\"list\">"<<endl; 
	sprintf(str, "<input type=\"hidden\" checked name=\"-file\" value=\"%s\">", dir);
	out<<str;
	out<<"<p><font color=red>Too</font> many or too few hits? Try a different probability threshold (0-0.2) <input type=\"text\" name=\"-prob\" value=\"0.1\" size=\"5\">"<<endl;
	out<<"<input type=\"submit\" name=\"send\" value=\"submit\"></p>"<<endl;
	out<<"<p> (The smaller the value is, the less hits you will get)</p>"<<endl;
	out<<"</form>"<<endl;

	WriteHtmlEnd(out);

	out.close();
}

//output a report of the database-searching result in one line per pair format
void FATCATDB::ExtractAlign(char *pro1, char *pro2, char *input, char *outfile)
{
	char	name1[100], name2[100];
	ifstream in(input);
	if(!in)	{
		printf("open file %s error\n", input);
		exit(1);
	}	
	ofstream out(outfile);
	if(!out)	{
		printf("open file %s error\n", outfile);
		exit(1);
	}	
	int	lab = 0;
	char	str[1000];
	while(!in.eof())	{
		str[0] = 0;
		in.getline(str, 200);
		if(!strncmp(str, "Align", 5))	{
			sscanf(str, "%*s%s%*s%*s%s", name1, name2);
			if((!strcmp(name1, pro1) && !strcmp(name2, pro2)) ||
			   (!strcmp(name1, pro2) && !strcmp(name2, pro1)))	{
				out<<str<<endl;
				lab = 1;
			}
			else	lab = 0;
		}
		else if(lab && !strncmp(str, "Note", 4))	{
			out<<str<<endl;
			break;
		}
		else if(lab)	{
			out<<str<<endl;
		}
	}
	out.close();
}

void FATCATDB::WriteProb(char *outfile)
{
	ofstream out(outfile);
	if(!out)	{
		printf("open file %s error\n", outfile);
		exit(1);
	}

	int	i, j;
	char	str[100];
	for(i = 0; i < itemNum; i ++)	{
		j = sortorder[i];
		sprintf(str, "%-20s %-20s %4d %4d %e\n", code1[j], code2[j], len1[j], len2[j], sortscore[i]);
		out<<str;
	}
	out.close();
}

//get the distribute of a given target
void FATCATDB::Distribute(int step_in)
{
	step = step_in;
	scalelist = NewArray <double> (step);
	distrib = NewArray <double> (step);

	int	i;
	distribute <double>  (itemNum, itemvalue, step, scalelist, distrib);
	double	windowlen = scalelist[1] - scalelist[0];

	for(i = 0; i < step; i ++)	{
		distrib[i] /= double(itemNum); //percentage
		distrib[i] /= windowlen; //divided by window-length
	}
}

//--------------------------------------
void FATCATDB::ExtTarget(double *value)
{
	int	i;
	for(i = 0; i < itemNum; i ++)	value[i] = itemvalue[i];
}

//--------------------------------------
int FATCATDB::GetTarget(char *target)
{
	if(itemvalue != NULL)	return 1;
	itemvalue = NewArray <double> (itemNum0);

	int	i;
	SIGEVA	sig;
	for(i = 0; i < itemNum; i ++)	{
		if(!strcmp(target, "oprobability"))	{
			itemvalue[i] = sig.calSigOldFlexi(len1[i], len2[i], score[i], twi[i]);
		}
		else if(!strcmp(target, "rprobability"))	{
			itemvalue[i] = sig.calSigRigid(len1[i], len2[i], score[i], optRms[i], optLen[i]);
		}
		else if(!strcmp(target, "fprobability"))	{
			itemvalue[i] = sig.calSigFlexi(len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "rprobabilitys2"))	{
			itemvalue[i] = sig.calSigAll(0, 2, len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "fprobabilitys2"))	{
			itemvalue[i] = sig.calSigAll(5, 2, len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "rprobabilitys1"))	{
			itemvalue[i] = sig.calSigAll(0, 1, len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "fprobabilitys1"))	{
			itemvalue[i] = sig.calSigAll(5, 1, len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "rprobabilitys3"))	{
			itemvalue[i] = sig.calSigAll(0, 3, len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "fprobabilitys3"))	{
			itemvalue[i] = sig.calSigAll(5, 3, len1[i], len2[i], score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "score"))		itemvalue[i] = double(score[i]);	
		else if(!strcmp(target, "alignlen"))	itemvalue[i] = double(optLen[i]);
		else if(!strcmp(target, "gaps"))	itemvalue[i] = double(gap[i]);		
		else if(!strcmp(target, "twists"))	itemvalue[i] = double(twi[i]);	
		else if(!strcmp(target, "rmsd"))	itemvalue[i] = double(optRms[i]);
		else if(!strcmp(target, "normscore"))	{
			itemvalue[i] = sig.normScore(score[i], optRms[i], optLen[i], twi[i]);
		}
		else if(!strcmp(target, "z-score"))	{
		}
		else	{
			printf("wrong target  %s\n", target);
			exit(1);
		}
		//printf("item %d, value %f\n", i, itemvalue[i]);
		//getchar();
	}
	return itemNum;
}

//write the distribute to a file
void FATCATDB::WriteDistribute(char *target, char *filename)
{
	ofstream out(filename);
	if(!out)	{
		printf("open file %s error\n", filename);
		exit(1);
	}
	out<<"# distribution of "<<target<<endl;
	int	i;
	char	info[100];
	for(i = 0; i < step; i ++)	{
		sprintf(info, "%-8.4f %-8.6f\n", scalelist[i], distrib[i]);
		out<<info;
	}
	out.close();
}

//draw a graph of distribute in postscript format 
void FATCATDB::PsDistribute(char *target, char *filename)
{
	PsShow	*ps = new PsShow(filename);
	ps->PageArea(100, 400, 200, 400);
	ps->setMarkFontSize(12);
	ps->setTextFontSize(12);

	double	**point = NewMatrix <double> (step, 2);
	int	ptnum = 0;
	int	i;
	for(i = 0; i < step; i ++)	{
		if(fabs(distrib[i]) < 1.0e-8)	continue;	
		point[ptnum][0] = scalelist[i];
		point[ptnum][1] = distrib[i];
		ptnum ++;
	}
	double	color[] = {0, 0, 0};

	int     ifcodtransform = 1;
	ps -> drawAxis(point, ptnum, target, "frequency");
	//ps -> setTransform(); //set the scale & translate
	ps -> prtPoint(point, ptnum, ifcodtransform, "cross", color);

	delete ps;
	DelMatrix <double> (point, step);
}

int FATCATDB::GetItemNum(void)
{
	return itemNum;
}
