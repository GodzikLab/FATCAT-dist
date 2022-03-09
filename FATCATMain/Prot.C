//--------------------------------------------------------------
//simple Protein class
//read CA coordinations and solvent accessible
//--------------------------------------------------------------
#include <cstring>
#include <cstdlib>
#include "Prot.h"
#include "Amino.h"
using namespace AMINO;
//--------------------------------------------------------------
PROT::PROT(char *file)
{
	Initial();
	AssignPdb(file, 0, 100000);
	SAA();
}

//--------------------------------------------------------------
PROT::PROT(char *file, int begres, int reslen)
{
	Initial();
	AssignPdb(file, begres, reslen);
	SAA();
}

//--------------------------------------------------------------
void PROT::SAA(void)
{
	saa = NewArray <char> (length + 1);
	int	i;
	for(i = 0; i < length; i ++)	{
		saa[i] = aa1Get(res[i]);
	}
	saa[i] = '\0';
}

//--------------------------------------------------------------
void PROT::GetSAA(char *saa_o)
{
	strcpy(saa_o, saa); 
}

//--------------------------------------------------------------
void PROT::GetIndex(int *index_o)
{
	int	i;
	for(i = 0; i < length; i ++)	{
		index_o[i] = 0;
		sscanf(index[i], "%d", &index_o[i]);
	}
}

//--------------------------------------------------------------
int  PROT::GetIndex(int index_i)
{
	int	i;
	sscanf(index[index_i], "%d", &i);
	return i;
}

//--------------------------------------------------------------
int PROT::ResIndex(char *index_i)
{
	int	i;
	for(i = 0; i < length; i ++)	{
		if(!strcmp(index_i, index[i]))	break; 
	}
	if(i == length)	return -1;
	else	return i;
}

//--------------------------------------------------------------
PROT::PROT(void)
{
	Initial();
}

//--------------------------------------------------------------
PROT::~PROT(void)
{
	if(res != NULL)		DelMatrix <char> (res, length0);
	if(chain != NULL)	DelArray <char> (chain);
	if(index != NULL)	DelMatrix <char> (index, length0);
	if(caCod != NULL)	DelArray <double> (caCod);	
	if(asa != NULL)		DelArray <double> (asa);	
	if(res2Atm != NULL)	DelArray <int> (res2Atm);
	if(atmCod != NULL)	DelArray <double> (atmCod);
	if(atmName != NULL)	DelMatrix <char> (atmName, atm0);
	if(atm2Res != NULL)	DelArray <int> (atm2Res);
	if(saa != NULL)		DelArray <char> (saa);
}

//--------------------------------------------------------------
void PROT::Initial(void)
{
	length = 0;
	res = NULL;
	chain = NULL;
	index = NULL;
	caCod = NULL;
	asa = NULL;
	res2Atm = NULL;

	totalAtm = 0;
	atmCod = NULL;
	atmName = NULL;
	atm2Res = NULL;
	saa = NULL;
}

//--------------------------------------------------------------
void PROT::Calloc(int len)
{
	if(caCod != NULL)	return; //already calloced
	length = length0 = len;
	res = NewMatrix <char> (len, 5);
	chain = NewArray <char> (len + 1);
	index = NewMatrix <char> (len, 5);
	caCod = NewArray <double> (3 * len);
	asa = NewArray <double> (len);
	res2Atm = NewArray <int> (len + 1);
}

//--------------------------------------------------------------
void PROT::CallocAtm(int atm)
{
	if(atmCod != NULL)	return;
	totalAtm = atm0 = atm;
	atmCod = NewArray <double> (3 * atm);
	atmName = NewMatrix <char> (atm, 6);
	atm2Res = NewArray <int> (atm);
}

//read processed-pdb files with only Ca coordination
//--------------------------------------------------------------
void PROT::AssignPdbCa(char *file)
{
	int	len1, len2;
	len1 = len2 = 0;
	char	str[50000];
	ifstream pdb(file);
	if(!pdb) {  printf("open pdb file error in AssignPdbCa\n"); exit(1); }
	while(!pdb.eof())	{
		pdb.getline(str, 5000);
		if(!strncmp(str, "#Length", 7))	{
			sscanf(str, "%*s%d", &len1);
			Calloc(len1);
		}
		if(str[0] == '#' || str[0] == '\0')	continue;
		if(caCod == NULL) { printf("caCod is not calloced in AssignPdbCa\n"); exit(1); }	
		if(len2 >= length) { printf("wrong in residue number in AssignPdbCa\n"); exit(1); }	
		sscanf(str, "%s", index[len2]);
		sscanf(str + 8, "%s%lf%lf%lf\n", res[len2], 
			&caCod[3 * len2], &caCod[3 * len2 + 1], &caCod[3 * len2 + 2]);
		len2 ++;
	}
	pdb.close();
}

//read original-pdb files with all atoms 
//--------------------------------------------------------------
void PROT::AssignPdb(char *file, int begres, int reslen)
{
//printf("file %s, begres %d, reslen %d\n", file, begres, reslen);

	int	i, len, atm, ca, k;
	char	str[50000], oldline[50000];
	ifstream pdb(file);
	if(!pdb) {
		printf("open pdb file %s error\n", file);
		exit(1);
	}

	//extract length of residues and atoms first
	len = atm = 0;
	while(!pdb.eof())	{
		str[0] = '\0';
		pdb.getline(str, 5000);
		if(!strncmp(str, "ENDMDL", 6))	break;
		if(!strncmp(str, "ATOM", 4))	{
			if(strncmp(str + 22, oldline + 22, 5) || len == 0)	len ++;	
			//bug reported by Martin for case d1tlca_
			atm ++;
		}	
		strcpy(oldline, str);
	}
	pdb.close();

	//calloc the varibles
	Calloc(len);
	CallocAtm(atm);

	//assign the varibles
	int	*caread = NewArray <int> (len);
	len = atm = ca = 0;
	int	add = -1;
	int	len_bf, atm_bf;
	len_bf = atm_bf = 0;
	ifstream pdb2(file);
	while(!pdb2.eof())	{
		str[0] = '\0';
		pdb2.getline(str, 5000);
		if(!strncmp(str, "ENDMDL", 6))	break;
		if(!strncmp(str, "ATOM", 4))	{
			if(add == -1 || strncmp(str + 22, oldline + 22, 5))	{
				add ++;
				if(add < begres)	{
					strcpy(oldline, str);
					continue;	
				}
				if(add >= reslen + begres)	break; 
					//note: only part structure when begres and reslen are used	
				if(len >= length0)	{
					printf("length0 %d\nthis line %s\n", length0, str);
					printf("wrong residue read\n");
					exit(1);
				}
				if(len >= 1 && caread[len - 1] == 0)	{
					//printf("Warning: pdb %s, the %s residue has no CA atom, removed!\n", file, index[len - 1]);
					len = len_bf;
					atm = atm_bf;
				} //the residues without CA atom is not considered, Y.Y, Dec 8, 03
				sscanf(str + 17, "%s", res[len]);
				sscanf(str + 22, "%s", index[len]);
				res2Atm[len] = atm;
				chain[len] = str[21];
				caread[len] = 0;
				len_bf = len;
				atm_bf = atm;
				len ++;	
			} //a new residue
			if(add < begres)	continue;	
			if(atm >= atm0)	{
				printf("wrong atm read\n");
				exit(1);
			}
			sscanf(str + 12, "%s", atmName[atm]);
			sscanf(str + 30, "%lf%lf%lf", &atmCod[3*atm], &atmCod[3*atm+1], &atmCod[3*atm+2]);
			if(!strcmp(atmName[atm], "CA"))	{
				if(caread[len - 1] == 0)	{
					for(k = 0; k < 3; k ++)	caCod[3*ca+k] = atmCod[3*atm+k];
					ca ++;
				} //only one CA is read for each residue, Y.Y, Dec,8,2003, annoyed by so many not-perfect pdb files
				caread[len - 1] ++;
			}
			atm2Res[atm] = len - 1;
			atm ++;
		}	
		strcpy(oldline, str);
	}
	if(len >= 1 && caread[len - 1] == 0)	{
		//printf("Warning: pdb %s, the %s residue has no CA atom, removed!\n", file, index[len - 1]);
		len = len_bf;
		atm = atm_bf;
	}
	chain[len] = '\0';
	pdb2.close();
	res2Atm[len] = atm;

	length = len;
	totalAtm = atm;


	if(ca != len)	{
		for(i = 0; i < len; i ++) {
			if(caread[i] == 0)	{
				printf("Warning: pdb %s, the %s residue has no CA atom!\n", file, index[i]); exit(1);
			}
			if(caread[i] > 1)	{
				printf("Warning: pdb %s, the %s residue has two CA atom!\n", file, index[i]); exit(1);
			}
		}
		//sprintf(str, "read pdb %s error: %d residues have no CA atoms, residue %d, ca %d!\n", file, len - ca, len, ca);
		//DelArray <int> (caread);
		//Y.Y, Dec,8,2003, annoyed by so many not-perfect pdb files
	}
	DelArray <int> (caread);
}

//--------------------------------------------------------------
void PROT::AssignPdb(PROT *pro)
{
	//calloc the varibles
	Calloc(pro->length);
	CallocAtm(pro->totalAtm);

	//assign coordinations and residues from input pdb
	int	i, k;
	strcpy(chain, pro->chain);
	for(i = 0; i < length; i ++)	{
		strcpy(res[i], pro->res[i]);
		strcpy(index[i], pro->index[i]);
		for(k = 0; k < 3; k ++)	{
			caCod[3 * i + k] = pro->caCod[3 * i + k];
		}
		asa[i] = pro->asa[i];
		res2Atm[i] = pro->res2Atm[i];
	}
	res2Atm[i] = pro->res2Atm[i];
	for(i = 0; i < totalAtm; i ++) {
		for(k = 0; k < 3; k ++)	{
			atmCod[3 * i + k] = pro->atmCod[3 * i + k];
		}
		strcpy(atmName[i], pro->atmName[i]);
		atm2Res[i] = pro->atm2Res[i];
	}
}

//output the pdb to the stream
//---------------------------------------------------------------------------
void PROT::PrintPdb(char forcechain, ofstream &stream)
{
	int	i, j;
	char	str[1000];
	for(i = 0; i < totalAtm; i ++)	{
		j = atm2Res[i];
	        sprintf(str, "ATOM%7d%2c%-4s%-4s%c%4s%4c%8.3f%8.3f%8.3f\n",
	                     i + 1, ' ', atmName[i], res[j], forcechain, index[j], ' ',
                             atmCod[3*i+0], atmCod[3*i+1], atmCod[3*i+2]);
		stream<<str;
	}
}

//return the coordinations for a given residue list
double* PROT::Cod4Res(int n, int *list)
{
	double *cod = NewArray <double> (n * 3);
	int	i, j, k;
	for(i = 0; i < n; i ++)	{
		j = list[i];
		for(k = 0; k < 3; k ++)	{
			cod[3 * i + k] = caCod[3 * j + k];
		}
	}
	return cod;
}

//return the rmsd of the CAs between the input pro and this protein, at given positions
double PROT::CalCaRmsd(PROT *pro, int resn, int *res1, int *res2)
{
	double	*cod1 = Cod4Res(resn, res1);
	double	*cod2 = pro->Cod4Res(resn, res2);
	double	r[9], t[3];
	double	rmsd = kearsay(resn, cod1, cod2, r, t);
	DelArray <double> (cod1);
	DelArray <double> (cod2);
	return rmsd;
}

//transformed pdb coordination according to input transformation
void PROT::TransCod(double *r, double *t)
{
	int	i;
	for(i = 0; i < length; i ++)	{
		tran_ord(r, t, caCod + 3 * i);
	}
	for(i = 0; i < totalAtm; i ++)	{
		tran_ord(r, t, atmCod + 3 * i);
	}
}

//transformed pdb coordination from r1 to r2, according to input transformation
//be cautious to use this function
//it will cause breaks in the protein
void PROT::TransCod(double *r, double *t, int r1, int r2)
{
	int	i;
	for(i = r1; i < r2; i ++)	{
		tran_ord(r, t, caCod + 3 * i);
	}
	//printf("r1 %d, atm %d, r2 %d, atm %d\n", r1, res2Atm[r1], r2, res2Atm[r2]);
	for(i = res2Atm[r1]; i < res2Atm[r2]; i ++)	{
		tran_ord(r, t, atmCod + 3 * i);
	}
}

int PROT::GetLength(void)
{
	return length;
}
 
double** PROT::DisTable(int maxlen)
{
	double	**dis = NewMatrix <double> (length, length);
	int	i, j, k;
	for(i = 0; i < length; i ++)	{
		dis[i][i] = 0;
		for(j = i + 1; j < length && j <= i + maxlen; j ++)	{
			dis[i][j] = 0;
			for(k = 0; k < 3; k ++)	{
				dis[i][j] += (caCod[3*i+k]-caCod[3*j+k]) * (caCod[3*i+k]-caCod[3*j+k]);
			}
			dis[i][j] = sqrt(dis[i][j]);
			dis[j][i] = dis[i][j];
		}
	}
	return dis;
}
