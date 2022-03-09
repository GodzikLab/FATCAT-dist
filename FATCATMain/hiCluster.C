//---------------------------------------------------------------------------------
//A class for clustering
//single-linkage clustering, average-linkage clustering are realized
//corresponding class: hiCluster.h
//by Y.Y, 10/22/02 latest update
//---------------------------------------------------------------------------------

#include "hiCluster.h"
#include "clustdisTemp.h"

using namespace STATIS;

hiCluster::hiCluster(void)
{
	initial();
}

//assign the matrix from input 
hiCluster::hiCluster(int len, double **dist)
{
	initial();
	assignCluster(len, dist);
}

//assign the matrix from input 
hiCluster::hiCluster(int len, char **name, double **dist)
{
	initial();
	assignCluster(len, name, dist);
}

void hiCluster::initial(void)
{
	itemn = 0;
	items = NULL;
	matrix = NULL;
	matrix2 = NULL;
	clustnum = 0;
	clustmemb = NULL;
	treestr = NULL;
	clsize = NULL;
	cllist = NULL;
	memb2clust = NULL;
	clustsize = NULL;
	clustscore = NULL;
	clustscore2 = NULL;
}

//----------------------------------------------
hiCluster::~hiCluster(void)
{
	int	i;
	for(i = 0; i < itemn; i ++)	{
		delete[] matrix[i];
		delete[] items[i];
		delete[] memb2clust[i];
		delete[] clustsize[i];
		delete[] clustscore[i];
		if(matrix2 != NULL)	delete[] matrix2[i];
		if(clustscore2 != NULL)	delete[] clustscore2[i];
	} 
	delete[] matrix;
	if(matrix2 != NULL)	delete[] matrix2;
	if(clustscore2 != NULL)	delete[] clustscore2;
	delete[] items;
	delete[] memb2clust;
	delete[] clustsize;
	delete[] clustscore;
	if(treestr != NULL)	delete[] treestr;

	if(cllist != NULL)	{
		for(i = 0; i < clnum; i ++)	delete[] cllist[i];
		delete[] cllist;
		delete[] clsize;
	}
	delClust(clustmemb[0]);
	delete[] clustmemb;
}

//------------------------------------------------
void hiCluster::delClust(clust *oneclust)
{
	if(oneclust->left != NULL)	delClust(oneclust->left);
	if(oneclust->right != NULL)	delClust(oneclust->right);
	delete oneclust;
}

//average-linkage-clustering
void hiCluster::averageLinkageClust()
{
	printf("into 1, itemn %d\n", itemn);
	runclust(averageLinkage);
}

//complete-linkage-clustering
void hiCluster::completeLinkageClust()
{
	runclust(completeLinkage);
}

void hiCluster::assignCluster(int len, double **dist)
{
	int	i;
	char	**name = new char*[len]; 
	for(i = 0; i < len; i ++)	{
		name[i] = new char[50];
		sprintf(name[i], "%d", i + 1);
	}
	assignCluster(len, name, dist);
	for(i = 0; i < len; i ++)	delete[] name[i];
	delete[] name;
}
	
void hiCluster::assignCluster(int len, char **name, double **dist)
{
	initial();
	initial(len);

	int	i, j;
	maxstrlen = 0;
	for(i = 0; i < itemn; i ++)	{
		strcpy(items[i], name[i]);
		maxstrlen += strlen(name[i]);
		for(j = 0; j < itemn; j ++)	{
			matrix[i][j] = dist[i][j];
		}
	}

	maxstrlen += 4 * itemn;

	printf("in assignCluster: item %d\n", itemn);
}

//calloc the variants
void hiCluster::initial(int len)
{
	itemn = len;
	matrix = new double*[itemn];
	items = new char*[itemn];
	int	i, j;
	for(i = 0; i < itemn; i ++)	{
		matrix[i] = new double[itemn];
		items[i] = new char[50];
	}	
	for(i = 0; i < itemn; i ++)	{
		for(j = 0; j < itemn; j ++)	matrix[i][j] = 0;
		items[i][0] = 0;
	}
	memb2clust = new int*[itemn];
	clustsize = new int*[itemn];
	clustscore = new double*[itemn];
	for(i = 0; i < itemn; i ++)	{
		memb2clust[i] = new int[itemn];
		clustsize[i] = new int[itemn];
		clustscore[i] = new double[itemn];
		for(j = 0; j < itemn; j ++)	{
			memb2clust[i][j] = 0;
			clustsize[i][j] = 0;
			clustscore[i][j] = 0;
		}
	}
}

//------------------------------------------------
void hiCluster::assignSecondDis(double **dis2)
{
	matrix2 = new double*[itemn];
	clustscore2 = new double*[itemn];
	int	i, j;
	for(i = 0; i < itemn; i ++)	{
		matrix2[i] = new double[itemn];
		clustscore2[i] = new double[itemn];
		for(j = 0; j < itemn; j ++)	{
			matrix2[i][j] = clustscore2[i][j] = dis2[i][j];
		}
	}	
}


//single-linkage-clustering
void hiCluster::singleLinkageClust()
{
	runclust(singleLinkage);
}

void hiCluster::runclust(char *method)
{
	if(!strncmp(method, "average", 7))	runclust(averageLinkage); 
	else if(!strncmp(method, "complete", 8))runclust(completeLinkage); 
	else if(!strncmp(method, "single", 6))	runclust(singleLinkage); 
	else	{
		printf("undefined clustering method\n");
		exit(1);
	}
}

//*funcs: the function of calculate cluster-distance
//could be average-linkage, single-linkage, complete-linkage
void hiCluster::runclust(double (funcs)(int len, double *dis))
{
	initClust();
	printf("initclust..itemn %d, clustnum %d\n", itemn, clustnum);

	printf("now clust..\n");
	int	i, j, mini, minj;
	mini = minj = 0;
	double	mindis, dis, mindis2;
	mindis = dis = mindis2 = 0;

	while(clustnum > 1)	{
		mindis = 10000;
		for(i = 0; i < clustnum - 1; i ++)	{
			for(j = i + 1; j < clustnum; j ++)	{
				dis = clustscore[i][j];
				if(mindis > dis)	{
					mindis = dis;
					if(matrix2 != NULL)	mindis2 = clustscore2[i][j];
					mini = i;
					minj = j;
				}	
				else if(fabs(mindis - dis) < 1e-10 && matrix2 != NULL && mindis2 > clustscore2[i][j])	{
					mindis = dis;
					mindis2 = clustscore2[i][j];
				       	mini = i;
					minj = j;	
				} //first distance equals, determined by the secondary distance value; special case FATCAT2Tree
			}
		}
		combineClust(mini, minj, mindis, funcs);
	}
}

//depth-first searching
void hiCluster::branchSearch(clust *root, int *n, int *list)
{
	if(root->left != NULL)	{
		branchSearch(root->left, n, list);
		branchSearch(root->right, n, list);
	}
	else	list[(*n) ++] = root->index;
}

//initial the clustering
void hiCluster::initClust(void)
{
	clustnum = itemn;
	printf("In initClust: itemn %d, clustnum %d\n", itemn, clustnum);
	clustmemb = new clust*[itemn];
	int	i, j;
	for(i = 0; i < itemn; i ++)	{
		clustmemb[i] = new clust[itemn];
		clustmemb[i]->index = i;
		clustmemb[i]->size = 1;
		clustmemb[i]->dis = 0;
		clustmemb[i]->left = NULL;
		clustmemb[i]->right = NULL;
		memb2clust[clustnum - 1][i] = i; //note: clustnum -1 instead of clustnum
		clustsize[clustnum - 1][i] = 1;
	}
	//each item forms a seperate a cluster initially
	for(i = 0; i < itemn; i ++)	{
		for(j = i + 1; j < itemn; j ++)	{
			clustscore[i][j] = clustscore[j][i] = matrix[i][j];
			if(matrix2 != NULL)	clustscore2[i][j] = clustscore2[j][i] = matrix2[i][j];
		}
	}
}

//combine two clusters
void hiCluster::combineClust(int c1, int c2, double dis, double (funcs)(int len, double *dlist))
{
	int	cmin = (c1 < c2)?c1:c2;
	int	cmax = (c1 > c2)?c1:c2;

	clust*  newclust = new clust;
	newclust->index = 0;
	newclust->dis = dis;
	newclust->size = clustmemb[c1]->size + clustmemb[c2]->size;
	newclust->left = clustmemb[c1];
	newclust->right = clustmemb[c2];
	//merge two clusters
	
	clustmemb[cmin] = newclust;
	//replace the cmin-cluster with the new cluster

	if(cmax != clustnum - 1)	{
		for(int i = cmax; i < clustnum - 1; i ++)	{
			clustmemb[i] = clustmemb[i + 1];
		}
	} //resort the clust list, remove the cmax

	int	i, j, n, m;
	for(i = 0; i < clustnum - 1; i ++)	{
		for(j = 0; j < clustnum - 1; j ++) {
			if(i >= cmax && j >= cmax)	{
				clustscore[i][j] = clustscore[i + 1][j + 1];		
				if(matrix2 != NULL)	clustscore2[i][j] = clustscore2[i + 1][j + 1];
			}
			else if(i >= cmax)		{
				clustscore[i][j] = clustscore[i + 1][j];
				if(matrix2 != NULL)	clustscore2[i][j] = clustscore2[i + 1][j];
			}
			else if(j >= cmax)		{
				clustscore[i][j] = clustscore[i][j + 1];	
				if(matrix2 != NULL)	clustscore2[i][j] = clustscore2[i][j + 1];	
			}
		}
	} //resort the clust-score

	int	*d1 = new int[itemn];
	int	*d2 = new int[itemn];
	double	*dlist = new double[itemn * itemn];
	int	add, bn, bm;
	double	d;
	bn = 0;
	branchSearch(clustmemb[cmin], &bn, d1);
	for(i = 0; i < clustnum - 1; i ++)	{
		if(i == cmin)	continue;
		bm = 0;
		branchSearch(clustmemb[i], &bm, d2);
		add = 0;
		for(n = 0; n < bn; n ++)	{
			for(m = 0; m < bm; m ++)	{
				dlist[add ++] = matrix[d1[n]][d2[m]];
			}
		}
		d = funcs(add, dlist);
		clustscore[i][cmin] = clustscore[cmin][i] = d;
		if(matrix2 != NULL)	{
			add = 0;
			for(n = 0; n < bn; n ++)	{
				for(m = 0; m < bm; m ++)	{
					dlist[add ++] = matrix2[d1[n]][d2[m]];
				}
			}
			d = funcs(add, dlist);
			clustscore2[i][cmin] = clustscore2[cmin][i] = d;
		}
	} //update the clust-score between the new cluster with the old clusters
	delete[] d1;
	delete[] d2;
	delete[] dlist;

	clustnum --;
	//clustnum-th memb2clust and clustsize recorded in memb2clust[clustnum - 1] and clustsize[clustnum - 1]
	//so excute clustnum -- first then update the memb2clust and clustsize

	for(i = 0; i < itemn; i ++)	{
		if(memb2clust[clustnum][i] == cmax)	memb2clust[clustnum - 1][i] = cmin;
		else if(memb2clust[clustnum][i] > cmax)	memb2clust[clustnum - 1][i] = memb2clust[clustnum][i] - 1;
		else					memb2clust[clustnum - 1][i] = memb2clust[clustnum][i];
	}
	for(i = 0; i < clustnum - 1; i ++)	{
		if(i == cmin)	clustsize[clustnum - 1][i] += clustsize[clustnum][cmax]; 
		else if(i >= cmax)	{
			clustsize[clustnum - 1][i] = clustsize[clustnum - 1][i + 1];
		}
	} //update the member-cluster information
}

//-------------------------------------------------------
int hiCluster::extclust(double discut)
{
	int	i, j, k;
	int	*clsize_t = new int[itemn];
	int	**cllist_t = new int*[itemn];
	for(i = 0; i < itemn; i ++)	{
		clsize_t[i] = 0;
		cllist_t[i] = new int[itemn];
	}

	int	n = 0;
	extclust(clustmemb[0], discut, &n, clsize_t, cllist_t);
	//search from the root
	clnum = n;
	//extract the clusters

	int	*sort = new int[clnum];
	for(i = 0; i < clnum; i ++)	sort[i] = i;
       	qksort <int> (n, clsize_t, sort);
       	reverse <int> (n, clsize_t, sort); 
	//sort the clusterers by decreasing size

	if(clsize == NULL)	{
		clsize = new int[clnum];
		cllist = new int*[clnum];
		for(i = 0; i < clnum; i ++)	{
			clsize[i] = clsize_t[i];
			cllist[i] = new int[clsize[i]];
		}
	}
	for(i = 0; i < clnum; i ++)	{
		clsize[i] = clsize_t[i];
		j = sort[i];
		for(k = 0; k < clsize[i]; k ++)	{
			cllist[i][k] = cllist_t[j][k];
		}
	}
	//assign the cluster information to global variants

	printf("discut %f, clustnum %d\n", discut, clnum);
	for(i = 0; i < clnum; i ++)	{
		printf("clust %d, member %d\n", i, clsize[i]);
	}

	delete[] clsize_t;
	for(i = 0; i < itemn; i ++)	{
		delete[] cllist_t[i];
	}
	delete[] cllist_t;
	delete[] sort;
	return clnum;
}

//-------------------------------------------------------
void hiCluster::extclust(clust* root, double discut, int *clnum, int *clsize, int **cllist)
{
	//printf("discut %f, root dis %f\n", discut, root->dis);
	printf("MARK!! this root dis %.3e\n", root->dis);
	if(root->dis <= discut || root->left == NULL)	{
		//stop if distance between clusters is larger than discut, or confront a leave
		if(root->left == NULL)	{
			cllist[*clnum][0] = root->index;
			clsize[(*clnum) ++] = 1;	
		}
		else	{
			branchSearch(root, &clsize[*clnum], cllist[*clnum]);
			//printf("cluster %d, size %d\n", *clnum, clsize[*clnum]);
			(*clnum) ++;
		} //depth-first searching
	} 
	else	{
		extclust(root->left, discut, clnum, clsize, cllist);
		extclust(root->right, discut, clnum, clsize, cllist);
	}
}

void hiCluster::extclust(int **mb2clust)
{
	int	i, j;
	for(i = 0; i < itemn; i ++)	{
		for(j = 0; j < itemn; j ++)	{
			mb2clust[i][j] = memb2clust[i][j];
		}
	} 
}

//derive cluster structure given number of clusters
//represented in {0,1} matrix, eg. i, j belong to a cluster, Si,j = 1, else = 0
//-------------------------------------------------------
int hiCluster::extclust(int clnum, int voidsigleton, int **nbmatrix)
{
	int	*size = new int[itemn];
	int	**list = new int*[itemn];
	int	i, j, k, c, nonsigleton;
	int	success = 1;
	for(i = 0; i < itemn; i ++)	{
		list[i] = new int[itemn];
	}
	if(voidsigleton)	{	
		for(i = clnum; i < itemn; i ++)	{
			getclust(i, size, list);
			nonsigleton = 0;
			for(k = 0; k < i; k ++)	{
				printf("search %d: cluster %d, size %d\n", i, k, size[k]);
				if(size[k] > 1)	nonsigleton ++;
			}
			if(nonsigleton >= clnum)	break;
		}
		c = i;
	}
	else	{
		getclust(clnum, size, list);
		c = clnum;
	}

	if(c < itemn)	{
		for(i = 0; i < itemn; i ++)	{
			for(j = 0; j < itemn; j ++)	{
				if(i == j)	nbmatrix[i][j] = 1;
				else	nbmatrix[i][j] = 0;
			}	
		}
		for(i = 0; i < c; i ++)	{
			for(j = 0; j < size[i]; j ++)	{
				for(k = j + 1; k < size[i]; k ++)	{
					nbmatrix[list[i][j]][list[i][k]] = nbmatrix[list[i][k]][list[i][j]] = 1;
				}
			}
		}
	}
	else	success = 0;
	delete[] size;
	for(i = 0; i < itemn; i ++)	delete[] list[i];
	delete[] list;

	return success;
}

//-------------------------------------------------------
void hiCluster::getclust(int clnum, int *csize, int **clist)
{
	int	i, c;
	for(i = 0; i < clnum; i ++)	csize[i] = 0;
	for(i = 0; i < itemn; i ++)	{
		c = memb2clust[clnum - 1][i]; //note: clnum - 1 instead of clnum
		clist[c][csize[c] ++] = i; 
	}
}

//------------------------------------------------------------
void hiCluster::exttree(clust *root, char *info)
{
	if(root->left != NULL)	{
		char  *info1 = new char[maxstrlen]; 
		char  *info2 = new char[maxstrlen]; 
		exttree(root->left, info1);
		exttree(root->right, info2);
		sprintf(info, "(%s,%s)", info1, info2);
		delete[] info1;
		delete[] info2;
	}
	else	{
		sprintf(info, "%s", items[root->index]);
	}	
}

//------------------------------------------------------------
void hiCluster::writetree(char *outfile)
{
	ofstream out(outfile);
	if(!out)	{
		printf("open tree outfile error\n");
		exit(1);
	}
	if(treestr == NULL)	{
		char	*info = new char[maxstrlen];
		exttree(clustmemb[0], info);
		out<<info;
		delete[] info;
	}
	else	{
		out<<treestr;
	}
	out.close();
}

//------------------------------------------------------------
char* hiCluster::exttree(void)
{
	char	*info;
	if(treestr == NULL)	{
		info = new char[maxstrlen];
		exttree(clustmemb[0], info);
	}
	else	{
		info = new char[strlen(treestr) + 1];
		strcpy(info, treestr);
	}
	return info;
}

//------------------------------------------------------------
void hiCluster::writeclust(int sizecut, ofstream &stream)
{
	char	info[1000];
	sprintf(info, "#total %d clusters\n", clnum);
	stream<<info;
	int	i, j, k;
	for(i = 0; i < clnum; i ++)	{
		if(clsize[i] < sizecut)	continue;
		sprintf(info, "#%d cluster\n", i + 1);
		stream<<info;
		for(j = 0; j < clsize[i]; j ++)	{
			k = cllist[i][j];
			stream<<items[k]<<endl;
		}
		stream<<endl;
	}
}

clust* hiCluster::getroot(void)
{
	return clustmemb[0];
}

//Creating subtree string from its parent full tree that only includes the molecules of the given list
//-----------------------------------------------------------------------------
char* hiCluster::subtreestr(int mnum, int *mlist)
{
	char	*sub = new char[maxstrlen];
	sub[0] = 0;
	int	*musage = new int[itemn];
	int	i;
	for(i = 0; i < itemn; i ++)	musage[i] = -1;
	for(i = 0; i < mnum; i ++)	musage[mlist[i]] = 1;

	for(i = 0; i < itemn; i ++)	{
		printf("musage: %d %d\n", i, musage[i]);
	}

	writesubtreestr(clustmemb[0], musage, sub);

	return sub;
}

//This function is used by SubTreeStr()
void hiCluster::writesubtreestr(clust *curr, int *musage, char *str)
{
	if(curr -> left == NULL)	{
		int	i = curr -> index;
		if(musage[i] != -1)	strcpy(str, items[i]);
		else	str[0] = 0;
	}
	else	{
		char	*str1 = new char[10000];
		char	*str2 = new char[10000];
		str1[0] = str2[0] = 0;
		writesubtreestr(curr -> left, musage, str1);
		writesubtreestr(curr -> right, musage, str2);
		if(strlen(str1) > 0 && strlen(str2) > 0)	{
			sprintf(str, "(%s,%s)", str1, str2); 
		}
		else if(strlen(str1) > 0)	{
			strcpy(str, str1);
		}
		else if(strlen(str2) > 0)	{
			strcpy(str, str2);
		}
		delete[] str1;
		delete[] str2;
	}
}
