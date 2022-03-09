//---------------------------------------------------------------------------------
//A class for clustering
//single-linkage clustering, average-linkage clustering are realized
//by Y.Y, 10/22/02 latest update
//---------------------------------------------------------------------------------

#ifndef HiCLUSTER_H
#define HiCLUSTER_H

#include "basic.h"
#include "tempkit.h"

typedef struct clustdef
{
	int	index;
	int	size;
	double 	dis;	
	struct clustdef* left;
	struct clustdef* right;
} clust; 


class hiCluster	{
	private:
		int	itemn;	 //number of items
		char	**items; //labels of each item
		double	**matrix;//the distance matrix of items
		double	**matrix2; //the secondary distance matrix (functions when the first distance is tie)
		int	clustnum;    //the number of clusters in clustering process
		clust	**clustmemb; //the member list of each cluster
		int	maxstrlen;	
		double	**clustscore;	//the stored score for each cluster pair
		double	**clustscore2;	//the stored score for each cluster pair

		void	initial(void);
		void	initial(int len);
		void	initClust(void);
		void	combineClust(int c1, int c2, double dis, double (funcs)(int len, double *dis));
		void	delClust(clust *oneclust);
	public:
		int	**memb2clust; //record the clustering result in different steps of clustering
		int	**clustsize; //record the clustering result in different steps of clustering
		char	*treestr;    //the tree format of the clustering
		int	clnum;       //clnum, clsize, cllist used in specific cluster dissection
		int	*clsize;     //according to a distance threshould
		int	**cllist;
	public:
		hiCluster(void);
		hiCluster(int len, double **matrix);
		hiCluster(int len, char **name, double **matrix);
		~hiCluster(void);
		void	assignSecondDis(double **dis2);
		void	assignCluster(int len, double **matrix);
		void	assignCluster(int len, char **name, double **matrix);
		void	averageLinkageClust(void); //average-linkage-clustering
		void	completeLinkageClust(void);//compelete-linkage-clustering
		void	singleLinkageClust(void);  //single-linkage-clustering
		void	runclust(double (funcs)(int len, double *dis));
		void	runclust(char *method);
		void	branchSearch(clust *root, int *n, int *list);
		void	extclust(int **mb2clust);
		int	extclust(double discut);
		void 	extclust(clust* root, double discut, int *clnum, int *clsize, int **cllist);
		int	extclust(int clnum, int voidsigleton, int **nbmatrix);
		void	getclust(int clnum, int *csize, int **clist);
		void	exttree(clust *root, char *info);
		char*	exttree(void); //output the clustering result in tree format, to be viewed
		void	writetree(char *out);
		void	writeclust(int sizecut, ofstream &out);
		clust*	getroot(void);
		char*	subtreestr(int mnum, int *mlist);
		void	writesubtreestr(clust *curr, int *musage, char *str);
};

#endif
