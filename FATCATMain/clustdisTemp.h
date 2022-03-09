#ifndef CLUSTDIS_H
#define CLUSTDIS_H

#include "basic.h"

//The dissimilarity between clusters is calculated using average values. 
//It is UPGMA - Unweighted Pair-Groups Method Average.
double averageLinkage(int len, double *dis)
{
	double	avedis = 0;
	for(int i = 0; i < len; i ++)	{
		avedis += dis[i];
	}
	return (avedis / double(len));
}

//COMPLETE LINKAGE CLUSTERING (Maximum or Furthest-Neighbour Method): 
//The dissimilarity between 2 groups is equal to the
//greatest dissimilarity between a member of cluster i and a member of cluster j. 
double completeLinkage(int len, double *dis)
{
	double	maxdis = dis[0];
	for(int i = 1; i < len; i ++)	{
		if(maxdis < dis[i])	maxdis = dis[i];
	}	
	return maxdis;
}

//SINGLE LINKAGE CLUSTERING (Minimum or Nearest-Neighbour Method): 
//The dissimilarity between 2 clusters is the minimum
//dissimilarity between members of the two clusters.
double singleLinkage(int len, double *dis)
{
	double	mindis = dis[0];
	for(int i = 1; i < len; i ++)	{
		if(mindis > dis[i])	mindis = dis[i];
	}	
	return mindis;
};

#endif
