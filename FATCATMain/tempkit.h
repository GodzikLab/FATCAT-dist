//-------------------------------------------------------------------------
//Y.Y Jun/03/04 latest update 
//The macro-like functionality of templates, forces a restriction for multi-file projects: 
//the implementation (definition) of a template class or function must be in the same file as the declaration. 
//That means we cannot separate the interface in a separate header file
//-------------------------------------------------------------------------
#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "basic.h"

using namespace std;

//-------------------------------------------------------------------------
//the following are templates for array calloc and deletion (one dim, two dims, and three dims)
//-------------------------------------------------------------------------
namespace ARRAY	{
	//-------------------------------------------------------------------
	//dim: the dimention of the array
	//return an array point
	//-------------------------------------------------------------------
	template <class Tdata>
	Tdata *NewArray(int dim)
	{
		Tdata *array = new Tdata[dim];
		for(int i = 0; i < dim; i ++)	array[i] = 0;
		return array;
	};

	//-------------------------------------------------------------------
	//delete an array
	//-------------------------------------------------------------------
	template <class Tdata>
	void DelArray(Tdata *array)
	{
		delete[] array; 
	};

	//-------------------------------------------------------------------
	//print an array
	//-------------------------------------------------------------------
	template <class Tdata>
	void PrintArray(int numEntries, Tdata *A)
	{   
    		for (int i = 0; i < numEntries; i++) {
            		cout<<A[i]<<endl;
        	}
        	cout<<endl; 
	};   

	//-------------------------------------------------------------------
	//new a matrix, give first dimention dim1, second dimention dim2
	//return a matrix point
	//-------------------------------------------------------------------
	template <class Tdata>
	Tdata **NewMatrix(int dim1, int dim2)
	{
		Tdata **matrix = new Tdata*[dim1];
		for(int i = 0; i < dim1; i ++)	{
			matrix[i] = new Tdata[dim2];
			for(int j = 0; j < dim2; j ++)	matrix[i][j] = 0;
		}
		return matrix;
	};

	//-------------------------------------------------------------------
	//delete a matrix
	//-------------------------------------------------------------------
	template <class Tdata>
	void DelMatrix(Tdata **matrix, int dim1)
	{
		for(int i = 0; i < dim1; i ++)	{ 
			delete[] matrix[i];
		}
		delete[] matrix; 
	};

	//-------------------------------------------------------------------
	//new a 3-dimentional array, give first, second and third dim
	//return a point
	//-------------------------------------------------------------------
	template <class Tdata>
	Tdata ***NewArray3(int dim1, int dim2, int dim3)
	{
		Tdata ***array3 = new Tdata**[dim1];
		for(int i = 0; i < dim1; i ++)	{
			array3[i] = new Tdata*[dim2];
			for(int j = 0; j < dim2; j ++)	{
				array3[i][j] = new Tdata[dim3];
				for(int k = 0; k < dim3; k ++)	array3[i][j][k] = 0;
			}
		}
		return array3;
	};

	//-------------------------------------------------------------------
	//delete a matrix
	//-------------------------------------------------------------------
	template <class Tdata>
	void DelArray3(Tdata ***array3, int dim1, int dim2)
	{
		for(int i = 0; i < dim1; i ++)	{ 
			for(int j = 0; j < dim2; j ++)	{
				delete[] array3[i][j];
			}
			delete[] array3[i];
		}
		delete[] array3; 
	};
}; //end of ARRAY


//----------------------------------------------------------------------------
//the following are templates for statistics 
//includes qksort, average-calculation, meansquare calculation, median value calculation
//rank a series, linearfit, ranksum test, maximum and minimum calculation et al.
//----------------------------------------------------------------------------
namespace STATIS	{
	template <class Tdata>
	void qksort(int ilo, int ihi, Tdata *A, int *sort) {
    		Tdata pivot;		// pivot value for partitioning array
    		int ulo, uhi;	// indices at ends of unpartitioned region
    		int ieq;		// least index of array entry with value equal to pivot
    		Tdata tempEntry;	// temporary entry used for swapping
    		int    tempSort;    // temporary entry used for recording sort

    		if (ilo >= ihi) {
			return;
    		}
    		// Select a pivot value.
    		pivot = A[(ilo + ihi)/2];
    		// Initialize ends of unpartitioned region and least index of entry
    		// with value equal to pivot.
    		ieq = ulo = ilo;
    		uhi = ihi;
    		// While the unpartitioned region is not empty, try to reduce its size.
    		while (ulo <= uhi) {
			if (A[uhi] > pivot) {
	    			// Here, we can reduce the size of the unpartitioned region and
	    			// try again.
	    			uhi--;
			} 
			else {
	    			// Here, A[uhi] <= pivot, so swap entries at indices ulo and uhi
	    			tempEntry = A[ulo];
	    			A[ulo] = A[uhi];
	    			A[uhi] = tempEntry;
	    			tempSort = sort[ulo];
	    			sort[ulo] = sort[uhi];
	    			sort[uhi] = tempSort;
	    			// After the swap, A[ulo] <= pivot.
	    			if (A[ulo] < pivot) {
					// Swap entries at indices ieq and ulo.
					tempEntry = A[ieq];
					A[ieq] = A[ulo];
					A[ulo] = tempEntry;
					tempSort = sort[ieq];
					sort[ieq] = sort[ulo];
					sort[ulo] = tempSort;
					// After the swap, A[ieq] < pivot, so we need to change
					// ieq.
					ieq++;
					// We also need to change ulo, but we also need to do
					// that when A[ulo] = pivot, so we do it after this if
					// statement.
	    			}
	    			// Once again, we can reduce the size of the unpartitioned
	    			// region and try again.
	    			ulo++;
			}
    		}
    		// Now, all entries from index ilo to ieq - 1 are less than the pivot
    		// and all entries from index uhi to ihi + 1 are greater than the
    		// pivot.  So we have two regions of the array that can be sorted
    		// recursively to put all of the entries in order.
    		qksort(ilo, ieq - 1, A, sort);
    		qksort(uhi + 1, ihi, A, sort);
	};

	//-------------------------------------------------------------------------
	template <class Tdata>
	void qksort(int numEntries, Tdata *A, int *sort) 
	{
    		qksort(0, numEntries - 1, A, sort);
	};

	//-------------------------------------------------------------------------
	template <class Tdata>
	void reverse(int numEntries, Tdata *A, int *sort)
	{
    		int	i;
    		int mid = int(numEntries / 2.0);
    		int tempSort;
    		Tdata tempEntry;
    		for(i = 0; i < mid; i ++)	{
			tempSort = sort[i];
        		sort[i] = sort[numEntries - 1 - i];	
			sort[numEntries - 1 - i] = tempSort;
			tempEntry = A[i];
        		A[i] = A[numEntries - 1 - i];	
			A[numEntries - 1 - i] = tempEntry;
    		}
	};

	template <class Tdata>
	Tdata max(int num, Tdata *data)
	{
		Tdata mv = data[0];
		for(int i = 1; i < num; i ++)	{
			if(data[i] > mv)	mv = data[i];
		}
		return mv;
	};

	template <class Tdata>
	Tdata min(int num, Tdata *data)
	{
		Tdata mv = data[0];
		for(int i = 1; i < num; i ++)	{
			if(data[i] < mv)	mv = data[i];
		}
		return mv;
	};

	template <class Tdata>
	void distribute(int num, Tdata *data, int step, double *scalelist, double *distrib)
	{
		Tdata minv = min(num, data);
		Tdata maxv = max(num, data);
		double	scale = (double(maxv) - double(minv))/double(step - 1);
		int	i, j;
		for(i = 0; i < step; i ++)	{
			scalelist[i] = minv + i * scale;
		}
		for(i = 0; i < step; i ++)	distrib[i] = 0;
		for(i = 0; i < num; i ++)	{
			j = int( (double(data[i]) - double(minv)) / scale); 
			distrib[j] += 1;
		}
	};

}; //end of STATIS

#endif
