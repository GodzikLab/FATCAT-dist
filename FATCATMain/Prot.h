//-----------------------------------------------
//The simple protein class 
//only with CA coordinations, and solvent accessible information
//by Y.Y, date 10/24/02
//latest update 10/24/02
//-----------------------------------------------
#ifndef PROT_H
#define PROT_H

#include "basic.h"

using namespace ARRAY;
using namespace GEOMETRY;

class MULTAFPCHAIN;
class POSTALIGN;
class LSVECTOR;

class PROT	{
	protected:
		int	length;	  //length of the protein
		char	**res;    //the residue list in AAA
		char	*saa;     //the residue list in A
		char	*chain;   //the chain of residues
		char	**index;  //the original index of each position in pdb file
		double	*caCod;   //the CA coordination, 3 sizes of the length (x,y,z)	
		double	*asa;     //the solvent accessible of residues
		int	*res2Atm; //residue to atom index

		int	totalAtm; //the number of total atoms
		double	*atmCod;  //coordinations of all atoms
		char	**atmName;//name of atoms
		int	*atm2Res; //atom to residue index

		int	length0;
		int	atm0;

		void	Initial(void);
		void	Calloc(int len);
		void	CallocAtm(int atm);
		void	AssignPdbCa(char *file);  //read processed pdb file with only CA coords.
	public:
		PROT(char *filename);
		PROT(char *filename, int begres, int reslen);
		PROT(void);
		~PROT(void);
		void	SAA(void);
		void	AssignPdb(char *file, int begres, int reslen); //read original pdb file
		void	AssignPdb(PROT *pro); //read original pdb file
		double* Cod4Res(int n, int *list); //return the coordinations for a given residue list
		void	PrintPdb(char forcechain, ofstream &stream);
		double	CalCaRmsd(PROT *pro, int resn, int *res1, int *res2);
		void 	TransCod(double *r, double *t);
		void 	TransCod(double *r, double *t, int r1, int r2);
		int	GetLength(void);
		void	GetSAA(char *saa_o);
		void	GetIndex(int *index_o);
		int	GetIndex(int index_i);
		int	ResIndex(char *index_i);
		double**	DisTable(int maxlen);
		friend class AFPCHAIN;
		friend class LSVECTOR;
		friend class POSTALIGN;
};

#endif
