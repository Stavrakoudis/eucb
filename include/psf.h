# ifndef __PSF__H
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <util.h>
# include <rad.h>
# include <vector>
using namespace std;

/*	
 *	A structure for the bonding between atoms.
 * */
typedef struct
{
	int first;
	int second;
}BondStruct;


typedef struct
{
	int 	first;
	int	second;
	int	third;
}AngleStruct;

typedef struct
{
	int	first;
	int	second;
	int	third;
	int	fourth;
}DihedralStruct;

/*	
 *	The basic eucb structure. It  holds all the needed 
 *	information to describe an atom.
 * */
typedef struct
{
	int	atom_id;
	string	chain_id;
	int	res_id;
	string  res_name;
	string  atom_name;
	string 	atom_type;
	Data 	dnum1;
	Data 	dnum2;
	int	lastnum;
	float	radious;
	vector<int>	bond;
	int	ishydro;
}atom;

extern vector<BondStruct> bondstruct;
extern vector<AngleStruct> anglestruct;
extern vector<DihedralStruct> dihedralstruct;
extern vector<DihedralStruct> improperstruct;
/*	
 *	Read the psffile and store the information for the atoms in 
 *	the PSF file to the table.
 * */
extern int	parse_psf(char *psffile,vector<atom> &table);

/*
	Write a portion of the table to the file fname.
*/
extern	void	write_psf(string fname,vector<atom> &table,IntVector &posAtom);

/*	
 *	Return 1 if the atoms with id1 and id2 are connected.
 * */
extern int 	isConnected(vector<atom> &table,int id1,int id2);

/*
 * 	Return the number of atoms that are connected with atom
 * 	with atom_id id.
 * */
extern int	countConnected(vector<atom> &table,int id);

/*
 * 	List all the chain_ids of the psf.
 * */
extern void	enumerateChains(vector<atom> &table,vector<string> &result);

/*
 * 	Return 1 if the atom is hydrogen.
 * */
extern	int	isHydro(atom a);

/*
 * 	Return 1 if the atom is backbone.
 * */
extern	int	isBackbone(atom a);

/*
 * 	Return 1 if the atom is backbone4.
 * */
extern  int	isBackbone4(atom a);

/*
 * 	Return 1 if the atom is backbone3.
 * */
extern  int 	isBackbone3(atom a);

/*
 * 	Return 1 if the atom is sidechain.
 * */
extern  int	isSidechain(atom a);

/*
 *	Return 1 if the atom is water atom.
 * */
extern  int	isWater(atom a);

/*
 *	Return 1 if the atom is positive atom.
 * */
extern int	isPositive(atom a);

/*
 *	Return 1 if the atom is positive atom.
 * */
extern int	isNegative(atom a);


/*
 *	Return 1 if the atom is  a hydrophobic atom.
 * */
extern	int 	isHydroPhobic(atom a);
/*	Print the atom information with various formats.
 *	This informations is usually used for the naming
 *	of the produced .dat and .stat files.
 * */
extern	 	string printAtom1(atom a);
extern 	 	string printAtom2(atom a);
extern	 	string printAtom3(atom a);
extern	 	string printAtom5(atom a);
extern		string printAtomHeader();
extern		string printAtom4(atom a);
extern		string printAtom6(atom a);
extern		string printAtomLatex(atom a,int f);
extern		string printAtomLatexGroup(atom a,int f);
# define __PSF__H
# endif
