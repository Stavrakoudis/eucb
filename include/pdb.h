# ifndef __PDB__H
# define __PDB__H
# include <util.h>

typedef struct
{
	int	atom_id;
	string	atom_name;
	string	res_name;
	string  pdb_chain;
	int	res_id;
	Data	x;
	Data	y;
	Data	z;
	Data	val1;
	Data	val2;
	string	chain_id;
	string	pdb_letter;	
}PdbAtom;

extern vector<PdbAtom> pdb_table;

extern	Data	getPdbDistance(PdbAtom a,PdbAtom b);

extern string  pdb_name;

/*
 * 	pdbpos: A pointer holding in (x,y,z) format
 * 	the position of the atoms in the pdb file.
 * */
extern PointVector pdbpos;

/*
 * 	Read and parse the pdb file.
 * */
extern int readPdb(string filename);
# endif
