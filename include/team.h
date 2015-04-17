# ifndef __TEAM__H
# define __TEAM__H
# include <psf.h>
# include <util.h>
# include <vector>
# include <string>
using namespace std;

/*	
 *	TeamStruct: Determines the members of a sequence. Sequences are used
 *	in most commands such as group, tors, hbonds etc. 
 *	mask: Determines which fields are used.
 *	chain_id: The chain of the sequence.
 *	res_id: The residue number of the sequence.
 *	res_name: The residue name of the sequence.
 *	atom_name: The atom name of the sequence.
 * */
struct TeamStruct
{
	int 	mask;
	string 	chain_id;
	int 	res_id;
	string 	res_name;
	string	atom_name;
};
	
class Team;

extern Team*	makeHydrophobic();
extern Team*	makeAromatic();
extern Team*	makeAliphatic();
extern Team*	makePositive();
extern Team*	makeNegative();

class Team
{
	private:
		/*
 		Variables:
		=======================================================
		table: All the parts that are composing the sequence.
		res_list: The residue names of the sequence.
		Methods:
		=======================================================
		Team(): The class can be constructed either from a file
		or from a list of strings.
		setSeqString(): Change the sequence string.
		setResList(): Set the names of the residues.
		find(): Return 1 or 0 if the given atom is 
		in the sequence.
		replace(): replace the sequence strings with s.
		print(): Print the sequence to the standard output.
		enumerateChains(): Return all the chain_ids 
		in the sequence.
		enumerateResId(): Return all the residue numbers
		in the chain specified by the first argument.
		enumerateAtoms(): Return all the atoms ids that
		are in the current sequence.
		~Team(): The constructor of the class.
 		* */
		vector<TeamStruct> table;
		vector<string> res_list;
		vector<string> res_name;
		vector<string> atom_name;
		vector<string> atom_type;
	public:
		Team(char *file,char *str1,char *str2);
		Team(char *file,vector<string> list);
		Team(vector<string> list);
		Team(vector<string> arg1,vector<string> arg2);
		void	setSeqString(vector<string> list);
		void	setResList(vector<string> &s);
		void	setAtomType(vector<string> &s);
		void	setResName(vector<string> &s);
		void	setAtomName(vector<string> &s);
		int 	find(atom atom,string &chain);
		void	replace(vector<string> &s);
		void 	print();
		void	enumerateChains(vector<string> &result);
		void	enumerateResId(string s,vector<int> &res);
		void	enumerateAtoms(vector<atom> &list,vector<int> &res);
		int	getAtomNameSize();
		~Team();
};

# endif
