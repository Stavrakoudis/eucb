# ifndef __TEX__H
# define __TEX__H
# include <psf.h>
# include <pdb.h>

typedef struct
{
	int first;
	int second;
	double percent;
}TexGroup;

class Tex
{
	private:
		vector<string> 		header;
		IntVector		isgroup;
		vector<IntVector> 	Atom;
		DoubleVector		percent;
		string			filename;
		string			caption;
		void			sort();
		int			enable_atoms;
	public:
		Tex(string name);
		void	enableAtoms(int f);
		void	setCaption(string s);
		void	setHeader(vector<string> s);
		void	addAtomList(IntVector a,Data p);
		void	addGroupList(vector<TexGroup> &g);
		void	print();
};

# endif
