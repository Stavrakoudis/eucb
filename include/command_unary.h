# ifndef __COMMAND_UNARY__H
# define __COMMAND_UNARY__H
# define UNARY_DISTANCE		1
# define UNARY_ANGLE		2
# define UNARY_DIHEDRAL		3
# include <command.h>

class CommandUnary :public Command
{
	private:
		/*
 		Variables:
		=======================================================
 		flag: Determine the unary opeation to be performed.
		unary_atom: The list of atoms in the unary operation.
		Methods:
		=======================================================
		CommandUnary():  The constructor of the class.
		discover_atom(): Find the atom id of the given atom.
		unary_print():   Print the needed dat, stat, sda and
		hist files for the unary operation.
		unary_distance:  It called when the unary operation is
		a distance.
		unary_angle:  It called when the unary operation is
		a torsion.
		unary_dihedral:  It called when the unary operation is
		a dihedral angle.
		Run():  Implement the unary operation.
 		* */
		int flag;
		vector<string> unary_atom;
		int	discover_atom(string s);
		void	unary_print(string filename,string statname,string histname,string smoothname,
			IntVector f,
			DoubleVector x,int is_angle);
		void	unary_distance();
		void	unary_angle();
		void	unary_dihedral();
	public:
		CommandUnary(int Flag,vector<string> &t);
		virtual void Run();
};
# endif
