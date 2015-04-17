# ifndef __COMMAND_NOE__H
# define __COMMAND_NOE__H
# include <command.h>

class CommandNoe :public Command
{
	private:
		/*
 		*	Variables:
 		*	===================================================
 		*	noe_atom: The list of atom names belonging to 
 		*	the noe list  (usually HA, HN)
 		*	average:
 		*	Used in the statistics part of noe. The default 
 		*	value for this parameter is 6.0
 		*	Methods:
 		*	================================================== 
 		*	CommandNoe(): The constructor of the class.
 		*	setAtoms(): Change the list of atom names.
 		*	setAverage(): Change the value of average.
 		*	Run():       Implement the command.
 		* */
		vector<string> noe_atom;
		Data average;
	public:
		CommandNoe();	
		void	setAtoms(vector<string> &s);
		void	setAverage(Data d);
		virtual void Run();
};
# endif
