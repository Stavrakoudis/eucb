# ifndef __COMMAND_SALT__H
# define __COMMAND_SALT__H
# include <command.h>

class 	CommandSalt: public Command
{
	private:
		/*
 			Variables:
			===============================================
			poscsalt: The positive salt pairs of (res,atom)
			negcsalt; The negative salt pairs of (res,atom)
			Methods:
			===============================================
			CommandSalt: The constructor of the class.
			Run:	     Implement the command.
 		* */
		vector<ComplexSaltStruct> poscsalt;
		vector<ComplexSaltStruct> negcsalt;
	public:
		CommandSalt();
		virtual void Run();
};
# endif
