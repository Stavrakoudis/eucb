# ifndef __COMMAND__COMPLEXSALT__H
# define __COMMAND__COMPLEXSALT__H
# include <command.h>


class	CommandComplexSalt : public Command
{
	private:
		/*
 			Variables:
			===============================================
			possalt: The positive c-salt pairs of (res,atom)
			negsalt; The negative c-salt pairs of (res,atom)
			Methods:
			===============================================
			CommandSalt: The constructor of the class.
			Run:	     Implement the command.
 		* */
		vector<ComplexSaltStruct> possalt;
		vector<ComplexSaltStruct> negsalt;
	public:
		CommandComplexSalt();
		virtual void Run();
};

# endif
