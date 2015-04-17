# ifndef __COMMAND_RMSD__H
# define __COMMAND_RMSD__H
# include <command.h>

class CommandRmsd : public Command
{
	private:
		/*
		CommandRmsd(): The constructor of the class.
		Run():  Implement the command.	
 		* */
	public:
		CommandRmsd();
		virtual  void Run();
};

# endif
