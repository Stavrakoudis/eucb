# ifndef __COMMAND_MAKE__RMSF__H
# define __COMMAND_MAKE__RMSF__H
# include <command.h>

class CommandRmsf
	:public Command
{
	public:
		CommandRmsf();
		void printRMSF(string chain);
		virtual void Run();
};

# endif
