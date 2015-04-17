# ifndef __COMMAND_ANAL__H
# define __COMMAND_ANAL__H
# include <command.h>

class 	CommandAnal: public Command
{
	public:
		CommandAnal();
		virtual void Run();
};
# endif
