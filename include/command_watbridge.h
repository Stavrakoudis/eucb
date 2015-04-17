# ifndef __COMMAND_WATBRIDGE__H
# define __COMMAND_WATBRIDGE__H
# include <command.h>

class CommandWatbridge: public Command
{
	public:
		CommandWatbridge();
		virtual void Run();
		~CommandWatbridge();
};
# endif
