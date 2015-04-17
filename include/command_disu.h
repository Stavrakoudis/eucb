# ifndef __COMMANDDISU__H
# define __COMMANDDISU__H
# include <command.h>

class CommandDisu :public Command
{
	public:
		CommandDisu();
		virtual void Run();
		~CommandDisu();
};

# endif
