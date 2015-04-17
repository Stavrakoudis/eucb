# ifndef __COMMAND_CLOSEWATER__H
# define __COMMAND_CLOSEWATER__H
# include <command.h>

class CommandCloseWater :public Command
{
	public:
		CommandCloseWater();
		virtual void Run();
		~CommandCloseWater();
};
# endif
