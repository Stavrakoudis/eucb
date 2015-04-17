# ifndef __COMMAND_ACCEPTOR__H
# define __COMMAND_ACCEPTOR__H
# include <command.h>
class CommandAcceptor :public Command
{
	public:
		CommandAcceptor();
		virtual void Run();
};
# endif
