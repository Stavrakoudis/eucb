# ifndef __COMMAND_CLUSTERRMSD__H
# define __COMMAND_CLUSTERRMSD__H
# include <command.h>

class CommandClusterRmsd : public Command
{
	private:
		double	limit;
	public:
		CommandClusterRmsd(double l);
		virtual  void Run();
};

# endif
