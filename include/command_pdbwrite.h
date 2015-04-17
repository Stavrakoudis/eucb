# ifndef __COMMAND__PDB__WRITE__H
# define  __COMMAND__PDB__WRITE__H

# include <command.h>

class CommandPdbWrite : public Command
{
	public:
		CommandPdbWrite();
		virtual void Run();
};

# endif

