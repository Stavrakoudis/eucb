# ifndef __COMMAND__PDO__H
# define __COMMAND__PDO__H
# include <command.h>

class CommandPdo : public Command
{
	private:
		int pdo_level;
	public:
		CommandPdo(int l);
		virtual void Run();
};
# endif
