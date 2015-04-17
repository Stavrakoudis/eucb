# ifndef __COMMANDSASA__H
# define __COMMANDSASA__H
# include <command.h>
class CommandSasa :public Command
{
	private:
		Data	srad;
	public:
		CommandSasa(Data s);
		virtual void Run();
		~CommandSasa();
};

# endif
