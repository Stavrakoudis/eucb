# ifndef __COMMANDENERGY__H
# define __COMMANDENERGY__H

# include <command.h>
class CommandEnergy :public Command
{
	private:
		vector<string> list;
	public:
		CommandEnergy(vector<string> &s);
		virtual void Run();
		~CommandEnergy();
};

# endif
