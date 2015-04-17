# ifndef __COMMANDWRITEDCD__H
# define __COMMANDWRITEDCD__H
# include <command.h>
class CommandWriteDcd :public Command
{
	private:
		string fname;
		int	fit;
	public:
		CommandWriteDcd(string s,int f);
		virtual void Run();
		~CommandWriteDcd();
};

# endif
