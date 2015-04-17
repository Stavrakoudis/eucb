# ifndef __COMMAND__CENTER_H
# define __COMMAND__CENTER_H
# include <command.h>

class CommandCenter :public Command
{
	private:
		/*
 		Variables:
		=============================================================
		Methods:
		=============================================================
		CommandCenter(): The default constructor.
		Run():		   Perform the execution of the command.
		~CommandCenter(): Free the allocated memory.
 		* */
		int	ca_flag1;
		int	noh_flag1;
		int	sidechain_flag1;
		int	all_flag1;

		int	ca_flag2;
		int	noh_flag2;
		int	sidechain_flag2;
		int	all_flag2;
	public:	
		CommandCenter();
		void	setFlags1(vector<string> &flags1);
		void	setFlags2(vector<string> &flags2);
		virtual void Run();
		~CommandCenter();
};
# endif
