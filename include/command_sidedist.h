# ifndef __COMMAND_SIDEDIST__H
# define __COMMAND_SIDEDIST__H
# include <command.h>

class CommandSideDist: public Command
{
	private:
		/*
 		Variables:
		=======================================================
		Methods:
		CommandSideDist(): The constructor of the class.
		setSearch1(): Enhance the first sequence.
		setSearch2(): Enhance the second sequence.
		setName1():  Set the name of the first sequence.
		setName2():  Set the name of the second sequence.
		Run():       Implement the command.
		~CommandSideDist(): The destructor of the class.
 		* */
	public:
		CommandSideDist();
		void	setName1(string s);
		void	setName2(string s);
		virtual void Run();
		~CommandSideDist();
};
# endif
