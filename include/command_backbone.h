# ifndef __COMMAND_BACKBONE__H
# define __COMMAND_BACKBONE__H
# include <command.h>
class CommandBackbone :public Command
{
	/*
 		Methods:
		=======================================================
		CommandBackbone(): The constructor of the class.
		Run(): Execute the command.
 	* */
	public:
		CommandBackbone();
		virtual void Run();
};
# endif
