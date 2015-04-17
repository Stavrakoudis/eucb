# ifndef __COMMAND_HBONDS1__H
# define __COMMAND_HBONDS1__H
# include <command.h>
# include <tex.h>
class CommandHbonds1 :public Command
{
	public:
		CommandHbonds1();
		virtual void Run();
		void	SubRun(Team *team1,Team *team2,Tex &myTex);
};
# endif
