# ifndef __COMMANDHYDROSOLVE__H
# define __COMMANDHYDROSOLVE__H

# include <command.h>
class CommandHydroSolve :public Command
{
	private:
		int ishydrophobic(atom a);
		Data	percent;
		Data	distance;
		Data	angle;
	public:
		CommandHydroSolve();
		void	setPda(Data p,Data d,Data a);
		virtual void Run();
		~CommandHydroSolve();
};

# endif
