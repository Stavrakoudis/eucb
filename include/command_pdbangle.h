# ifndef __COMMAND_PDBANGLE__H
# define __COMMAND_PDBANGLE__H
# include <command.h>

class	CommandPdbAngle :public Command
{
	private:
		int	phi_flag;
		int	psi_flag;
		int	omega_flag;
		int	chi1_flag;
	public:
		CommandPdbAngle();
		void	setAngleList(vector<string> &s);
		virtual void Run();
};
# endif
