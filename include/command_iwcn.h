# ifndef __COMMAND_IWCN__H
# define __COMMAND_IWCN__H
# include <command.h>

class CommandIwcn : public Command
{
	private:
		double	distance,percent;
		int	maxn;
	public:
		CommandIwcn(double p,double d,int m);
		virtual void Run();
};
# endif
