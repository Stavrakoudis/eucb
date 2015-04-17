# ifndef __COMMAND__BTURN__H
# define __COMMAND__BTURN__H
# include <command.h>

# define THETA1		30
# define THETA2		45
# define CUTOFF_THETA	0.0

typedef struct
{
	string type;
	Data phi2;
	Data psi2;
	Data phi3;
	Data psi3;
	Data	distance;	
	int	counter;
}BetaStruct;

class CommandBturn
	:public Command
{
	private:
		vector<BetaStruct> Beta;
		Data	bturn_distance;
		Data  bturn_angle;
		Data 	bturn_percent;
	public:
		CommandBturn(Data P,Data D,Data A);
		virtual void Run();
};
# endif
