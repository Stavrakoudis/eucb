#ifndef __COMMAND_STACK__H
#define __COMMAND_STACK__H
# include <command.h>

typedef struct
{
	string res_name;
	string atom_name;
}HeavySideChain;

class	CommandStack: public Command
{
	private:
		vector<HeavySideChain> SideChain;
		Data	stack_percent;
		Data	stack_distance;
		Data	stack_angle;
	public:
		CommandStack(Data sp,Data sd,Data sa);
		virtual void Run();
};

#endif
