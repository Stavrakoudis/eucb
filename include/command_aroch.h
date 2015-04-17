#ifndef __COMMAND_AROCH__H
#define __COMMAND_AROCH__H
# include <command.h>
# include <command_stack.h>

class	CommandAroch: public Command
{
	private:
		vector<HeavySideChain> SideChain;
		Data	stack_percent;
		Data	stack_distance;
		Data	stack_angle;
	public:
		CommandAroch(Data sp,Data sd,Data sa);
		virtual void Run();
};

#endif
