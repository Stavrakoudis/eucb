#ifndef __COMMAND_ARONH__H
#define __COMMAND_ARONH__H
# include <command.h>
# include <command_stack.h>

class	CommandAronh: public Command
{
	private:
		vector<HeavySideChain> SideChain;
		Data	stack_percent;
		Data	stack_distance;
		Data	stack_angle;
	public:
		CommandAronh(Data sp,Data sd,Data sa);
		virtual void Run();
};

#endif
