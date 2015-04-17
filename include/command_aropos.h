#ifndef __COMMAND_AROPOS__H
#define __COMMAND_AROPOS__H
# include <command.h>
# include <command_stack.h>

typedef struct
{
	string chain;
	string nname;
	vector<string> hname;
}AroposNames;

class	CommandAropos: public Command
{
	private:
		vector<HeavySideChain> SideChain;
		vector<AroposNames>    posnames;
		Data	stack_percent;
		Data	stack_distance;
		Data	stack_angle;
	public:
		CommandAropos(Data sp,Data sd,Data sa);
		virtual void Run();
};

#endif
