#ifndef __COMMAND_AHELIX__H
#define __COMMAND_AHELIX__H
# include <command.h>

class	CommandHelix: public Command
{
	private:
		Data	percent;
	public:
		CommandHelix(Data p);
		virtual void Run();
};

#endif
