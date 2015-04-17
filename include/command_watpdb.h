# ifndef __COMMAND__WATPDB__H
# define  __COMMAND__WATPDB__H

# include <command.h>

class CommandWatpdb : public Command
{
	private:
		Data distance;
		void PrintAtomInPdb(FILE *fp,DoubleVector2 &D,int id);
	public:
		CommandWatpdb(Data d);
		virtual void Run();
};

# endif

