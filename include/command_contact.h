# ifndef __COMMAND_CONTACT__H
# define __COMMAND_CONTACT__H
# include <command.h>

class CommandContact: public Command
{
	private:
		string selection1;
		string selection2;
	public:
		CommandContact();
		void	setSelection1(string s);
		void	setSelection2(string s);
		virtual void Run();
		~CommandContact();
};
# endif
