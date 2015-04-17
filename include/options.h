# ifndef __OPTIONS__H
# define __OPTIONS__H
# include <util.h>
# include <string>
# include <vector>
using namespace std;
typedef void (*OptCallBack)(string);

class Options
{
	private:
		/*
 		Variables:
		=======================================================
		option: The name of the supported options.
		option_help: The help strings of each option.
		optarg: The functions that will be called when the 
		corresponding option will be activated.
		Methods:
		=======================================================
		Options(): The constructor of the class.
		addOption(): Insert a new option to the list of 
		supported options.
		parse(): Parse the command line and search for 
		available options.
		printHelp(): Print help information about the 
		supported options.
		~Options(): The destructor of the class.
 		* */
		vector<string> option;
		vector<string> option_help;
		vector<OptCallBack> optarg;
	public:
		Options();
		void	addOption(string s,OptCallBack opt,string help_string);
		int	parse(char *cmdline);
		void	printHelp();
		void	printHelp(string s);
		~Options();
};

# endif
