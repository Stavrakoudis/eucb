# include <options.h>
# include <psf.h>

/*	Options: A genetic class for the handling of command
 *	line options.
 * */


/*	The default constructor.
 * */
Options::Options()
{
}

/*	Add a new option to the option list. The first argument is the name of the option,
 *	the second argument is the function that will called when the option is activated
 *	and the last argument is the string that will printed when the user asks 
 *	information about the argument.
 * */
void	Options::addOption(string s,OptCallBack opt,string help_string)
{
	int k=option.size();
	option.push_back(s);
	optarg.push_back(opt);
	option_help.push_back(help_string);
}

/*	Print a help screen.
 * */
void	Options::printHelp()
{
	for(int i=0;i<option.size();i++)
	{
		printf("%s\n",option_help[i].c_str());
	}
}

/*	Print a help about a command.
 * */
void	Options::printHelp(string s)
{
	for(int i=0;i<option.size();i++)
		if(option[i]==s) 
		{
			printf("%s\n",option_help[i].c_str());
			return ;
		}
	psf_error("Undefined option "+s);
}

/*	Parse the command line given and store the arguments into the internal array.
 * */
int	Options::parse(char *cmdline)
{
	int index=0;
	char word[1024];
	int icount=0;
	string lastoption="";
	string arg="";
	while(getword(cmdline,word,index))
	{
		if(word[0]=='-')
		{
			if(lastoption!="")
			{
			arg=word;
			int found=0;
			for(int i=0;i<option.size();i++)
			{
				if(option[i]==lastoption)
				{
					if(optarg.size()<=i)
					{
						psf_error("THE NUMBER OF CALLBACKS IS NOT SUFFICIENT ");
					}
					optarg[i](arg);
					found=1;
					break;
				}
			}
			if(!found)
			{
				string s1="THE OPTION "+lastoption+" IS NOT DEFINED";
				psf_error(s1);
			}
			}
			lastoption=word;
		}
		else
		{
			arg=word;
			int found=0;
			for(int i=0;i<option.size();i++)
			{
				if(option[i]==lastoption)
				{
					if(optarg.size()<=i)
					{
						psf_error("THE NUMBER OF CALLBACKS IS NOT SUFFICIENT ");
					}
					optarg[i](arg);
					found=1;
					break;
				}
			}
			if(!found)
			{
				string s1="THE OPTION "+lastoption+" IS NOT DEFINED";
				psf_error(s1);
			}
			lastoption="";
		}
	}	
	if(lastoption!="")
	{
		arg=word;
		int found=0;
		for(int i=0;i<option.size();i++)
		{
			if(option[i]==lastoption)
			{
				if(optarg.size()<=i)
				{
					psf_error("THE NUMBER OF CALLBACKS IS NOT SUFFICIENT ");
				}
				optarg[i](arg);
				found=1;
				break;
			}
		}
		if(!found)
		{
				string s1="THE OPTION "+lastoption+" IS NOT DEFINED";
				psf_error(s1);
		}
	}
}

Options::~Options()
{
}
