# ifndef __COMMAND__H
# define __COMMAND__H

# include <vector>
# include <psf.h>
# include <dcd.h>
# include <team.h>
# include <time.h>
# include <string>
using namespace std;

class Command
{
	protected:
		/*
 		Variables:
		=============================================================
 			name   :The name of the command.
			logname:The name of the log file.
			fplog  :The file pointter for the log file.
			team   :The sequence used for the current command. 
		Methods:
		=============================================================
		Command(): The construtor of the class. It sets the name of 
		the class and it opens the log file.
		setSearch(): Specify the sequence for the team of the class.
		setFile(): Read the elements of the team from a file.
		setRes():  Specify the residues for the team of the class.
		setTeam(): Change the team of the class.
		Run():    A virtual function with the code that will be 
		implemented by the class.
		getName(): Return the name of the class.
		~Command(): The destructor of the class. It frees the 
		memory of the team and it closes the log file.
 		* */
		string name;
		string logname;
		FILE   *fplog;
		Team	*team;
		Team	*team1;
		Team	*team2;
		time_t	time_start,time_end;
		void	Log(string message);
	public:
		Command(string s);
		void	setSearch(vector<string> &s);
		void	setSearch1(vector<string> &s);
		void	setSearch2(vector<string> &s);
		void  	setFile(string s);
		void  	setRes(vector<string> &s);
		void	setTeam(Team *t);
		Team	*getTeam();
		Team	*getTeam1();
		Team	*getTeam2();
		virtual void Run()=0;
		string  getName();
		void	Error(string s);
		~Command();
};

# endif
