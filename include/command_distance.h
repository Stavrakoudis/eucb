# ifndef __COMMAND__DISTANCE__H
# define __COMMAND__DISTANCE__H
# include <command.h>

class CommandDistance :public Command
{
	private:
		/*
 		Variables:
		=============================================================
		Methods:
		=============================================================
		CommandDistance(): The default constructor.
		setFile1():        Read the elements of team1 from a file.
		setFile2():	   Read the elements of team2 from a file.
		setTeam1():	   Change the first team.
		setTeam2():	   Change the second team.
		setSearch1():	   Improve the first team.
		setSearch2():	   Improve the second team.
		setRes1():	   Select the residues for the first team.
		setRes2():	   Select the residues for the second team.
		Run():		   Perform the execution of the command.
		~CommandDistance(): Free the allocated memory.
 		* */
	public:	
		CommandDistance();
		void  setFile1(string s);
		void  setFile2(string s);
		void  setTeam1(Team *t);
		void  setTeam2(Team *t);
		void  setRes1(vector<string> &s);
		void  setRes2(vector<string> &s);
		virtual void Run();
		~CommandDistance();
};
# endif
