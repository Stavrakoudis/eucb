# include <command.h>
# include <time.h>
# include <unistd.h>

Command::Command(string s)
{
	name = s;
	logname = s+".log";
	fplog = fopen(logname.c_str(),"w");
	if(!fplog) psf_error("Can not open "+logname+" for writting");
	team = NULL;
	team1=NULL;
	team2=NULL;
	time_t rawtime;
  	struct tm * timeinfo;

 	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
	char str[200];
  	sprintf(str, "#Current local time and date: %s", asctime (timeinfo) );	
	str[strlen(str)-1]=0;
	Log(str);
	char currentPath[FILENAME_MAX];
	getcwd(currentPath,sizeof(currentPath));
	sprintf(str,"#Current directory: %s",currentPath);	
	Log(str);
	extern int Argc;
	extern char **Argv;
	string st=Argv[0];
	for(int i=1;i<Argc;i++) st=st+" "+Argv[i];
	Log("#Command line: "+st);
	time(&time_start);
}

void	Command::Log(string message)
{
	PrintLine(fplog,message);
}

void	Command::setSearch(vector<string> &s)
{
	if(team!=NULL) 
		team->setSeqString(s);
	else	
		team = new Team(s);
}

void	Command::setSearch1(vector<string> &s)
{
	if(team1!=NULL) 
		team1->setSeqString(s);
	else	
		team1 = new Team(s);
}

void	Command::setSearch2(vector<string> &s)
{
	if(team2!=NULL) 
		team2->setSeqString(s);
	else	
		team2 = new Team(s);
}

void  	Command::setFile(string s)
{
	vector<string> st;
	char file[1024];
	strcpy(file,s.c_str());
	if(team!=NULL) delete team;
	team = new Team(file,st);
}

void  	Command::setRes(vector<string> &s)
{
	if(team!=NULL) team->setResList(s);
}

void  Command::setTeam(Team *t)
{
	team = t;
}

string  Command::getName()
{
	return name;
}

Team	*Command::getTeam()
{
	return team;
}

Team	*Command::getTeam1()
{
	return team1;
}

Team	*Command::getTeam2()
{
	return team2;
}

void	Command::Error(string s)
{
	Log(s);
	psf_error(s);
}

Command::~Command()
{
	if(team!=NULL) delete team;
	if(team1!=NULL) delete team1;
	if(team2!=NULL) delete team2;
	time(&time_end);
	Data dif_time = difftime(time_end,time_start);
	Log("\n#Execution took "+printNumber(dif_time,10,2)+" seconds to complete");
	fclose(fplog);
}
