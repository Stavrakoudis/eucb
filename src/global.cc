# include <command.h>
# include <limits.h>

int	donor_act=0;
int 	acceptor_act=0;
int	donor_or_acceptor_first=0;//1 for donor and 2 for acceptor

vector<atom> table;
vector<string> res_list;

string 	psffile="";
string 	dcdfile="";
string	pdbfile="";
Data	critical_distance=5.5;//A^0
Data  critical_percent=0.05;//5%
Data	critical_angle=120.0;

int	backbone_flag=0;
Data	backbone_critical_distance=3.4;
Data	backbone_critical_percent=0.01;
Data	backbone_critical_angle=120.0;

int	smart_skip = 10;
Data	smart_distance = 10.0;

int	first=1;
int	last=INT_MAX;
int	step=1;
Dcd	*dcd;
int	isbackbone=0;
string    donor_selection="all";
string    acceptor_selection="all";
Team 	*donorTeam=NULL;
Team	*acceptorTeam=NULL;

/*	Definitions for the -watbridge parameter
 * */
int	watbridge_flag=0;
int 	watbridge_level=1;

/*	Definitions for the -noe command
 * */
Data	noe_critical_distance=5.0,noe_critical_percent=0.2;
/*	Declarations for the hist command.
 * */
int histflag=1;
Data bindist=0.5;
Data binangle=10.0;
Data bindihe=10.0;

class Command;

vector<Command*> command_list;

/*	Declarations for the smooth command
 * */
int smooth_flag=0;
int smooth_start=50;
int smooth_step=-1;
int weight_flag = 0;
int contact_flag =0;
Data contact_distance=4.0;
Data contact_percent=0.05;
int	res_diff=2;

int astart  = 0;
int kstep   = 1;
int astep   = 1;

int statcount=10;
Data	diel=1.0;

string	parfile="";

int	zero_angle_flag = 0;	
int	degrees_angle_flag = 0;
int	ref=0;
int 	getAtomPos(int res_id,string atom_name)
{
	for(int i=0;i<table.size();i++) 
		if(table[i].res_id==res_id && table[i].atom_name==atom_name)
			return table[i].atom_id;
	return 0;
}

