# include <global.h>
# include <options.h>
# include <command_distance.h>
# include <command_bturn.h>
# include <command_complexsalt.h>
# include <command_rmsf.h>
# include <command_angle.h>
# include <command_noe.h>
# include <command_stack.h>
# include <command_salt.h>
# include <command_sidedist.h>
# include <command_anal.h>
# include <command_backbone.h>
# include <command_donor.h>
# include <command_watbridge.h>
# include <command_unary.h>
# include <command_aronh.h>
# include <command_aropos.h>
# include <command_aroch.h>
# include <command_rmsd.h>
# include <command_center.h>
# include <command_pdbwrite.h>
# include <command_pdbangle.h>
# include <command_pdo.h>
# include <command_contact.h>
# include <command_hydrosolve.h>
# include <command_closewater.h>
# include <command_iwcn.h>
# include <command_ahelix.h>
# include <command_watpdb.h>
# include <command_clusterrmsd.h>
# include <command_entropy.h>
# include <command_writedcd.h>
# include <command_energy.h>
# include <command_sasa.h>
# include <command_disu.h>
# include <command_hbonds1.h>
# include <helpfiles.h>
# include <iostream>
Team	*LastTeam=NULL;
string	last_team="";

typedef struct
{
	string name;
	string description;
}HelpStruct;

vector<HelpStruct> HelpTable;


void	parse_help()
{
	for(int icount=0;icount<eucbhelp_table.size();icount++)
	{
		char line[1024];
		strcpy(line,eucbhelp_table[icount].c_str());
		char word[1024];
		int index=0;
		getword(line,word,index);
		int s=HelpTable.size();
		if(!strcmp(word,"OPTION")) 
		{
			getword(line,word,index);
			HelpTable.resize(s+1);
			HelpTable[s].name = word;
			HelpTable[s].description="";
		}
		else
		{
			string ss=line;
			if(HelpTable[s-1].description=="")
				HelpTable[s-1].description="\t"+ss;
			else
				HelpTable[s-1].description+="\n\t"+ss;
		}
	}
}

int	helpPos(string name)
{
	for(int i=0;i<HelpTable.size();i++)
		if(HelpTable[i].name == name) return i;
	return 0;
}

string	helpDescription(string name)
{
	int p=helpPos(name);
	return "OPTION: "+name+"\n"+"\tDESCRIPTION: "+HelpTable[p].description;
}

Options opt;
int     ilast=0;
int 	team1_replace=0;
int	team2_replace=0;
int	angle_flag;
int	isnoe=0;

int 	commandPos(string name)
{
	for(int i=0;i<command_list.size();i++)
		if(command_list[i]->getName() == name) return i;
	return -1;
}


void	BreakCritical(string s)
{
	vector<string> str;
	getCommaElements(s,str);
	for(int i=0;i<str.size();i++)
	{
		string ss=str[i];
		if(i==1)
		{
			if(isbackbone) backbone_critical_distance=atof(ss.c_str());
			else
			if(isnoe)      noe_critical_distance=atof(ss.c_str());
			else
			if(contact_flag) contact_distance=atof(ss.c_str());
			else
				       critical_distance=atof(ss.c_str());
		}
		else
		if(i==0)
		{
			if(isbackbone) backbone_critical_percent=atof(ss.c_str());
			else
			if(isnoe)      noe_critical_percent = atof(ss.c_str());
			else
			if(contact_flag) contact_percent = atof(ss.c_str());
			else
			critical_percent=atof(ss.c_str());
		}
		else
		{
			if(isbackbone) backbone_critical_angle=atof(ss.c_str());
			else	       critical_angle=atof(ss.c_str());
		}
	}
	isbackbone=0;
	isnoe = 0;
	contact_flag=0;
}

string	getLast(string s,int n)
{
	string ret="";
	int ipos=0;

	for(ipos=0;ipos<s.size()-n;ipos++);
	for(int i=ipos;i<s.size();i++) ret=ret+s[i];
	return ret;
}

/*	The callback function for -help option.

 * */
void	help_cb(string s)
{
	if(s=="" || s=="-help")
	opt.printHelp();
	else
	{
		printf("========== HELP FOR OPTION %s ===============\n",s.c_str());
		s="-"+s;
		opt.printHelp(s);	
	}
	exit(EXIT_SUCCESS);
}

/*	The callback function for -group1 option.
 * */
void	group1_cb(string s)
{
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandDistance();

	if(s=="aromatic") ((CommandDistance*)command_list[st])->setTeam1(makeAromatic());
	else
	if(s=="aliphatic") ((CommandDistance*)command_list[st])->setTeam1(makeAliphatic());
	else
	if(s=="positive") ((CommandDistance*)command_list[st])->setTeam1(makePositive());
	else
	if(s=="negative") 
			((CommandDistance*)command_list[st])->setTeam1(makeNegative());
	else
	if(s=="hydrophobic") 
			((CommandDistance*)command_list[st])->setTeam1(makeHydrophobic());
	else              ((CommandDistance*)command_list[st])->setFile1(s);
	last_team="DISTANCE1";
}

/*	The callback function for -group2 option.
 * */
void	group2_cb(string s)
{
	int st = commandPos("DISTANCE");
	if(st == -1) psf_error("YOU MUST SPECIFY FIRST THE GROUP1 TEAM");
	if(s=="aromatic")  ((CommandDistance*)command_list[st])->setTeam2(makeAromatic());
	else
	if(s=="aliphatic") ((CommandDistance*)command_list[st])->setTeam2(makeAliphatic());
	else
	if(s=="positive")  ((CommandDistance*)command_list[st])->setTeam2(makePositive());
	else
	if(s=="negative")  ((CommandDistance*)command_list[st])->setTeam2(makeNegative());
	else
	if(s=="hydrophobic") 
			((CommandDistance*)command_list[st])->setTeam2(makeHydrophobic());
	else               ((CommandDistance*)command_list[st])->setFile2(s);
	last_team="DISTANCE2";
}

void	psf_cb(string s)
{
	psffile = s;
}

void	dcd_cb(string s)
{
	dcdfile = s;
}

Team	*setTeam(char *optarg,string command_name)
{
	vector<string> str;
	BreakList(optarg,str);
	int p=commandPos(command_name);
	command_list[p]->setSearch(str);
	return command_list[p]->getTeam();
}

void	search_cb(string s)
{
	char optarg[1024];
	strcpy(optarg,s.c_str());
	vector<string> str;
	LastTeam=NULL;
	
	if(s=="-seq" ||
		strlen(optarg)==0 || s.size() == 0 || s=="" ) 
			psf_error("You must provide an argument for the -seq command");

	if(backbone_flag==1 && donor_act==0 && acceptor_act==0)
	{
		BreakList(optarg,str);
		int p=commandPos("HBONDS");
		command_list[p]->setSearch(str);
		LastTeam = command_list[p]->getTeam();
	}
	else
	if(backbone_flag==1 && ilast==2)
	{
		BreakList(optarg,str);
		donorTeam=new Team(str);
		LastTeam = donorTeam;
	}
	else
	if(backbone_flag==1 &&  ilast==1)
	{
		BreakList(optarg,str);
		acceptorTeam=new Team(str);
		LastTeam=acceptorTeam;
	}
	else
	if(last_team[last_team.size()-1] == '1')
	{
		string st = last_team.substr(0,last_team.size()-1);
		BreakList(optarg,str);
		int p=commandPos(st);
		command_list[p]->setSearch1(str);
		LastTeam=command_list[p]->getTeam1();
	}
	else
	if(last_team[last_team.size()-1] == '2')
	{
		string st = last_team.substr(0,last_team.size()-1);
		BreakList(optarg,str);
		int p=commandPos(st);
		command_list[p]->setSearch2(str);
		LastTeam=command_list[p]->getTeam2();
	}
	else
	if(last_team!="")
	LastTeam =setTeam(optarg,last_team);

	last_team = "";
}

void	res_cb(string s)
{
	char optarg[1024];
	strcpy(optarg,s.c_str());
	vector<string> res_list;
	BreakList(optarg,res_list);
	int p=commandPos("DISTANCE");
	if(p!=-1)
	{
		((CommandDistance*)command_list[p])->setRes1(res_list);
		((CommandDistance*)command_list[p])->setRes2(res_list);
	}
	p=commandPos("RMSF");
	if(p!=-1) ((CommandRmsf*)command_list[p])->setRes(res_list);
	p=commandPos("TORSION");
	if(p!=-1) ((CommandAngle*)command_list[p])->setRes(res_list);
	if(backbone_flag==1 && donor_act==0 && acceptor_act==0)
	{
		p=commandPos("HBONDS");
		command_list[p]->setRes(res_list);
	}
	else
	if(backbone_flag==1 && ilast==2)
	{
		if(donorTeam) 	donorTeam->setResList(res_list);
	}
	else
	if(backbone_flag==1 &&  ilast==1)
	{
		if(acceptorTeam) acceptorTeam->setResList(res_list);
	}
}

void	critical_cb(string s)
{
	char optarg[1024];
	strcpy(optarg,s.c_str());
	BreakCritical(optarg);
}

void	first_cb(string s)
{
	first = atoi(s.c_str());
}

void	last_cb(string s)
{
	last = atoi(s.c_str());
}

void	skip_cb(string s)
{
	step = atoi(s.c_str());
}

void	angle_cb(string s)
{
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandAngle();
	vector<string> list;
	getCommaElements(s,list);
	((CommandAngle*)command_list[st])->setAngleList(list);
	last_team = "TORSION";
}

void	pdbtors_cb(string s)
{
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandPdbAngle();
	vector<string> list;
	getCommaElements(s,list);
	((CommandPdbAngle*)command_list[st])->setAngleList(list);
	last_team = "PDBTORSION";
}
void	backbone_cb(string s)
{
	last_team="";
	angle_flag=0;
	isbackbone=1;
	backbone_flag = 1;
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandBackbone();
}

void	smart_cb(string s)
{
	if(s=="") return;
	string ss1="";
	string ss2="";
	int index=0;
	ss1=getNextWord(s,index,',');
	ss2=getNextWord(s,index,',');
	if(ss1!="") smart_skip=atoi(ss1.c_str());
	if(ss2!="") smart_distance=atof(ss2.c_str());
}

void	donor_cb(string s)
{
	if(s=="") return;
	last_team="";
	angle_flag=0;
	isbackbone=1;
	backbone_flag = 1;
	ilast=2;
	donor_act=1;
	donor_selection=s;
	if(donor_or_acceptor_first==0)
		donor_or_acceptor_first=1;
}

void	acceptor_cb(string s)
{
	if(s=="") return;
	last_team="";
	angle_flag=0;
	isbackbone=1;
	backbone_flag = 1;
	acceptor_act=1;
	acceptor_selection=s;
	ilast=1;
	if(donor_or_acceptor_first==0)
		donor_or_acceptor_first=2;
}

void	rmsf_cb(string s)
{
	last_team = "RMSF";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandRmsf();
}

void	bturn_cb(string s)
{
	int st=command_list.size();
	command_list.resize(st+1);

	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data bturn_percent = 0.1;
	Data bturn_distance = 7.0;
	Data bturn_angle = 60.0;
	if(s1.size() && isdigit(s1[0]))  bturn_percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  bturn_distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  bturn_angle   =atof(s3.c_str());

	command_list[st]=new CommandBturn(bturn_percent,bturn_distance,bturn_angle);
	last_team = "BTURN";
}

void	unary_distance_cb(string s)
{
	string s1="",s2="";
	int index=0;
	if(s.size()==0) psf_error("You muste specify arguments for distance command");

	string list1,list2;
	list1=getNextWord(s,index,'-');
	list2=getNextWord(s,index,'-');
	vector<string> comma_list1;
	vector<string> comma_list2;
	vector<string> t;
	getCommaElements(list1,comma_list1);
	for(int i=0;i<comma_list1.size();i++) 
		t.push_back(comma_list1[i]);

	string seperator="-";

	t.push_back(seperator);

	getCommaElements(list2,comma_list2);
	for(int i=0;i<comma_list2.size();i++) 
		t.push_back(comma_list2[i]);

	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandUnary(UNARY_DISTANCE,t);
}

void	unary_angle_cb(string s)
{
	string s1="",s2="",s3="";
	if(s.size()==0) psf_error("You muste specify arguments for angle command");
	string list1,list2,list3;
	int index=0;
	list1=getNextWord(s,index,'-');
	list2=getNextWord(s,index,'-');
	list3=getNextWord(s,index,'-');
	vector<string> comma_list1;
	vector<string> comma_list2;
	vector<string> comma_list3;
	vector<string> t;
	getCommaElements(list1,comma_list1);
	for(int i=0;i<comma_list1.size();i++) 
		t.push_back(comma_list1[i]);

	string seperator="-";

	t.push_back(seperator);

	getCommaElements(list2,comma_list2);
	for(int i=0;i<comma_list2.size();i++) 
		t.push_back(comma_list2[i]);

	t.push_back(seperator);

	getCommaElements(list3,comma_list3);
	for(int i=0;i<comma_list3.size();i++) 
		t.push_back(comma_list3[i]);

	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandUnary(UNARY_ANGLE,t);
}

void 	unary_dihedral_cb(string s)
{
	string s1="",s2="",s3="",s4="";
	if(s.size()==0) psf_error("You muste specify arguments for dihedral command");
	int index=0;
	string list1,list2,list3,list4;
	list1=getNextWord(s,index,'-');
	list2=getNextWord(s,index,'-');
	list3=getNextWord(s,index,'-');
	list4=getNextWord(s,index,'-');
	vector<string> comma_list1;
	vector<string> comma_list2;
	vector<string> comma_list3;
	vector<string> comma_list4;
	vector<string> t;
	getCommaElements(list1,comma_list1);
	for(int i=0;i<comma_list1.size();i++) 
		t.push_back(comma_list1[i]);

	string seperator="-";

	t.push_back(seperator);

	getCommaElements(list2,comma_list2);
	for(int i=0;i<comma_list2.size();i++) 
		t.push_back(comma_list2[i]);

	t.push_back(seperator);

	getCommaElements(list3,comma_list3);
	for(int i=0;i<comma_list3.size();i++) 
		t.push_back(comma_list3[i]);

	t.push_back(seperator);

	getCommaElements(list4,comma_list4);
	for(int i=0;i<comma_list4.size();i++) 
		t.push_back(comma_list4[i]);


	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandUnary(UNARY_DIHEDRAL,t);
}

void	stack_cb(string s)
{
	last_team  = "STACK";
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data stack_percent = 0.1;
	Data stack_distance = 5.5;
	Data stack_angle = 30.0;
	if(s1.size() && isdigit(s1[0]))  stack_percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  stack_distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  stack_angle   =atof(s3.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandStack(stack_percent,stack_distance,stack_angle);
}

void	watbridge1_cb(string s)
{
	watbridge_flag = 1;
	isbackbone=1;
	last_team = "WATBRIDGE1";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandWatbridge();
}

void	watbridge2_cb(string s)
{
	last_team = "WATBRIDGE2";
	isbackbone=1;
	int p=commandPos("WATBRIDGE");
	if(p==-1) psf_error("You must use first the option -watbridge1");
}

void	psfanal_cb(string s)
{
	last_team ="PSFANAL";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandAnal();
}

void	salt_cb(string s)
{
	last_team = "SALT";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandSalt();
}

void	complexsalt_cb(string s)
{
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandComplexSalt();
	last_team = "COMPLEXSALT";
}

void	noe_cb(string s)
{
	last_team ="NOE";
	isnoe = 1;
	vector<string> noe_atom;
	noe_atom.resize(0);
	if(s.size()==0 || s=="-noe") 
	{
		noe_atom.resize(2);
		noe_atom[0]="HA";
		noe_atom[1]="HN";
	}
	else
	{
		getCommaElements(s,noe_atom);
	}
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandNoe();
	((CommandNoe*)command_list[st])->setAtoms(noe_atom);
}

void	hist_cb(string s)
{
	int k=atoi(s.c_str());
	if(k) histflag=1;
	else  histflag=0;
}

void	side1_cb(string s)
{
	last_team = "SIDEDIST1";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandSideDist();
	((CommandSideDist*)command_list[st])->setName1(s);
}

void	side2_cb(string s)
{
	int p=commandPos("SIDEDIST");
	if(p==-1)
		psf_error("You must specify the first side team in order to use this facility");
	((CommandSideDist*)command_list[p])->setName2(s);
	last_team = "SIDEDIST2";
}


void	hbonds1_cb(string s)
{
	last_team="IHBONDS1";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandHbonds1();
}


void	hbonds2_cb(string s)
{
	last_team = "IHBONDS2";
	int p=commandPos("IHBONDS");
	if(p==-1)
	{
		int st=command_list.size();
		command_list.resize(st+1);
		command_list[st]=new CommandHbonds1();
	}
}

void	center1_cb(string s)
{
	last_team = "CENTER1";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandCenter();
	vector<string> flags;
	if(s=="" || s[0]=='-') flags.resize(0);
	else
	getCommaElements(s,flags);
	((CommandCenter*)command_list[st])->setFlags1(flags);
}
	

void	center2_cb(string s)
{
	last_team = "CENTER2";
	int p=commandPos("CENTER");
	if(p==-1)
	{
		int st=command_list.size();
		command_list.resize(st+1);
		command_list[st]=new CommandCenter();
		p=st;
	}
	vector<string> flags;
	if(s=="" || s[0]=='-') flags.resize(0);
	else
	getCommaElements(s,flags);
	((CommandCenter*)command_list[p])->setFlags2(flags);
}


void	bindist_cb(string s)
{
	bindist = atof(s.c_str());
}

void	binangle_cb(string s)
{
	binangle = atof(s.c_str());
}

void	bindihe_cb(string s)
{
	bindihe = atof(s.c_str()); 
}

void	aronh_cb(string s)
{
	last_team  = "ARONH";
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data stack_percent = 0.1;
	Data stack_distance = 5.5;
	Data stack_angle = 30.0;
	if(s1.size() && isdigit(s1[0]))  stack_percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  stack_distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  stack_angle   =atof(s3.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandAronh(stack_percent,stack_distance,stack_angle);
}

void	aroch_cb(string s)
{
	last_team  = "AROCH";
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data stack_percent = 0.1;
	Data stack_distance = 5.5;
	Data stack_angle = 30.0;
	if(s1.size() && isdigit(s1[0]))  stack_percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  stack_distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  stack_angle   =atof(s3.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandAroch(stack_percent,stack_distance,stack_angle);
}

void	aropos_cb(string s)
{
	last_team  = "AROPOS";
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data stack_percent = 0.1;
	Data stack_distance = 5.5;
	Data stack_angle = 30.0;
	if(s1.size() && isdigit(s1[0]))  stack_percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  stack_distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  stack_angle   =atof(s3.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandAropos(stack_percent,stack_distance,stack_angle);
}

void	rmsd_cb(string s)
{
	last_team = "RMSD";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandRmsd();
}

void	avenoe_cb(string s)
{
	Data avg = atof(s.c_str());
	int pos=commandPos("NOE");
	if(pos!=-1)
	{
		((CommandNoe*)command_list[pos])->setAverage(avg);
	}
}

void	jhnha_cb(string s)
{
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data a=6.4;
	Data b=-1.4;
	Data c=1.9;
	if(s1.size() && isdigit(s1[0]))  a =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  b=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  c=atof(s3.c_str());
	int pos=commandPos("TORSION");
	if(pos!=-1)
	{
		((CommandAngle*)command_list[pos])->setJHNHA(a,b,c);
	}
}

void	pdb_cb(string s)
{
	pdbfile=s;
}

void	smooth_cb(string s)
{
	smooth_flag = 1;
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	if(s1.size() && isdigit(s1[0]))  smooth_start = atoi(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  smooth_step   = atoi(s2.c_str());
}

void	pdbwrite_cb(string s)
{
	last_team = "PDBWRITE";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandPdbWrite();
}

void	weight_cb(string s)
{
	weight_flag = 1;
}

void	pdo_cb(string s)
{
	last_team = "PDO";
	int level=2;
	if(s.size() && isdigit(s[0])) level=atoi(s.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandPdo(level);
}


void	contact1_cb(string s)
{
	last_team = "CONTACT1";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandContact();
	if(s[0]!='-') ((CommandContact*)command_list[st])->setSelection1(s);
	contact_flag=1;
}

void	contact2_cb(string s)
{
	last_team = "CONTACT2";
	int p=commandPos("CONTACT");
	if(p==-1) psf_error("You must use  the -contact1 option  before the -contact2 option");
	if(s[0]!='-') ((CommandContact*)command_list[p])->setSelection2(s);
	contact_flag=1;
}

void	resdiff_cb(string s)
{
	res_diff=atoi(s.c_str());
}

void	resname_cb(string s)
{
	vector<string> str;
	getCommaElements(s,str);
	if(LastTeam!=NULL)
		LastTeam->setResName(str);
}

void	atomname_cb(string s)
{
	vector<string> str;
	getCommaElements(s,str);
	if(LastTeam!=NULL)
		LastTeam->setAtomName(str);
}


void	atomtype_cb(string s)
{
	vector<string> str;
	getCommaElements(s,str);
	if(LastTeam!=NULL)
		LastTeam->setAtomType(str);
}

void	hydrosolve_cb(string s)
{
	last_team = "HYDROSOLVE";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandHydroSolve();
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	Data percent = 0.05;
	Data distance = 3.0;
	Data angle = 145.0;
	if(s1.size() && isdigit(s1[0]))  percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  angle   =atof(s3.c_str());
	((CommandHydroSolve*)command_list[st])->setPda(percent,distance,angle);
}

void	closewater_cb(string s)
{
	last_team = "HPC";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandCloseWater();
}


void	acor_cb(string s)
{
	int index=0;
	string s1=getNextWord(s,index,',');
	if(s1.size() && isdigit(s1[0])) kstep=atoi(s1.c_str());
	string s2=getNextWord(s,index,',');
	if(s2.size() && isdigit(s2[0])) astep=atoi(s2.c_str());
}

void	stat_cb(string s)
{
	statcount=atoi(s.c_str());
}

void	iwcn_cb(string s)
{
	double distance=3.5,percent=0.05;
	int maxn=1;
	int index=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	string s3=getNextWord(s,index,',');
	if(s1.size() && isdigit(s1[0]))  percent =atof(s1.c_str());
	if(s2.size() && isdigit(s2[0]))  distance=atof(s2.c_str());
	if(s3.size() && isdigit(s3[0]))  maxn = atoi(s3.c_str());

	last_team = "IWCN";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandIwcn(percent,distance,maxn);
}

void	mol_cb(string s)
{
	string pfile=s+".psf";
	string dfile=s+".dcd";
	string pdfile=s+".pdb";
	psf_cb(pfile);
	dcd_cb(dfile);
	pdb_cb(pdfile);
}

void	ahelix_cb(string s)
{
	double def_percent=0.5;
	if(s.size() && isdigit(s[0]))  def_percent = atof(s.c_str());
	last_team = "AHELIX";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandHelix(def_percent);
}

void	watpdb_cb(string s)
{
	double	def_distance=3.5;
	if(s.size() && isdigit(s[0])) def_distance = atof(s.c_str());
	last_team = "WATPDB";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandWatpdb(def_distance);
}

void	clusterrmsd_cb(string s)
{
	last_team = "CLUSTERRMSD";
	double l=0.3;
	if(s.size() && isdigit(s[0])) l=atof(s.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandClusterRmsd(l);
}

void	entropy_cb(string s)
{
	last_team="ENTROPY";
	int n=10;
	Data temp=310.0;
	int index=0;
	int flag=0;
	int fitflag=1;
	int pairflag=0;
	string s1=getNextWord(s,index,',');
	if(s1.size() && isdigit(s1[0])) n=atoi(s1.c_str());
	string s2=getNextWord(s,index,',');
	if(s2.size() && isdigit(s2[0])) temp=atof(s2.c_str());
	string s3=getNextWord(s,index,',');
	if(s3=="res") flag=1;
	if(s3=="nofit") fitflag=0;
	if(s3=="fit")   fitflag=1;
	if(s3=="bbfit") fitflag=2;
	if(s3=="cafit") fitflag=3;
	if(s3=="pair")  pairflag=1;
	string s4=getNextWord(s,index,',');
	if(s4=="res") flag=1;
	if(s4=="nofit") fitflag=0;
	if(s4=="fit")   fitflag=1;
	if(s4=="bbfit") fitflag=2;
	if(s4=="cafit") fitflag=3;
	if(s4=="pair")  pairflag=1;
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandEntropy(n,temp,flag,fitflag,pairflag);
}

void	writedcd_cb(string s)
{
	last_team="WRITEDCD";
	int index=0;
	int fit=0;
	string s1=getNextWord(s,index,',');
	string s2=getNextWord(s,index,',');
	if(s2=="fit")   fit=1;
	if(s2=="bbfit") fit=2;
	if(s2=="cafit") fit=3;
	string s3=s1.size()>=4?s1.substr(s1.size()-4,4):"";
	if(s3!=".dcd") s1=s1+".dcd";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandWriteDcd(s1,fit);
}

void	par_cb(string s)
{
	parfile = s;
}

void	energy_cb(string s)
{
	last_team="ENERGY";
	vector<string> str;
	getCommaElements(s,str);
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandEnergy(str);
}

void	elect_cb(string s)
{
	if(s.size() && isdigit(s[0])) diel=atof(s.c_str());
}

void	disu_cb(string s)
{
	last_team="DISU";
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandDisu();
}

void	version_cb(string s)
{
	printVersion();
}

void	sasa_cb(string s)
{
	last_team="SASA";
	Data p=1.4;
	if(isdigit(s[0]))  p=atof(s.c_str());
	int st=command_list.size();
	command_list.resize(st+1);
	command_list[st]=new CommandSasa(p);
}

void	ref_cb(string s)
{
	if(!isdigit(s[0])) psf_error("You must provide a reference number for -ref option");
	ref=atoi(s.c_str());
}

void parse_cmd_line(int argc,char **argv)
{
	parse_help();
	char cmdline[2048];
	strcpy(cmdline,"");
	for(int i=1;i<argc;i++)
		sprintf(cmdline,"%s %s",cmdline,argv[i]);

	opt.addOption("-help",help_cb,helpDescription("-help"));
	opt.addOption("-psf",psf_cb,helpDescription("-psf"));
	opt.addOption("-dcd",dcd_cb,helpDescription("-dcd"));

	opt.addOption("-acceptor",acceptor_cb,helpDescription("-acceptor"));
	opt.addOption("-ahelix",ahelix_cb,helpDescription("-ahelix"));
	opt.addOption("-block",stat_cb,helpDescription("-block"));
	opt.addOption("-acor",acor_cb,helpDescription("-acor"));
	opt.addOption("-angle",unary_angle_cb,helpDescription("-angle"));
	opt.addOption("-aroHN",aronh_cb,helpDescription("-aroHN"));
	opt.addOption("-aroHC",aroch_cb,helpDescription("-aroHC"));
	opt.addOption("-aropos",aropos_cb,helpDescription("-aropos"));
	opt.addOption("-atomname",atomname_cb,helpDescription("-atomname"));
	opt.addOption("-atomtype",atomtype_cb,helpDescription("-atomname"));
	opt.addOption("-avenoe",avenoe_cb,helpDescription("-avenoe"));

	opt.addOption("-bindist",bindist_cb,helpDescription("-bindist"));
	opt.addOption("-binangle",binangle_cb,helpDescription("-binangle"));
	opt.addOption("-bindihe",bindihe_cb,helpDescription("-bindihe"));
	opt.addOption("-bturn",bturn_cb,helpDescription("-bturn"));


	opt.addOption("-center1",center1_cb,helpDescription("-center1"));
	opt.addOption("-center2",center2_cb,helpDescription("-center2"));
	opt.addOption("-hpc",closewater_cb,helpDescription("-closewater"));
	opt.addOption("-clusterrmsd",clusterrmsd_cb,helpDescription("-clusterrmsd"));
	opt.addOption("-contact1",contact1_cb,helpDescription("-contact1"));
	opt.addOption("-contact2",contact2_cb,helpDescription("-contact2"));
	opt.addOption("-cutoff",critical_cb,helpDescription("-cutoff"));

	opt.addOption("-donor",donor_cb,helpDescription("-donor"));
	opt.addOption("-distance",unary_distance_cb,helpDescription("-distance"));
	opt.addOption("-dihedral",unary_dihedral_cb,helpDescription("-dihedral"));
	opt.addOption("-disu",disu_cb,helpDescription("-disu"));
	opt.addOption("-elect",elect_cb,helpDescription("-elect"));
	opt.addOption("-energy",energy_cb,helpDescription("-energy"));
	opt.addOption("-entropy",entropy_cb,helpDescription("-entropy"));

	opt.addOption("-first",first_cb,helpDescription("-first"));

	opt.addOption("-group1",group1_cb,helpDescription("-group1"));
	opt.addOption("-group2",group2_cb,helpDescription("-group2"));

	opt.addOption("-hbonds",backbone_cb,helpDescription("-hbonds"));
	opt.addOption("-hbonds1",hbonds1_cb,helpDescription("-hbonds1"));
	opt.addOption("-hbonds2",hbonds2_cb,helpDescription("-hbonds2"));

	opt.addOption("-hist",hist_cb,helpDescription("-hist"));
	opt.addOption("-hydrosolve",hydrosolve_cb,helpDescription("-hydrosolve"));

	opt.addOption("-iwcn",iwcn_cb,helpDescription("-iwcn"));

	opt.addOption("-JHNHA",jhnha_cb,helpDescription("-JHNHA"));

	opt.addOption("-last",last_cb,helpDescription("-last"));

	opt.addOption("-mol",mol_cb,helpDescription("-mol"));

	opt.addOption("-noe",noe_cb,helpDescription("-noe"));

	opt.addOption("-par",par_cb,helpDescription("-par"));
	opt.addOption("-pdb",pdb_cb,helpDescription("-pdb"));
	opt.addOption("-pdbtors",pdbtors_cb,helpDescription("-pdbtors"));
	opt.addOption("-pdbwrite",pdbwrite_cb,helpDescription("-pdbwrite"));
	opt.addOption("-pdo",pdo_cb,helpDescription("-pdo"));
	opt.addOption("-psfanal",psfanal_cb,helpDescription("-psfanal"));

	//opt.addOption("-res",res_cb,helpDescription("-res"));
	opt.addOption("-ref",ref_cb,helpDescription("-ref"));
	opt.addOption("-resdiff",resdiff_cb,helpDescription("-resdiff"));
	opt.addOption("-resname",resname_cb,helpDescription("-resname"));
	opt.addOption("-rmsd",rmsd_cb,helpDescription("-rmsd"));
	opt.addOption("-rmsf",rmsf_cb,helpDescription("-rmsf"));

	opt.addOption("-salt2",salt_cb,helpDescription("-salt2"));
	opt.addOption("-salt3",complexsalt_cb,helpDescription("-salt3"));
	opt.addOption("-sasa",sasa_cb,helpDescription("-sasa"));
	opt.addOption("-seq",search_cb,helpDescription("-seq"));
	opt.addOption("-side1",side1_cb,helpDescription("-side1"));
	opt.addOption("-side2",side2_cb,helpDescription("-side2"));
	opt.addOption("-skip",skip_cb,helpDescription("-skip"));
	opt.addOption("-smart",smart_cb,helpDescription("-smart"));
	opt.addOption("-smooth",smooth_cb,helpDescription("-smooth"));
	opt.addOption("-stack",stack_cb,helpDescription("-stack"));

	opt.addOption("-tors",angle_cb,helpDescription("-tors"));
	opt.addOption("-watbridge1",watbridge1_cb,helpDescription("-watbridge1"));
	opt.addOption("-watbridge2",watbridge2_cb,helpDescription("-watbridge2"));
	opt.addOption("-watpdb",watpdb_cb,helpDescription("-watpdb"));
	opt.addOption("-weight",weight_cb,helpDescription("-weight"));
	opt.addOption("-dcdwrite",writedcd_cb,helpDescription("-dcdwrite"));
	opt.addOption("-v",version_cb,helpDescription("-v"));
	

	opt.parse(cmdline);
	if(argc==1) {help_cb("");}

	
	while(psffile == "")	 
	{
		printf("PLEASE PROVIDE THE PSF FILE: ");
		cin>>psffile;
	}
	while(dcdfile == "")
	{
		printf("PLEASE PROVIDE THE DCD FILE: ");
		cin>>dcdfile;
	}
	if(donor_act)
	{
		int st=command_list.size();
		command_list.resize(st+1);
		command_list[st]=new CommandDonor();
		if(donor_or_acceptor_first==1) 
		{
			((CommandDonor*)command_list[st])->setFlag(1);
		}
		else 
		{
			((CommandDonor*)command_list[st])->setFlag(2);
		}
		((CommandDonor*)command_list[st])->setTeams(donorTeam,acceptorTeam);
	}
}
