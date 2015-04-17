# include <command_center.h>		
# include <global.h>
# include <unistd.h>


CommandCenter::CommandCenter()
	:Command("CENTER")
{
	ca_flag1=0;
	noh_flag1=1;
	all_flag1=0;
	ca_flag2=0;
	noh_flag2=1;
	sidechain_flag1=0;
	sidechain_flag2=0;
	all_flag2=0;
}

void	CommandCenter::setFlags1(vector<string> &flags1)
{
	for(int i=0;i<flags1.size();i++) 
	{
		if(flags1[i]=="ca") ca_flag1=1;
		if(flags1[i]=="noh") noh_flag1=1;
		if(flags1[i]=="all") all_flag1=1;
		if(flags1[i]=="sidechain") sidechain_flag1=1;
	}	
}

void	CommandCenter::setFlags2(vector<string> &flags2)
{
	for(int i=0;i<flags2.size();i++) 
	{
		if(flags2[i]=="ca") ca_flag2=1;
		if(flags2[i]=="noh") noh_flag2=1;
		if(flags2[i]=="all") all_flag2=1;
		if(flags2[i]=="sidechain") sidechain_flag2=1;
	}	
}


void	CommandCenter::Run()
{
	string suffix="";
	int icount=0;
	string chain1="",chain2="";
	IntVector to_team1;
	IntVector to_team2;
	for(int i=0;i<table.size();i++)
	{
		int flag1=0;
		int flag2=0;
		if(!noh_flag1 && isHydro(table[i]) && (!team1 || team1->find(table[i],chain1))) flag1=1;
		if(ca_flag1 && table[i].atom_name=="CA" && (!team1 || team1->find(table[i],chain1))) flag1=1;
		if(sidechain_flag1 && isSidechain(table[i]) && (!team1 || team1->find(table[i],chain1))) flag1=1;
		if(all_flag1 && (!team1 || team1->find(table[i],chain1))) flag1=1;
		if(!noh_flag2 && isHydro(table[i]) && (!team2 || team2->find(table[i],chain2))) flag2=1;
		if(ca_flag2 && table[i].atom_name=="CA" && (!team2 || team2->find(table[i],chain2))) flag2=1;
		if(sidechain_flag2 && isSidechain(table[i]) && (!team2 || team2->find(table[i],chain2))) flag2=1;
		if(all_flag2 && (!team2 || team2->find(table[i],chain2))) flag2=1;
		if(flag1 ) to_team1.push_back(table[i].atom_id);
		if(flag2) to_team2.push_back(table[i].atom_id);

	}
	
	
	IntVector f;
	DoubleVector Distance;
	Log("Center 1 size = "+printNumber((int)to_team1.size()));
	Log("Center 2 size = "+printNumber((int)to_team2.size()));
	dcd->getNoeCenter(to_team1,to_team2,f,Distance);
	string filename="dcenter.dat";
	string smoothname="dcenter.sda";
	string histname="dcenter.hist";
	string statname="dcenter.stat";
	PrintDistances(filename,f,Distance);
	PrintStatFile(statname,"w",f.size()/statcount,"#Distance statistics",f,Distance);
	if(smooth_flag)
	{
		IntVector sf;
		DoubleVector sd;
		makeSmooth(f,Distance,sf,sd,smooth_start,smooth_step);
		PrintDistances(smoothname,sf,sd);
	}
	if(histflag)
	{
		vector<HistStruct> st;
		Data dmin,dmax,davg,dstd;
		getVectorStatistics(Distance,dmin,dmax,davg,dstd);
		makeHist(Distance,st,dmin,dmax,bindist);
		PrintHist(histname,st);
	}
}


CommandCenter::~CommandCenter()
{
}
