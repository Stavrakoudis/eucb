# include <command_distance.h>		
# include <global.h>
# include <unistd.h>
# include <tex.h>
# include <counter.h>


CommandDistance::CommandDistance()
	:Command("DISTANCE")
{
}

void	CommandDistance::Run()
{
	Counter myCounter;
	string suffix="";
	int icount=0;
	string chain1="",chain2="";
	vector<int> to_team1;
	vector<int> to_team2;
	vector<NoeStruct> noe;

	Tex myTex("distance.tex");
	vector<string> s;
	s.resize(2);
	s[0]="Atom1";
	s[1]="Atom2";
	myTex.setHeader(s);
	myTex.setCaption("Atom distance");

	for(int i=0;i<table.size();i++)
	{
		if(isHydro(table[i])) continue;
		if(isBackbone(table[i])) continue;

		if(team1->find(table[i],chain1))
		{
			to_team1.push_back(i);
		}
		if(team2 && team2->find(table[i],chain2)) 
		{
			to_team2.push_back(i);
		}
	}
	
	NoeStruct n;
	if(!team2)
	{
		for(int i=0;i<to_team1.size();i++)
		{
			for(int j=0;j<i;j++)
			{
				if(table[to_team1[i]].res_id==table[to_team1[j]].res_id) continue;
				n.atom1=to_team1[i]+1;
				n.atom2=to_team1[j]+1;
				noe.push_back(n);
			}
		}
	}
	else
	{
		for(int i=0;i<to_team1.size();i++)
		{
			for(int j=0;j<to_team2.size();j++)
			{
				if(to_team1[i]==to_team2[j]) continue;
				if(table[to_team1[i]].res_id==table[to_team2[j]].res_id) continue;
				n.atom1=to_team1[i]+1;
				n.atom2=to_team2[j]+1;
				noe.push_back(n);
			}
		}
	}
	IntVector f;
	dcd->getNoe(f,noe);
	Log("#Critical distance = "+printDistance(critical_distance));
	Log("#Critical percent  = "+printPercent(critical_percent));
	for(int i=0;i<noe.size();i++)
	{
		Data p=getBelowPercent(noe[i].D,critical_distance);
		IntVector Acor;
		Acor.resize(f.size());
		for(int j=0;j<f.size();j++)
		{
			if(noe[i].D[j]<=critical_distance) Acor[j]=1; else Acor[j]=0;
			myCounter.add(f[j],"",noe[i].D[j]<=critical_distance?1:0);
		}

		if(p>=critical_percent)
		{
			string filename="dist_"+printAtom1(table[noe[i].atom1-1])+"_"+
				printAtom1(table[noe[i].atom2-1])+".dat";
			string smoothname="dist_"+printAtom1(table[noe[i].atom1-1])+"_"+
				printAtom1(table[noe[i].atom2-1])+".sda";
			string histname="dist_"+printAtom1(table[noe[i].atom1-1])+"_"+
				printAtom1(table[noe[i].atom2-1])+".hist";
			
			string acorname="dist_"+printAtom1(table[noe[i].atom1-1])+"_"+
				printAtom1(table[noe[i].atom2-1])+".cor";
			PrintAcor(acorname,f,Acor,kstep,astep);

			IntVector atom_pos;
			atom_pos.resize(2);
			atom_pos[0]=noe[i].atom1;
			atom_pos[1]=noe[i].atom2;
			myTex.addAtomList(atom_pos,p);

			PrintDistances(filename,f,noe[i].D);
			Log(printAtom3(table[noe[i].atom1-1])+" "+printAtom3(table[noe[i].atom2-1])+
				" "+printPercent(p));
			if(smooth_flag)
			{
				IntVector sf;
				DoubleVector sd;
				makeSmooth(f,noe[i].D,sf,sd,smooth_start,smooth_step);
				PrintDistances(smoothname,sf,sd);
			}
			if(histflag)
			{
				Data dmin,dmax,davg,dstd;
				vector<HistStruct> st;
				getVectorStatistics(noe[i].D,dmin,dmax,davg,dstd);
				makeHist(noe[i].D,st,dmin,dmax,bindist);
				PrintHist(histname,st);
			}
		}
	}
	myTex.print();
	myCounter.print("group_count.dat");
	myCounter.printStat("group_count.stat");
	if(smooth_flag) myCounter.printSmooth("group_count.sda");
	if(histflag)    myCounter.printHist("group_count.hist");
}

void  CommandDistance::setFile1(string s)
{
	vector<string> st;
	char file[1024];
	strcpy(file,s.c_str());
	team1 = new Team(file,st);
}

void  CommandDistance::setFile2(string s)
{
	vector<string> st;
	char file[1024];
	strcpy(file,s.c_str());
	team2 = new Team(file,st);
}


void  CommandDistance::setTeam1(Team *t)
{
	team1 = t;
}

void  CommandDistance::setTeam2(Team *t)
{
	team2 = t;
}

void  CommandDistance::setRes1(vector<string> &s)
{
	if(team1!=NULL) team1->setResList(s);
}

void  CommandDistance::setRes2(vector<string> &s)
{
	if(team2!=NULL) team1->setResList(s);
}

CommandDistance::~CommandDistance()
{
}
