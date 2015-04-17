# include <global.h>
# include <math.h>
# include <command_unary.h>

CommandUnary::CommandUnary(int Flag,vector<string> &t)
	:Command("UNARY")
{
	flag = Flag;
	unary_atom.resize(t.size());
	for(int i=0;i<t.size();i++) unary_atom[i]=t[i];
}

int CommandUnary::discover_atom(string s)
{
	string chainid="";
	string resid="";
	string atomname="";
	int index=0;
	chainid=getNextWord(s,index,':');
	resid=getNextWord(s,index,':');
	atomname=getNextWord(s,index,':');
	int iresid=atoi(resid.c_str());
	for(int i=0;i<table.size();i++)
	{
		if(table[i].chain_id==chainid && table[i].res_id==iresid && table[i].atom_name==atomname)
			return table[i].atom_id;
	}
	return -1;
}

void	CommandUnary::unary_print(string filename,string statname,string histname,string smoothname,IntVector f,
		DoubleVector x,int is_angle)
{
	if(is_angle) PrintAngles(filename,f,x); 
        else         PrintDistances(filename,f,x);
	if(smooth_flag)
	{
		IntVector sf;
		DoubleVector sd;
		if(is_angle)
			makeAngleSmooth(f,x,sf,sd,smooth_start,smooth_step);
		else
			makeSmooth(f,x,sf,sd,smooth_start,smooth_step);
		if(is_angle) PrintAngles(smoothname,sf,sd);
		else         PrintDistances(smoothname,sf,sd);
	}
	Data delta=bindist;
	if(is_angle==1) delta=binangle;
	if(is_angle==2) delta=bindihe;
	Data dmin,dmax,davg,dstd;
	if(is_angle)
		getAngleVectorStatistics(x,dmin,dmax,davg,dstd);
	else
		getVectorStatistics(x,dmin,dmax,davg,dstd);
	if(is_angle==2) {dmin=-180.0;dmax=180.0;}
	if(!is_angle)
		PrintStatFile(statname,"w",f.size()/statcount,"#Statistics for distance command",f,x);
	else	PrintAngleStatFile(statname,"w",f.size()/statcount,"#Statistics for angle command",f,x);
	if(histflag)
	{
		vector<HistStruct> st;
		makeHist(x,st,dmin,dmax,delta);
		if(is_angle) PrintAngleHist(histname,st);
		else         PrintHist(histname,st);
	}
}

/*	Print the distance between two particular atoms.
 * */
void CommandUnary::unary_distance()
{
	IntVector frame;
	DoubleVector distance;
	IntVector Atom1;
	IntVector Atom2;
	int i;
	Log("Creating distances ");
	Log("First group of atoms ");
	for(i=0;unary_atom[i]!="-";i++)
	{
		int d=discover_atom(unary_atom[i]);
		Log(printAtom3(table[d-1]));
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Atom1.push_back(d);
	}
	i++;
	Log("Second group of  atoms ");
	for(;i<unary_atom.size();i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom2.push_back(d);
	}
	dcd->getNoeCenter(Atom1,Atom2,frame,distance);
	string filename,statname,histname,smoothname;
	if(unary_atom.size()==3)
	{
		int atom1,atom2;
		atom1=Atom1[0];

		atom2=Atom2[0];
		filename="dist_"+printAtom1(table[atom1-1])+"_"+
			printAtom1(table[atom2-1])+".dat";
		statname="dist_"+printAtom1(table[atom1-1])+"_"+
			printAtom1(table[atom2-1])+".stat";
		histname="dist_"+printAtom1(table[atom1-1])+"_"+
			printAtom1(table[atom2-1])+".hist";
		smoothname="dist_"+printAtom1(table[atom1-1])+"_"+
			printAtom1(table[atom2-1])+".sda";
	}
	else
	{
		filename="dist.dat";
		statname="dist.stat";
		histname="dist.hist";
		smoothname="dist.sda";
	}
	unary_print(filename,statname,histname,smoothname,frame,distance,0);
}

/*	Print the angle between two particular atoms.
 * */
void CommandUnary::unary_angle()
{
	IntVector frame;
	DoubleVector angle;
	IntVector Atom1;
	IntVector Atom2;
	IntVector Atom3;
	int i;
	Log("Creating angles");
	Log("First group of atoms ");
	for(i=0;unary_atom[i]!="-";i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom1.push_back(d);
	}
	i++;
	Log("Second group of atoms ");
	for(;unary_atom[i]!="-";i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom2.push_back(d);
	}
	i++;
	Log("Third group of atoms ");
	for(;i<unary_atom.size();i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom3.push_back(d);
	}
	dcd->getNoeCenterAngle(Atom1,Atom2,Atom3,frame,angle);	
	string filename,statname,histname,smoothname;
	if(unary_atom.size()==5)
	{
		int atom1,atom2,atom3;
		atom1=Atom1[0];

		atom2=Atom2[0];

		atom3=Atom3[0];
		filename="angle_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+".dat";
		statname="angle_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+".stat";
		histname="angle_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+".hist";
		smoothname="angle_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+".sda";
	}
	else
	{
		filename="angle.dat";
		statname="angle.stat";
		histname="angle.hist";
		smoothname="angle.sda";
	}
	unary_print(filename,statname,histname,smoothname,frame,angle,1);
}

/*	Print the dihedral angle between two particular atoms.
 * */
void CommandUnary::unary_dihedral()
{
	IntVector frame;
	DoubleVector angle;
	IntVector Atom1;
	IntVector Atom2;
	IntVector Atom3;
	IntVector Atom4;
	int i;
	Log("Creating dihedrals ");
	Log("First group of atoms ");
	for(i=0;unary_atom[i]!="-";i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom1.push_back(d);
	}
	i++;
	Log("Second group of atoms ");
	for(;unary_atom[i]!="-";i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom2.push_back(d);
	}
	i++;
	Log("Third group of atoms ");
	for(;unary_atom[i]!="-";i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom3.push_back(d);
	}
	i++;
	Log("Fourth group of atoms ");
	for(;i<unary_atom.size();i++)
	{
		int d=discover_atom(unary_atom[i]);
		if(d==-1) psf_error("Undefined atom "+unary_atom[i]);
		Log(printAtom3(table[d-1]));
		Atom4.push_back(d);
	}
	dcd->getNoeCenterDihedral(Atom1,Atom2,Atom3,Atom4,frame,angle);
	string filename,statname,histname,smoothname;
	if(unary_atom.size()==7)
	{
		int atom1,atom2,atom3,atom4;
		atom1=Atom1[0];

		atom2=Atom2[0];

		atom3=Atom3[0];

		atom4=Atom4[0];

		filename="dihe_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+"_"+
			printAtom1(table[atom4-1])+".dat";
		statname="dihe_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+"_"+
			printAtom1(table[atom4-1])+".stat";
		histname="dihe_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+"_"+
			printAtom1(table[atom4-1])+".hist";
		smoothname="dihe_"+printAtom1(table[atom1-1])+"_"+printAtom1(table[atom2-1])+"_"+
			printAtom1(table[atom3-1])+"_"+
			printAtom1(table[atom4-1])+".sda";
	}
	else
	{
		filename="dihe.dat";
		statname="dihe.stat";
		histname="dihe.hist";
		smoothname="dihe.sda";
	}
	unary_print(filename,statname,histname,smoothname,frame,angle,2);
}


void	CommandUnary::Run()
{
	if(flag == UNARY_DISTANCE) unary_distance();
	else
	if(flag == UNARY_ANGLE)    unary_angle();
	else                       unary_dihedral();
}

