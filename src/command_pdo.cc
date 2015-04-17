# include <command_pdo.h>
# include <math.h>
# include <global.h>

CommandPdo::CommandPdo(int l)
	:Command("PDO")
{
	if(l<=0) psf_error("Only positive values are allowed for pdo ");
	pdo_level = l;
}

typedef struct
{
	int atom1;
	int atom2;
	int atom3;
}PStruct;

void	CommandPdo::Run()
{
	Team *pdo_team = team;
	vector<PdoStruct> Pdo;
	vector<PStruct> ArgPos;
	vector<PStruct> LysPos;
	vector<PStruct> GluPos;
	vector<PStruct> AspPos;
	string chain;
	Log("Pdo level = "+printNumber(pdo_level));
	for(int i=0;i<table.size();i++)
	{
		if(!pdo_team || pdo_team->find(table[i],chain)) 
		{
			if(table[i].res_name == "ARG") 
			{
				PStruct p;
				while(i<table.size() && table[i].res_name=="ARG") 
				{
					if(table[i].atom_name=="CB") p.atom1=i+1;
					if(table[i].atom_name=="CA") p.atom2=i+1;
					if(table[i].atom_name=="CZ") p.atom3=i+1;
					i++;
				}
				ArgPos.push_back(p);
			}
			else
			if(table[i].res_name == "LYS")
			{
				PStruct p;
				while(i<table.size() && table[i].res_name=="LYS") 
				{
					if(table[i].atom_name=="CB") p.atom1=i+1;
					if(table[i].atom_name=="CA") p.atom2=i+1;
					if(table[i].atom_name=="NZ") p.atom3=i+1;
					i++;
				}
				LysPos.push_back(p);
			}
			else
			if(table[i].res_name == "GLU")
			{
				PStruct p;
				while(i<table.size() && table[i].res_name == "GLU") 
				{
					if(table[i].atom_name=="CB") p.atom1=i+1;
					if(table[i].atom_name=="CA") p.atom2=i+1;
					if(table[i].atom_name=="CD") p.atom3=i+1;
					i++;
				}
				GluPos.push_back(p);
			}
			else
			if(table[i].res_name == "ASP")
			{
				PStruct p;
				while(i<table.size() && table[i].res_name == "ASP") 
				{
					if(table[i].atom_name=="CB") p.atom1=i+1;
					if(table[i].atom_name=="CA") p.atom2=i+1;
					if(table[i].atom_name=="CG") p.atom3=i+1;
					i++;
				}
				AspPos.push_back(p);
			}
		}
	}
	for(int i=0;i<ArgPos.size();i++)
	{
		for(int j=0;j<GluPos.size();j++)
		{
			if(table[ArgPos[i].atom1-1].chain_id == table[GluPos[j].atom1-1].chain_id
			 && abs(table[ArgPos[i].atom1-1].res_id-table[GluPos[j].atom1-1].res_id)==pdo_level)
			{
				PdoStruct p;
				p.dist_atom1=ArgPos[i].atom1;
				p.dist_atom2=GluPos[j].atom1;
				p.dist_atom3=ArgPos[i].atom3;
				p.dist_atom4=GluPos[j].atom3;
				p.dihe_atom1=ArgPos[i].atom1;
				p.dihe_atom2=ArgPos[i].atom2;
				p.dihe_atom3=GluPos[j].atom2;
				p.dihe_atom4=GluPos[j].atom1;
				p.dihe_atom5=ArgPos[i].atom3;
				p.dihe_atom6=ArgPos[i].atom2;
				p.dihe_atom7=GluPos[j].atom2;
				p.dihe_atom8=GluPos[j].atom3;
				Pdo.push_back(p);
			}
	
		}
		for(int j=0;j<AspPos.size();j++)
		{
			if(table[ArgPos[i].atom1-1].chain_id == table[AspPos[j].atom1-1].chain_id
			 && abs(table[ArgPos[i].atom1-1].res_id-table[AspPos[j].atom1-1].res_id)==pdo_level)
			{
				PdoStruct p;
				p.dist_atom1=ArgPos[i].atom1;
				p.dist_atom2=AspPos[j].atom1;
				p.dist_atom3=ArgPos[i].atom3;
				p.dist_atom4=AspPos[j].atom3;
				p.dihe_atom1=ArgPos[i].atom1;
				p.dihe_atom2=ArgPos[i].atom2;
				p.dihe_atom3=AspPos[j].atom2;
				p.dihe_atom4=AspPos[j].atom1;
				p.dihe_atom5=ArgPos[i].atom3;
				p.dihe_atom6=ArgPos[i].atom2;
				p.dihe_atom7=AspPos[j].atom2;
				p.dihe_atom8=AspPos[j].atom3;
				Pdo.push_back(p);
			}
		}
	}

	for(int i=0;i<LysPos.size();i++)
	{
		for(int j=0;j<GluPos.size();j++)
		{
			if(table[LysPos[i].atom1-1].chain_id == table[GluPos[j].atom1-1].chain_id
			 && abs(table[LysPos[i].atom1-1].res_id-table[GluPos[j].atom1-1].res_id)==pdo_level)
			{
				PdoStruct p;
				p.dist_atom1=LysPos[i].atom1;
				p.dist_atom2=GluPos[j].atom1;
				p.dist_atom3=LysPos[i].atom3;
				p.dist_atom4=GluPos[j].atom3;
				p.dihe_atom1=LysPos[i].atom1;
				p.dihe_atom2=LysPos[i].atom2;
				p.dihe_atom3=GluPos[j].atom2;
				p.dihe_atom4=GluPos[j].atom1;
				p.dihe_atom5=LysPos[i].atom3;
				p.dihe_atom6=LysPos[i].atom2;
				p.dihe_atom7=GluPos[j].atom2;
				p.dihe_atom8=GluPos[j].atom3;
				Pdo.push_back(p);
			}
	
		}
		for(int j=0;j<AspPos.size();j++)
		{
			if(table[LysPos[i].atom1-1].chain_id == table[AspPos[j].atom1-1].chain_id
			 && abs(table[LysPos[i].atom1-1].res_id-table[AspPos[j].atom1-1].res_id)==pdo_level)
			{
				PdoStruct p;
				p.dist_atom1=LysPos[i].atom1;
				p.dist_atom2=AspPos[j].atom1;
				p.dist_atom3=LysPos[i].atom3;
				p.dist_atom4=AspPos[j].atom3;
				p.dihe_atom1=LysPos[i].atom1;
				p.dihe_atom2=LysPos[i].atom2;
				p.dihe_atom3=AspPos[j].atom2;
				p.dihe_atom4=AspPos[j].atom1;
				p.dihe_atom5=LysPos[i].atom3;
				p.dihe_atom6=LysPos[i].atom2;
				p.dihe_atom7=AspPos[j].atom2;
				p.dihe_atom8=AspPos[j].atom3;
				Pdo.push_back(p);
			}
		}
	}
	IntVector F;
	dcd->getPdo(F,Pdo);
	IntVector2 Acor;
	Acor.resize(3);
	Acor[0].resize(F.size());
	Acor[1].resize(F.size());
	Acor[2].resize(F.size());
	for(int i=0;i<Pdo.size();i++)
	{
		int icount=0;
		for(int j=0;j<F.size();j++)
		{
			if(Pdo[i].D2[j]<=critical_distance)
			{
				Acor[1][j]=1;
			}
			else	Acor[1][j]=0;
			if(fabs(Pdo[i].A2[j])<=critical_angle)
			{
				Acor[2][j]=1;
			}
			else	Acor[2][j]=0;
			if(Acor[1][j]&&Acor[2][j]) {icount++;Acor[0][j]=1;} else Acor[0][j]=0;
		}
		if(icount*1.0/F.size()<critical_percent) continue;
		atom a,b;
		a=table[Pdo[i].dist_atom1-1];
		b=table[Pdo[i].dist_atom2-1];
		Log("Chain  = "+a.chain_id+" "+
			"First pdo = "+a.res_name+" "+printNumber(a.res_id)+" "+
			"Second pdo = "+b.res_name+" "+printNumber(b.res_id));
		string filename=
			"pdo_"+a.chain_id+"_"+
			a.res_name+"_"+printNumber(a.res_id)+"_"+
			b.res_name+"_"+printNumber(b.res_id)+".dat";

		string statname=
			"pdo_"+a.chain_id+"_"+
			a.res_name+"_"+printNumber(a.res_id)+"_"+
			b.res_name+"_"+printNumber(b.res_id)+".stat";

		string smoothname=
			"pdo_"+a.chain_id+"_"+
			a.res_name+"_"+printNumber(a.res_id)+"_"+
			b.res_name+"_"+printNumber(b.res_id)+".sda";


		string histnamea=
			"pdo_"+a.chain_id+"_"+
			a.res_name+"_"+printNumber(a.res_id)+"_"+
			b.res_name+"_"+printNumber(b.res_id)+".hista";


		string histnamed=
			"pdo_"+a.chain_id+"_"+
			a.res_name+"_"+printNumber(a.res_id)+"_"+
			b.res_name+"_"+printNumber(b.res_id)+".histd";

		string acorname=
			"pdo_"+a.chain_id+"_"+
			a.res_name+"_"+printNumber(a.res_id)+"_"+
			b.res_name+"_"+printNumber(b.res_id)+".cor";

		PrintAcor(acorname,F,Acor[0],kstep,astep);
		FILE *fp;
		fp=fopen(filename.c_str(),"w");
		if(!fp) Error(WriteError(filename));
		PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
			printDistanceHeader("Dist1")+" "+printDistanceHeader("Dist2")+" "+
			printAngleHeader("Dihe1")+" "+printAngleHeader("Dihe2"));
		for(int j=0;j<F.size();j++)
		{
			PrintLine(fp,printFrame(F[j])+" "+printDistance(Pdo[i].D1[j])+" "+
				printDistance(Pdo[i].D2[j])+" "+printAngle(Pdo[i].A1[j])+" "+
				printAngle(Pdo[i].A2[j]));
		}
		fclose(fp);
		atom d1=table[Pdo[i].dist_atom1-1];
		atom d2=table[Pdo[i].dist_atom2-1];
		atom d3=table[Pdo[i].dist_atom3-1];
		atom d4=table[Pdo[i].dist_atom4-1];
		string header1="#Distance statistics "+printAtom6(d1)+"-"+printAtom6(d2);
		PrintStatFile(statname,"w",F.size()/statcount,header1,F,Pdo[i].D1);
		atom a1=table[Pdo[i].dihe_atom1-1];
		atom a2=table[Pdo[i].dihe_atom2-1];
		atom a3=table[Pdo[i].dihe_atom3-1];
		atom a4=table[Pdo[i].dihe_atom4-1];
		string header3="#Angle statistics "+printAtom6(a1)+"-"+printAtom6(a2)+"-"+
			printAtom6(a3)+"-"+printAtom6(a4);
		PrintAngleStatFile(statname,"a",F.size()/statcount,header3,F,Pdo[i].A1);

		string header2="#Distance statistics "+printAtom6(d3)+"-"+printAtom6(d4);
		PrintStatFile(statname,"a",F.size()/statcount,header2,F,Pdo[i].D2);

		a1=table[Pdo[i].dihe_atom5-1];
		a2=table[Pdo[i].dihe_atom6-1];
		a3=table[Pdo[i].dihe_atom7-1];
		a4=table[Pdo[i].dihe_atom8-1];
		string header4="#Angle statistics "+printAtom6(a1)+"-"+printAtom6(a2)+"-"+
			printAtom6(a3)+"-"+printAtom6(a4);
		PrintAngleStatFile(statname,"a",F.size()/statcount,header4,F,Pdo[i].A2);

		string header5="#Existence statistics";
		string header6="#"+printString("Block",FRAME_WIDTH-1)+" "+
			printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
			printFrameHeader("Total")+" "+printFrameHeader("Distance")+" "+
			printFrameHeader("Angle");
		PrintIntStatFile(statname,"a",F.size()/statcount,header5,header6,F,Acor);

		Data dmin1,dmax1,davg1,dstd1;
		getVectorStatistics(Pdo[i].D1,dmin1,dmax1,davg1,dstd1);
		Data dmin2,dmax2,davg2,dstd2;
		getVectorStatistics(Pdo[i].D2,dmin2,dmax2,davg2,dstd2);
		Data amin1,amax1,aavg1,astd1;
		getAngleVectorStatistics(Pdo[i].A1,amin1,amax1,aavg1,astd1);
		Data amin2,amax2,aavg2,astd2;
		getAngleVectorStatistics(Pdo[i].A2,amin2,amax2,aavg2,astd2);

		if(smooth_flag)
		{
			IntVector sf;
			DoubleVector sd1,sd2,sa1,sa2;
			makeSmooth(F,Pdo[i].D1,sf,sd1,smooth_start,smooth_step);
			makeSmooth(F,Pdo[i].D2,sf,sd2,smooth_start,smooth_step);
			makeAngleSmooth(F,Pdo[i].A1,sf,sa1,smooth_start,smooth_step);
			makeAngleSmooth(F,Pdo[i].A2,sf,sa2,smooth_start,smooth_step);
			fp=fopen(smoothname.c_str(),"w");
			if(!fp) Error(WriteError(smoothname));
			PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
				printDistanceHeader("Dist1")+" "+printDistanceHeader("Dist2")+" "+
				printAngleHeader("Dihe1")+" "+printAngleHeader("Dihe2"));
			for(int j=0;j<sf.size();j++)
			{
				PrintLine(fp,printFrame(sf[j])+" "+printDistance(sd1[j])+" "+
				printDistance(sd2[j])+" "+printAngle(sa1[j])+" "+
				printAngle(sa2[j]));
			}
			fclose(fp);
		}
		if(histflag)
		{
			vector<HistStruct> st_dist1;
			vector<HistStruct> st_dist2;
			vector<HistStruct> st_angle1;
			vector<HistStruct> st_angle2;
			makeHist(Pdo[i].D1,st_dist1,(dmin1<dmin2)?dmin1:dmin2,
				(dmax1>dmax2)?dmax1:dmax2,bindist);
			makeHist(Pdo[i].D2,st_dist2,(dmin1<dmin2)?dmin1:dmin2,
				(dmax1>dmax2)?dmax1:dmax2,bindist);
			makeHist(Pdo[i].A1,st_angle1,-180.0,180.0,bindihe);
			makeHist(Pdo[i].A2,st_angle2,-180.0,180.0,bindihe);

			fp=fopen(histnamed.c_str(),"w");	
			if(!fp) Error(WriteError(histnamed));
			PrintLine(fp,"#"+printString("Center",DISTANCE_WIDTH-1)+" "+
			printPercentHeader("F(%)")+" "+printFrameHeader("F(#)")+" "+
			printPercentHeader("F(%)")+" "+printFrameHeader("F(#)"));
			for(int j=0;j<st_dist1.size();j++)
			{
				PrintLine(fp,printDistance(st_dist1[j].value)+" "+
					printPercent(st_dist1[j].percent)+" "+
					printFrame(st_dist1[j].count)+" "+
					printPercent(st_dist2[j].percent)+" "+
					printFrame(st_dist2[j].count));
			}
			fclose(fp);
			fp=fopen(histnamea.c_str(),"w");	
			if(!fp) Error(WriteError(histnamea));
			PrintLine(fp,"#"+printString("Center",ANGLE_WIDTH-1)+" "+
			printPercentHeader("F(%)")+" "+printFrameHeader("F(#)")+" "+
			printPercentHeader("F(%)")+" "+printFrameHeader("F(#)"));
			for(int j=0;j<st_angle1.size();j++)
			{
				PrintLine(fp,printAngle(st_angle1[j].value)+" "+
					printPercent(st_angle1[j].percent)+" "+
					printFrame(st_angle1[j].count)+" "+
					printPercent(st_angle2[j].percent)+" "+
					printFrame(st_angle2[j].count));
			}
			fclose(fp);
		}
	}
}
