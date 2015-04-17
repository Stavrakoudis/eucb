# include <global.h>
# include <math.h>
# include <command_aronh.h>

CommandAronh::CommandAronh(Data sp,Data sd,Data sa )
	:Command("ARONH")
{
	stack_percent = sp;
	stack_distance = sd;
	stack_angle = sa;
	SideChain.resize(9);
	/*
	SideChain[0].res_name="ARG";SideChain[0].atom_name="NE";
	SideChain[1].res_name="ARG";SideChain[1].atom_name="NH1";
	SideChain[2].res_name="ARG";SideChain[2].atom_name="NH2";
	*/
	SideChain[0].res_name="TYR";SideChain[0].atom_name="CD1";
	SideChain[1].res_name="TYR";SideChain[1].atom_name="CD2";
	SideChain[2].res_name="TYR";SideChain[2].atom_name="CZ";

	SideChain[3].res_name="PHE";SideChain[3].atom_name="CD1";
	SideChain[4].res_name="PHE";SideChain[4].atom_name="CD2";
	SideChain[5].res_name="PHE";SideChain[5].atom_name="CZ";

	SideChain[6].res_name="TRP";SideChain[6].atom_name="CD1";
	SideChain[7].res_name="TRP";SideChain[7].atom_name="CE3";
	SideChain[8].res_name="TRP";SideChain[8].atom_name="CZ2";
	/*
	SideChain[12].res_name="HSD";SideChain[12].atom_name="CG";
	SideChain[13].res_name="HSD";SideChain[13].atom_name="CE1";
	SideChain[14].res_name="HSD";SideChain[14].atom_name="NE2";

	SideChain[15].res_name="HSE";SideChain[15].atom_name="CG";
	SideChain[16].res_name="HSE";SideChain[16].atom_name="CE1";
	SideChain[17].res_name="HSE";SideChain[17].atom_name="NE2";

	SideChain[18].res_name="HSP";SideChain[18].atom_name="CG";
	SideChain[19].res_name="HSP";SideChain[19].atom_name="CE1";
	SideChain[20].res_name="HSP";SideChain[20].atom_name="NE2";

	SideChain[21].res_name="PRO";SideChain[21].atom_name="CB";
	SideChain[22].res_name="PRO";SideChain[22].atom_name="CG";
	SideChain[23].res_name="PRO";SideChain[23].atom_name="CD";
	*/
}

void CommandAronh::Run()
{
	Team *stack_team = team;
	string chain;
	int icount=0;
	vector<stack_table> stable;
	vector<NHStruct> nhpos;

	for(int i=0;i<table.size();i++)
	{
		if(table[i].atom_name == "N" && (stack_team == NULL || stack_team->find(table[i],chain)))
		{
			int s=nhpos.size();
			if(!s)
			{
				nhpos.resize(s+1);
				nhpos[s].npos = i;
				nhpos[s].hpos = -1;
			}
			else
			{
				if(nhpos[s-1].npos==-1) nhpos[s-1].npos=i;
				else
				{
					nhpos.resize(s+1);
					nhpos[s].npos = i;
					nhpos[s].hpos = -1;
				}
			}
		}
		if(table[i].atom_name == "HN" && (stack_team == NULL || stack_team->find(table[i],chain)))
		{
			int s=nhpos.size();
			if(!s)
			{
				nhpos.resize(s+1);
				nhpos[s].hpos = i;
				nhpos[s].npos = -1;
			}
			else
			{
				if(nhpos[s-1].hpos==-1) nhpos[s-1].hpos=i;
				else
				{
					nhpos.resize(s+1);
					nhpos[s].hpos = i;
					nhpos[s].npos = -1;
				}
			}
		}
		for(int j=0;j<SideChain.size();j++)
		{
			if(table[i].res_name==SideChain[j].res_name && 
			  table[i].atom_name==SideChain[j].atom_name &&
			  (stack_team==NULL || stack_team->find(table[i],chain)))
			{
				icount++;
				if(icount==1) {int s=stable.size();stable.resize(s+1);}
				if(icount==1) stable[stable.size()-1].id1=i+1;
				if(icount==2) stable[stable.size()-1].id2=i+1;
				if(icount==3) stable[stable.size()-1].id3=i+1;
				if(icount==3) icount=0;
			}
		}
	}
	IntVector f;
	DoubleVector2 Dist;
	DoubleVector2 Dist2;
	DoubleVector2 Angle;

	dcd->setFrameStep(smart_skip);
	dcd->getAronh(f,Dist,Dist2,Angle,stable,nhpos,stack_distance,stack_angle);
	vector<stack_table> Stable;
	vector<NHStruct> Nhpos;
	for(int i=0;i<nhpos.size();i++)
	{
		if(nhpos[i].npos==-1 || nhpos[i].hpos==-1) continue;
		for(int j=0;j<stable.size();j++)
		{
			int found_flag=0;
			for(int k=0;k<f.size();k++)
			{
			if(Dist2[i*stable.size()+j][k]>Dist[i*stable.size()+j][k]) continue;
			if(Dist[i*stable.size()+j][k]<=smart_distance 
			 && (Angle[i*stable.size()+j][k])>=stack_angle
			 && Dist2[i*stable.size()+j][k]<=Dist[i*stable.size()+j][k]) 
			{found_flag=1;break;}
			}
			if(found_flag) 
			{
				Stable.push_back(stable[j]);
				Nhpos.push_back(nhpos[i]);
			}
		}
	}

	
	dcd->setFrameStep(step);
	dcd->getAronh(f,Dist,Dist2,Angle,Stable,Nhpos,stack_distance,stack_angle);
	IntVector2 Acor;
	Acor.resize(3);
	Acor[0].resize(f.size());
	Acor[1].resize(f.size());
	Acor[2].resize(f.size());
	for(int i=0;i<Nhpos.size();i++)
	{
		if(Nhpos[i].npos==-1 || Nhpos[i].hpos==-1) continue;
		for(int k1=0;k1<3;k1++)
			for(int k2=0;k2<f.size();k2++)
				Acor[k1][k2]=0;
		for(int j=0;j<Stable.size();j++)
		{
			int counter=0;
			int icount1=0;
			int icount2=0;
			for(int k=0;k<f.size();k++)
			{
				if(Dist2[i*Stable.size()+j][k]>Dist[i*Stable.size()+j][k]) continue;
				if(Dist[i*Stable.size()+j][k]<=stack_distance 
				 && Angle[i*Stable.size()+j][k]>=stack_angle
				 && Dist2[i*Stable.size()+j][k]<=Dist[i*Stable.size()+j][k]) 
				{
					counter++;
					Acor[0][k]=1;
				}
				else Acor[0][k]=0;
				if(Dist[i*Stable.size()+j][k]<=stack_distance) 
				{
					Acor[1][k]=1;
					icount1++;
				}
				else
					Acor[1][k]=0;
				if((Angle[i*Stable.size()+j][k])>=stack_angle) 
				{
					icount2++;
					Acor[2][k]=1;
				}
				else
					Acor[2][k]=0;
			}
			if(counter*1.0/f.size()>=stack_percent)
			{
				string filename="aroHN_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".dat";
				string acorname="aroHN_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".cor";

				PrintLine(fplog,"Create data file "+filename);


				FILE *fp=fopen(filename.c_str(),"w");	
				if(!fp) Error(WriteError(filename));
				PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
					printString("Exist",5)+" "+
					printDistanceHeader("N-Ring")+" "+
					printDistanceHeader("H-Ring")+" "+
					printAngleHeader("Angle"));
				for(int k=0;k<f.size();k++)
				{
					int exist=Acor[0][k];
				if(Dist[i*Stable.size()+j][k]<=stack_distance 
				 && Angle[i*Stable.size()+j][k]>=stack_angle
				 && Dist2[i*Stable.size()+j][k]<=Dist[i*Stable.size()+j][k])  exist=1; else exist=0;
					PrintLine(fp,printFrame(f[k])+" "+
						printNumber(exist,5)+" "+
						printDistance(Dist[i*Stable.size()+j][k])+" "+
						printDistance(Dist2[i*Stable.size()+j][k])+" "+
						printAngle(Angle[i*Stable.size()+j][k]));
				}
				fclose(fp);
				if(smooth_flag)
				{
					filename="aroHN_"+printAtom2(table[Stable[j].id1-1])+"_"+
							  printAtom2(table[Nhpos[i].npos])+".sda";
					PrintLine(fplog,"Create smooth data file "+filename);
					IntVector sf;
					DoubleVector sd,sa,sd2;
					makeSmooth(f,Dist[i*Stable.size()+j],sf,sd,smooth_start,smooth_step);
					makeSmooth(f,Dist2[i*Stable.size()+j],sf,sd2,smooth_start,smooth_step);
					makeAngleSmooth(f,Angle[i*Stable.size()+j],sf,sa,smooth_start,smooth_step);
					fp=fopen(filename.c_str(),"w");
					if(!fp) Error(WriteError(filename));
					for(int k=0;k<sf.size();k++)
						PrintLine(fp,printFrame(sf[k])+" "+printDistance(sd[k])+" "+
							printDistance(sd2[k])+" "+
							printAngle(sa[k]));
					fclose(fp);
				}
				Data dmin,dmax,davg,dstd;
				Data amin,amax,aavg,astd;
				getVectorStatistics(Dist[i*Stable.size()+j],dmin,dmax,davg,dstd);
				getAngleVectorStatistics(Angle[i*Stable.size()+j],amin,amax,aavg,astd);
				string statname="aroHN_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".stat";
				PrintLine(fplog,"Create statistics file "+statname);
				string header1="#Distance statistics";
				string header2="#Angle statistics";
				PrintStatFile(statname,"w",f.size()/statcount,header1,f,Dist[i*Stable.size()+j]);
				string header1_1="#Distance2 statistics";
				PrintStatFile(statname,"a",f.size()/statcount,header1_1,f,Dist2[i*Stable.size()+j]);

				PrintAngleStatFile(statname,"a",f.size()/statcount,header2,f,Angle[i*Stable.size()+j]);
				string header3="#Existance";
				string header4="#"+printString("Block",FRAME_WIDTH-1)+" "+
					printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
					printFrameHeader("Icount")+" "+printFrameHeader("Dcount")+" "+
					printFrameHeader("Acount");
				PrintIntStatFile(statname,"a",f.size()/statcount,header3,header4,f,Acor);
				PrintAcor(acorname,f,Acor[0],kstep,astep);
				if(histflag)
				{
				vector<HistStruct> st_dist;
				vector<HistStruct> st_angle;
				makeHist(Dist[i*Stable.size()+j],st_dist,dmin,dmax,bindist);
				makeHist(Angle[i*Stable.size()+j],st_angle,0.0,180.0,binangle);
				string histnamea="aroHN_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".hista";
				PrintAngleHist(histnamea,st_angle);
				string histnamed="aroHN_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".histd";
				PrintHist(histnamed,st_dist);
				string histname2="aroHN_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".hist2";
				fp=fopen(histname2.c_str(),"w");
				if(!fp) Error(WriteError(histname2));
				for(int k=0;k<st_dist.size();k++)
				{
					for(int m=0;m<st_angle.size();m++)
					{
					int icount=0;
					for(int l=0;l<f.size();l++)
					{
					if(Dist[i*Stable.size()+j][l]>=st_dist[k].value-bindist/2.0 
					 && Dist[i*Stable.size()+j][l]<=st_dist[k].value+bindist/2.0 
					 && Angle[i*Stable.size()+j][l]>=st_angle[m].value-binangle/2.0
					 && Angle[i*Stable.size()+j][l]<=st_angle[m].value+binangle/2.0)
							icount++;
						}
					PrintLine(fp,printDistance(st_dist[k].value)+" "+
					 printAngle(st_angle[m].value)+" "+printFrame(icount));
					}
				   PrintLine(fp,"");
				}
				fclose(fp);
				}
			}
		}
	}
}
