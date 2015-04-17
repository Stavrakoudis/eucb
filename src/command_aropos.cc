# include <global.h>
# include <command_aropos.h>

CommandAropos::CommandAropos(Data sp,Data sd,Data sa )
	:Command("AROPOS")
{
	stack_percent = sp;
	stack_distance = sd;
	stack_angle = sa;
	SideChain.resize(24);

	SideChain[0].res_name="TYR";SideChain[0].atom_name="CD1";
	SideChain[1].res_name="TYR";SideChain[1].atom_name="CD2";
	SideChain[2].res_name="TYR";SideChain[2].atom_name="CZ";

	SideChain[3].res_name="PHE";SideChain[3].atom_name="CD1";
	SideChain[4].res_name="PHE";SideChain[4].atom_name="CD2";
	SideChain[5].res_name="PHE";SideChain[5].atom_name="CZ";

	SideChain[6].res_name="TRP";SideChain[6].atom_name="CD1";
	SideChain[7].res_name="TRP";SideChain[7].atom_name="CE3";
	SideChain[8].res_name="TRP";SideChain[8].atom_name="CZ2";

	posnames.resize(6);
	posnames[0].chain = "ARG";
	posnames[0].nname = "NE";
	posnames[0].hname.resize(1);
	posnames[0].hname[0]="HE";

	posnames[1].chain = "ARG";
	posnames[1].nname = "NH1";
	posnames[1].hname.resize(2);
	posnames[1].hname[0]="HH11";
	posnames[1].hname[1]="HH12";

	posnames[2].chain = "ARG";
	posnames[2].nname = "NH2";
	posnames[2].hname.resize(2);
	posnames[2].hname[0]="HH21";
	posnames[2].hname[1]="HH22";

	posnames[3].chain = "LYS";
	posnames[3].nname = "NZ";
	posnames[3].hname.resize(3);
	posnames[3].hname[0]="HZ1";
	posnames[3].hname[1]="HZ2";
	posnames[3].hname[2]="HZ3";

	posnames[4].chain = "HSP";
	posnames[4].nname = "ND1";	
	posnames[4].hname.resize(1);
	posnames[4].hname[0]="HD1";

	posnames[5].chain = "HSP";
	posnames[5].nname = "NE2";	
	posnames[5].hname.resize(1);
	posnames[5].hname[0]="HE2";
}

void CommandAropos::Run()
{
	Team *stack_team = team;
	string chain;
	int icount=0;
	vector<stack_table> stable;
	vector<AroposStruct> nhpos;

	for(int i=0;i<table.size();i++)
	{
		for(int j=0;j<posnames.size();j++)
		{
			if(table[i].res_name == posnames[j].chain && (stack_team==NULL ||
			stack_team->find(table[i],chain)))
			{
				int s = nhpos.size();
				if(table[i].atom_name == posnames[j].nname)
				{
					nhpos.resize(s+1);
					nhpos[s].npos = i;
					nhpos[s].hpos.resize(0);
				}
				else
				{
					for(int k=0;k<posnames[j].hname.size();k++)
					{
						if(table[i].atom_name == posnames[j].hname[k])
						{
							nhpos[s-1].hpos.push_back(i);
						}
					}
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
	vector<DoubleVector> Dist;
	vector<DoubleVector> Angle;
	vector<stack_table> Stable;
	vector<AroposStruct> Nhpos;
	dcd->setFrameStep(smart_skip);
	dcd->getAropos(f,Dist,Angle,stable,nhpos,smart_distance,stack_angle);
	for(int i=0;i<nhpos.size();i++)
	{
		for(int j=0;j<stable.size();j++)
		{
			int found_flag=0;
			for(int k=0;k<f.size();k++)
			{
				if(Dist[i*stable.size()+j][k]<=smart_distance 
				 && Angle[i*stable.size()+j][k]<=stack_angle) 
				{
					found_flag=1;
					break;
				}
			}
			if(found_flag)
			{
				Stable.push_back(stable[j]);
				Nhpos.push_back(nhpos[i]);
			}
		}
	}

	
	dcd->setFrameStep(step);
	dcd->getAropos(f,Dist,Angle,Stable,Nhpos,smart_distance,stack_angle);

	IntVector2 Acor;
	Acor.resize(3);
	Acor[0].resize(f.size());
	Acor[1].resize(f.size());
	Acor[2].resize(f.size());
	for(int i=0;i<Nhpos.size();i++)
	{
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
				if(Dist[i*Stable.size()+j][k]<=stack_distance 
				 && Angle[i*Stable.size()+j][k]<=stack_angle) 
				{
					counter++;
					Acor[0][k]=1;
				}
				else	Acor[0][k]=0;
				if(Dist[i*Stable.size()+j][k]<=stack_distance) 
				{
					icount1++;
					Acor[1][k]=1;
				}
				else	Acor[1][k]=0;
				if(Angle[i*Stable.size()+j][k]<=stack_angle) 
				{
					icount2++;
					Acor[2][k]=1;
				}
				else	
					Acor[2][k]=0;
			}
			if(counter*1.0/f.size()>=stack_percent)
			{
				string filename="aropos_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".dat";

				string acorname="aropos_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".cor";
				PrintAcor(acorname,f,Acor[0],kstep,astep);


				FILE *fp=fopen(filename.c_str(),"w");	
				for(int k=0;k<f.size();k++)
				{
					PrintLine(fp,printFrame(f[k])+" "+
						printDistance(Dist[i*Stable.size()+j][k])+" "+
						printAngle(Angle[i*Stable.size()+j][k]));
				}
				fclose(fp);
				Data dmin,dmax,davg,dstd;
				if(smooth_flag)
				{
					filename="aropos_"+printAtom2(table[Stable[j].id1-1])+"_"+
							  printAtom2(table[Nhpos[i].npos])+".sda";
					IntVector sf;
					DoubleVector sd,sa;
					makeSmooth(f,Dist[i*Stable.size()+j],sf,sd,smooth_start,smooth_step);
					makeAngleSmooth(f,Angle[i*Stable.size()+j],sf,sa,smooth_start,smooth_step);
					fp=fopen(filename.c_str(),"w");
					for(int k=0;k<sf.size();k++)
						PrintLine(fp,printFrame(sf[k])+" "+printDistance(sd[k])+" "+
							printAngle(sa[k]));
					fclose(fp);
				}
				Data amin,amax,aavg,astd;
				getVectorStatistics(Dist[i*Stable.size()+j],dmin,dmax,davg,dstd);
				getAngleVectorStatistics(Angle[i*Stable.size()+j],amin,amax,aavg,astd);
				string statname="aropos_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".stat";

				string header1="#Distance statistics";
				string header2="#Angle statistics";
				PrintStatFile(statname,"w",f.size()/statcount,header1,f,Dist[i*Stable.size()+j]);
				PrintAngleStatFile(statname,"a",f.size()/statcount,header2,f,Angle[i*Stable.size()+j]);
				string header3="#Existance";
				string header4="#"+printString("Block",FRAME_WIDTH-1)+" "+
					printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
					printFrameHeader("Icount")+" "+printFrameHeader("Dcount")+" "+
					printFrameHeader("Acount");
				PrintIntStatFile(statname,"a",f.size()/statcount,header3,header4,f,Acor);

				vector<HistStruct> st_dist;
				vector<HistStruct> st_angle;
				makeHist(Dist[i*Stable.size()+j],st_dist,dmin,dmax,bindist);
				makeHist(Angle[i*Stable.size()+j],st_angle,0.0,90.0,binangle);
				string histnamea="aropos_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".hista";
				fp=fopen(histnamea.c_str(),"w");
				for(int k=0;k<st_angle.size();k++)
				{
					PrintLine(fp,printAngle(st_angle[k].value)+" "+
					 printPercent(st_angle[k].percent)+" "+
					 printFrame(st_angle[k].count));
				}
				fclose(fp);
				string histnamed="aropos_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".histd";
				fp=fopen(histnamed.c_str(),"w");
				for(int k=0;k<st_dist.size();k++)
				{
					PrintLine(fp,printDistance(st_dist[k].value)+" "+
					 printPercent(st_dist[k].percent)+" "+
					 printFrame(st_dist[k].count));
				}
				fclose(fp);
				string histname2="aropos_"+
					printAtom2(table[Stable[j].id1-1])+"_"+
					printAtom2(table[Nhpos[i].npos])+".hist2";
				fp=fopen(histname2.c_str(),"w");
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
