# include <global.h>
# include <donors.h>
# include <math.h>
# include <command_watbridge.h>


CommandWatbridge::CommandWatbridge()
	:Command("WATBRIDGE")
{
}

static	int isDonor(int id)
{
	for(int i=0;i<Donor.size();i++) if(Donor[i]==id) return i+1;
	return -1;
}

static	int isAcceptor(int id)
{
	for(int i=0;i<Acceptor.size();i++) if(Acceptor[i]==id) return i+1;
	return -1;
}


typedef struct
{
	int		id;
	int	 	water;
}wstruct;

/*	Implement the -watbridge option.
 * */
void	CommandWatbridge::Run()
{
	IntVector frame;
	Log("Critical percent  =  "+printPercent(backbone_critical_percent));
	Log("Critical distance =  "+printDistance(backbone_critical_distance));
	Log("Critical angle    =  "+printAngle(backbone_critical_angle));
	
	IntVector List;
	IntVector2 List2;
	IntVector2 Water2;
	IntVector Water;
	string chain;
	IntVector List3;
	Team *watbridge_team1=team1;
	Team *watbridge_team2=team2;
	for(int i=0;i<Donor.size();i++)
	{
		if(!watbridge_team1 || watbridge_team1->find(table[Donor[i]-1],chain)) 
			if(isin(List,Donor[i])==-1)
				List.push_back(Donor[i]);
		if(!watbridge_team2 || watbridge_team2->find(table[Donor[i]-1],chain)) 
			if(isin(List3,Donor[i])==-1)
				List3.push_back(Donor[i]);
	}
	for(int i=0;i<Acceptor.size();i++)
	{
		if(!watbridge_team1 || watbridge_team1->find(table[Acceptor[i]-1],chain))
			if(isin(List,Acceptor[i])==-1) 
				List.push_back(Acceptor[i]);
		if(!watbridge_team2 || watbridge_team2->find(table[Acceptor[i]-1],chain))
			if(isin(List3,Acceptor[i])==-1)
				List3.push_back(Acceptor[i]);
	}
	for(int i=0;i<table.size();i++)
	{
		if(table[i].atom_name=="OH2" && table[i].res_name=="TIP3")
			Water.push_back(table[i].atom_id);
	}
	List2.resize(List.size());
	Water2.resize(List.size());
	dcd->setFrameStep(smart_skip);
	
	vector<NoeStruct> noe;
	for(int i=0;i<List.size();i++)
	{
		for(int j=0;j<List3.size();j++)
		{
			NoeStruct n;
			n.atom1=List[i];
			n.atom2=List3[j];
			noe.push_back(n);
		}
	}
	dcd->getNoe(frame,noe);
	int icount=0;
	for(int i=0;i<List.size();i++)
	{
		for(int j=0;j<List3.size();j++)
		{
			DoubleVector D=noe[icount].D;
			Data p=getBelowPercent(D,smart_distance);
		
			if(p>=backbone_critical_percent)
				List2[i].push_back(List3[j]);
			icount++;
		}
	}

	for(int i=0;i<List.size();i++)
	{
		if(List2[i].size()==0)
		{
			for(int j=i;j<List.size()-1;j++)
			{
				List[j]=List[j+1];
				List2[j]=List2[j+1];
			}
			int s=List.size();
			List.resize(s-1);
			List2.resize(s-1);
			Water2.resize(s-1);
		}
	}
	
	printf ("List  :   %d\n", List.size() );  // here
	printf ("Water :   %d\n", Water.size() ); // here

	noe.resize(0);
	for(int i=0;i<List.size();i++)
	{
		for(int j=0;j<Water.size();j++)
		{
			NoeStruct n;
			n.atom1=List[i];
			n.atom2=Water[j];
			noe.push_back(n);
		}
	}
	printf ("noe    :   %d\n", noe.size() ); // here
	printf ("frame  :   %d\n", frame.size() ); // here
	dcd->getNoe(frame,noe);
	printf ("NOE\n"); //here
	icount=0;
	for(int i=0;i<List.size();i++)
	{
		for(int j=0;j<Water.size();j++)
		{
			DoubleVector D=noe[icount].D;
			Data p=getBelowPercent(D,1.5 * smart_distance);
			if(fabs(p)>1e-7) // what?
			//if(p>=backbone_critical_percent)
			{
				Water2[i].push_back(Water[j]);
			}
			icount++;
		}
	}
	printf ("OK 1\n"); //here
	noe.resize(0);

	IntVector List4;
	IntVector2 List4_2;
	IntVector2 Water4_2;
	for(int i=0;i<List.size();i++)
	{
		for(int j=0;j<List2[i].size();j++)
		{
			int id=List2[i][j];
			List4.push_back(id);
			IntVector empty;
			empty.resize(0);
			List4_2.push_back(empty);
			Water4_2.push_back(Water2[i]);
		}
	}
	for(int i=0;i<List4.size();i++)
	{
		List.push_back(List4[i]);
		List2.push_back(List4_2[i]);
		Water2.push_back(Water4_2[i]);
	}
	printf ("OK 2\n"); //here
	
	dcd->setFrameStep(step);
	dcd->setFrameEnd(last);

	vector<HbondStruct> hb;
	vector<AcceptorStruct> ha;

	IntVector DonorAcceptor;
	DonorAcceptor.resize(List.size());
	for(int i=0;i<List.size();i++)
	{
		DonorAcceptor[i]=0;
//		if(List2[i].size()==0) continue;
		int donor_pos = isDonor(List[i]);
		if(donor_pos!=-1)
		{
			DonorAcceptor[i]=1;
			HbondStruct st;
			st.donor_id=List[i];
			st.Acceptor=Water2[i];
			st.Hydro=DonorCon[donor_pos-1];
			hb.push_back(st);
		}
		int acceptor_pos=isAcceptor(List[i]);	
		if(acceptor_pos!=-1)
		{
			DonorAcceptor[i]=2;
			AcceptorStruct st;
			st.acceptor_id=List[i];
			st.Donor.resize(Water2[i].size());
			for(int j=0;j<st.Donor.size();j++)
			{
				st.Donor[j]=isDonor(Water2[i][j])-1;
			}
			st.Hydro.resize(Water2[i].size());
			for(int j=0;j<st.Hydro.size();j++)
				st.Hydro[j]=DonorCon[st.Donor[j]];
			ha.push_back(st);
		}
	}
	printf ("OK 3\n"); //here
	dcd->hbonds(frame,hb);
	printf ("OK 3a\n"); //here
	dcd->acceptors(frame,ha,Donor);
	vector<wstruct> ws;
	for(int i=0;i<hb.size();i++)
	{
		DoubleVector2 TotalDist;
		DoubleVector2 TotalAngle;
		TotalDist=hb[i].dist;
		TotalAngle=hb[i].angle;
		for(int j=0;j<hb[i].Acceptor.size();j++)
		{
			int icount=0;
			for(int k=0;k<frame.size();k++)
			{
				double d=TotalDist[j][k];
				double a=TotalAngle[j][k];
				if(d<=backbone_critical_distance && 
				   a>=backbone_critical_angle)
					icount++;
			}
			if(icount*1.0/frame.size()>=backbone_critical_percent) 
			{
				wstruct st;
				st.id=hb[i].donor_id;
				st.water=hb[i].Acceptor[j];		
				ws.push_back(st);
			}
		}
	}
	printf ("OK 4\n"); //here

	for(int i=0;i<ha.size();i++)
	{
		for(int j=0;j<ha[i].Donor.size();j++)
		{
			int icount=0;
			for(int k=0;k<frame.size();k++)
			{
				double d=ha[i].dist[j][k];
				double a=ha[i].angle[j][k];
                                if(d<=backbone_critical_distance &&a>=backbone_critical_angle)
                                        icount++;
			}
			if(icount*1.0/frame.size()>=backbone_critical_percent)
			{
				wstruct st;
				st.id=ha[i].acceptor_id;
				st.water=Donor[ha[i].Donor[j]];
				ws.push_back(st);
			}
		}
	}

	IntVector Found;
	IntVector Found2;
	IntVector2 WaterFound;
	WaterFound.resize(0);

	for(int i=0;i<ws.size();i++)
	{
		int id1=ws[i].id;
		for(int j=0;j<ws.size();j++)
		{
			if(i==j) continue;
			int id2=ws[j].id;
			int lpos=isin(List,id1);
			if(isin(List2[lpos],id2)!=-1 && ws[i].water==ws[j].water)
			{
				int ifound=-1;
				for(int k=0;k<Found.size();k++)
				{
					if(Found[k]==id1 && Found2[k]==id2) {ifound=k;break;}
				}
				if(ifound==-1)
				{
					Found.push_back(id1);	
					Found2.push_back(id2);
					int s=WaterFound.size();
					WaterFound.resize(s+1);
					WaterFound[s].push_back(ws[i].water);
				}
				else
				{
					int wpos=isin(WaterFound[ifound],ws[i].water);
					if(wpos==-1)
						WaterFound[ifound].push_back(ws[i].water);
				}
			}	
		}
	}

	for(int i=0;i<Found.size();i++)
	{
		atom atom1 = table[Found[i]-1];
		atom atom2 = table[Found2[i]-1];
		string filename ="wb_"+printAtom1(atom1)+"_"+printAtom1(atom2)+".dat";
		string statname ="wb_"+printAtom1(atom1)+"_"+printAtom1(atom2)+".stat";
		string acorname ="wb_"+printAtom1(atom1)+"_"+printAtom1(atom2)+".cor";
		int id1 = Found[i];
		int id2 = Found2[i];
		int ipos=0;
		int jpos=0;

		for(int j=0;j<List.size();j++)
		{
			if(List[j]==id1) {ipos=j;break;}
		}
		for(int j=0;j<List.size();j++)
		{
			if(List[j]==id2) {jpos=j;break;}
		}
		IntVector2 TotalAcor;
		TotalAcor.resize(WaterFound[i].size());
		for(int j=0;j<TotalAcor.size();j++)
		{
			TotalAcor[j].resize(frame.size());
			for(int k=0;k<frame.size();k++) TotalAcor[j][k]=0;
		}
		for(int j=0;j<WaterFound[i].size();j++)
		{
			int water_id=WaterFound[i][j];
			IntVector2 Acor;
			Acor.resize(3);
			Acor[0].resize(frame.size());
			Acor[1].resize(frame.size());
			Acor[2].resize(frame.size());

			if(DonorAcceptor[ipos]==1)
			{
				for(int k=0;k<hb.size();k++)
				{
					if(hb[k].donor_id==Found[i])
					{
						for(int m=0;m<hb[k].Acceptor.size();m++)
						{
							if(hb[k].Acceptor[m]!=water_id) continue;
							for(int n=0;n<frame.size();n++)
							{
							  double d=hb[k].dist[m][n];
							  double a=hb[k].angle[m][n];
							  if(d<=backbone_critical_distance &&
							     a>=backbone_critical_angle) 
							    Acor[1][m]=1; else Acor[1][m]=0;
							}
						}
					}
				}
			}
			else
			{
				for(int k=0;k<ha.size();k++)
				{
					if(ha[k].acceptor_id==Found[i])
					{
						for(int m=0;m<ha[k].Donor.size();m++)
						{
							if(Donor[ha[k].Donor[m]]!=water_id) continue;
							for(int n=0;n<frame.size();n++)
							{
							  double d=ha[k].dist[m][n];
							  double a=ha[k].angle[m][n];
							  if(d<=backbone_critical_distance &&
							     a>=backbone_critical_angle) 
							    Acor[1][n]=1; else Acor[1][n]=0;
							}
						}
					}
				}
			}
			if(DonorAcceptor[jpos]==1)
			{
				for(int k=0;k<hb.size();k++)
				{
					if(hb[k].donor_id==Found2[i])
					{
						for(int m=0;m<hb[k].Acceptor.size();m++)
						{
							if(hb[k].Acceptor[m]!=water_id) continue;
							for(int n=0;n<frame.size();n++)
							{
							  double d=hb[k].dist[m][n];
							  double a=hb[k].angle[m][n];
							  if(d<=backbone_critical_distance &&
							     a>=backbone_critical_angle) 
							    Acor[2][n]=1; else Acor[2][n]=0;
							}
						}
					}
				}
			}
			else
			{
				for(int k=0;k<ha.size();k++)
				{
					if(ha[k].acceptor_id==Found2[i])
					{
						for(int m=0;m<ha[k].Donor.size();m++)
						{
							if(Donor[ha[k].Donor[m]]!=water_id) continue;
							for(int n=0;n<frame.size();n++)
							{
							  double d=ha[k].dist[m][n];
							  double a=ha[k].angle[m][n];
							  if(d<=backbone_critical_distance &&
							     a>=backbone_critical_angle) 
							    Acor[2][n]=1; else Acor[2][n]=0;
							}
						}
					}
				}
			}
			for(int m=0;m<frame.size();m++) 
			{
				Acor[0][m]=Acor[1][m]*Acor[2][m];
				TotalAcor[j][m]+=Acor[0][m];
			}
		}
		int dcount=0;
		for(int k=0;k<frame.size();k++)
		{
			for(int m=0;m<TotalAcor.size();m++)
				if(TotalAcor[m][k]) {dcount++;break;}
		}
		if(dcount*1.0/frame.size()<backbone_critical_percent) continue;
		Log("Found water bridge "+printAtom3(atom1)+" "+printAtom3(atom2)+" "+printFrame(dcount));
		FILE *fp=fopen(filename.c_str(),"w");
		if(!fp) psf_error("Can not open "+filename+" for writing");

		PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+printFrameHeader("Count")+" "+
			printFrameHeader("AtomId"));	
		for(int m=0;m<frame.size();m++)
		{
			int d=0;
			for(int n=0;n<TotalAcor.size();n++) d+=TotalAcor[n][m];
			Print(fp,printFrame(frame[m])+" "+printFrame(d)+" ");
			for(int n=0;n<TotalAcor.size();n++)
			{
				if(TotalAcor[n][m]) 
				Print(fp,printFrame(table[WaterFound[i][n]-1].atom_id)+" ");
			}
			PrintLine(fp,"");
			
		}
		fclose(fp);
		IntVector2 fAcor;
		fAcor.resize(1);
		fAcor[0].resize(frame.size());
		for(int m=0;m<frame.size();m++) fAcor[0][m]=0;
		for(int m=0;m<frame.size();m++) 
			for(int n=0;n<TotalAcor.size();n++)
				fAcor[0][m]+=TotalAcor[n][m];
		string header1="#Existence statistics";
		
		string header2="#"+printString("Block",FRAME_WIDTH-1)+" "+
			printFrameHeader("From")+" "+printFrameHeader("To")+" "+
			printFrameHeader("Count");
		PrintAcor(acorname,frame,fAcor[0],kstep,astep);
		PrintIntStatFile(statname,"w",frame.size()/statcount,header1,header2,frame,fAcor);
	}
}

CommandWatbridge::~CommandWatbridge()
{
}

