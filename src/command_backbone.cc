# include <global.h>
# include <donors.h>
# include <command_backbone.h>
# include <counter.h>
# include <tex.h>

CommandBackbone::CommandBackbone() 
	:Command("HBONDS")
{
}

void	CommandBackbone::Run()
{
	vector<int> CopyDonor;
	vector<int> CopyAcceptor;
	vector<IntVector> CopyDonorCon;
	vector<int> f;

	vector<IntVector> ResAcc;
	string donor_chain="",acceptor_chain="";
	dcd->setFrameStep(smart_skip);
	Team *backBoneTeam = team;
		
	Tex myTex("hbonds.tex");
	vector<string> s;
        s.resize(2);
        s[0]="Acceptor";
        s[1]="Donor";
        myTex.setHeader(s);
        myTex.setCaption("Hydrogen bonds");
	Counter myCounter;

	for(int i=0;i<Donor.size();i++)
	{
		if(backBoneTeam == NULL || 
		backBoneTeam->find(table[Donor[i]-1],donor_chain))
		{
			int s=CopyDonor.size();
			CopyDonor.push_back(Donor[i]);
			CopyDonorCon.resize(s+1);
			CopyDonorCon[s].resize(DonorCon[i].size());
			for(int j=0;j<DonorCon[i].size();j++)
				CopyDonorCon[s][j]=DonorCon[i][j];
		}
	}
	for(int i=0;i<Acceptor.size();i++)
	{
		if(backBoneTeam==NULL ||  backBoneTeam->find(table[Acceptor[i]-1],acceptor_chain))
		{
			CopyAcceptor.push_back(Acceptor[i]);
		}
	}
	int total=0;

	IntVector X1;
	IntVector X2;
	IntVector HydroFlag;
	IntVector	  ResDonor;
	IntVector2 ResAcceptor;
	DoubleVector2 TotalDist,TotalAngle;
	DoubleVector2 TotalDist2;
	vector<IntVector> Flag;
	dcd->hbondsCheck(smart_distance,backbone_critical_angle,backbone_critical_percent,
		CopyDonor,CopyDonorCon,CopyAcceptor,Flag);
	vector<HbondStruct> st;
	for(int i=0;i<Flag.size();i++)
	{
		for(int j=0;j<Flag[i].size();j++)
		{
			if(Flag[i][j])
			{
				int pos=isin(ResDonor,i); //bazoume tin thesi tou donor kai oxi ayton.
				if(pos==-1)
				{
				int s = ResDonor.size();
				ResDonor.resize(s+1);		
				ResDonor[s]=i;
				ResAcceptor.resize(s+1);
				ResAcceptor[s].resize(1);
				
				ResAcceptor[s][0]=CopyAcceptor[j];
				}
				else
				ResAcceptor[pos].push_back(CopyAcceptor[j]);
			}
		}
	}
	dcd->setFrameStep(step);
	st.resize(ResDonor.size());

        for(int i=0;i<st.size();i++)
        {
                st[i].donor_id=CopyDonor[ResDonor[i]];
                st[i].Acceptor.resize(ResAcceptor[i].size());
                for(int j=0;j<ResAcceptor[i].size();j++)
                        st[i].Acceptor[j]=ResAcceptor[i][j];
                st[i].Hydro=CopyDonorCon[ResDonor[i]];
        }
        dcd->hbonds(f,st);


	int bb_bb_counter=0,bb_sc_counter=0,sc_sc_counter=0;
	int bb_bb_frame_counter=0,bb_sc_frame_counter=0,sc_sc_frame_counter=0;
	IntVector2 Acor;
	int donor_id;
	Acor.resize(3);
	Acor[0].resize(f.size());
	Acor[1].resize(f.size());
	Acor[2].resize(f.size());


	for(int i=0;i<st.size();i++)
	{
		donor_id=st[i].donor_id;
		IntVector countAcc;
		countAcc.resize(f.size());
		for(int j=0;j<f.size();j++) countAcc[j]=0;

		X1=CopyDonorCon[ResDonor[i]];
		X2=ResAcceptor[i];
		TotalDist=st[i].dist;
		TotalAngle=st[i].angle;
		TotalDist2=st[i].HydroAcc;
		HydroFlag = st[i].HydroFlag;
		int count1;
		for(int k1=0;k1<3;k1++) 
			for(int k2=0;k2<f.size();k2++)
				Acor[k1][k2]=0;
		for(int j=0;j<X2.size();j++)
		{
			int acceptor_id=X2[j];
		 	count1=0;
			for(int k=0;k<f.size();k++)
			{
				if(TotalDist[j][k]<=backbone_critical_distance) Acor[1][k]=1; else Acor[1][k]=0;
				if(TotalAngle[j][k]>=backbone_critical_angle) Acor[2][k]=1; else Acor[2][k]=0;
				
				if(TotalDist[j][k]<=backbone_critical_distance && 
					TotalAngle[j][k]>=backbone_critical_angle) 
				{
					Acor[0][k]=1;
					count1++;
				}
				else 	Acor[0][k]=0;
				
			}
			Data p1=count1*1.0/f.size();
			if(p1>=backbone_critical_percent)
			{
				string filename,statname,histname,acorname;
				atom a,b;
				a = table[donor_id-1];
				b = table[acceptor_id-1];	
				int s1=printAtom4(a).size();
				int s2=printAtom4(b).size();
				string ss=printAtom4(a);
				int kk;
				for(kk=0;kk<18-s1;kk++)
					ss+=" ";
				ss+=printAtom4(b);
				for( kk=0;kk<18-s2;kk++)
					ss+=" ";
				Log(ss+" "+printPercent(p1));
				if(isBackbone(a) && isBackbone(b)) 
				{
					bb_bb_counter++;	
					bb_bb_frame_counter=p1;
				}
				if(isBackbone(a) && isSidechain(b)) 
				{
					bb_sc_counter++;
					bb_sc_frame_counter=p1;
				}
				if(isBackbone(b) && isSidechain(a)) 
				{
					bb_sc_counter++;
					bb_sc_frame_counter=p1;
				}
				if(isSidechain(a) && isSidechain(b))
				{
					sc_sc_counter++;
					sc_sc_frame_counter+=p1;
				}	
				filename = "hb_"+printAtom1(a)+"_"+printAtom1(b)+".dat";
				statname = "hb_"+printAtom1(a)+"_"+printAtom1(b)+".stat";
				histname = "hb_"+printAtom1(a)+"_"+printAtom1(b)+".hist";
				acorname=  "hb_"+printAtom1(a)+"_"+printAtom1(b)+".cor";
				IntVector atom_pos;
				atom_pos.resize(2);
				atom_pos[0]=donor_id;
				atom_pos[1]=acceptor_id;
				myTex.addAtomList(atom_pos,p1);
				
				FILE *fout=fopen(filename.c_str(),"w");
				if(!fout) Error(WriteError(filename));
				int oneflag = X1.size() == 1;
				PrintLine(fout,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
				  printFrameHeader("Exist")+" "+
				  printDistanceHeader("Dist(D-A)")+" "+printString("Angle(D-H-A)",14)+" "+
				  printString("Dist(H-A)",14)+" "+((oneflag)?" ":printString("Atom",6)));
				for(int k=0;k<f.size();k++)
				{
					int exist;
					if(TotalDist[j][k]<=backbone_critical_distance && 
						TotalAngle[j][k]>=backbone_critical_angle) 
					/**/
					{
						exist=1;
					}
					/**/
					 else exist=0;
					countAcc[k]+=exist;
					Print(fout,printFrame(f[k])+" "+
						printFrame(exist)+" "+
						printDistance(TotalDist[j][k])+" "+
					  printNumber(TotalAngle[j][k],14,1)+" ");
					int hydro_id =st[i].HydroFlag[k];
					Data dmin=st[i].HydroAcc[j][k];
					PrintLine(fout,printNumber(dmin,14,3)+" "+
					((oneflag)?" ":printString(table[hydro_id-1].atom_name,6)));
					/**/
					/**/
				}
				fclose(fout);
				if(smooth_flag)
				{
				   string smoothname = "hb_"+printAtom1(a)+"_"+printAtom1(b)+".sda";
				   IntVector sf;
				   DoubleVector sd,sa;
				   makeSmooth(f,TotalDist[j],sf,sd,smooth_start,smooth_step);
				   makeAngleSmooth(f,TotalAngle[j],sf,sd,smooth_start,smooth_step);
				   PrintDistAndAngle(smoothname,sf,sd,sa);
				}
				string header1="#Distance statistics";
				string header2="#Angle statistics";
			PrintStatFile(statname,"w",f.size()/statcount,header1,f,TotalDist[j]);
			PrintAngleStatFile(statname,"a",f.size()/statcount,header2,f,TotalAngle[j]);
			string header3="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Hbcount")+" "+printFrameHeader("Dcount")+" "+
				printFrameHeader("Acount");
			PrintIntStatFile(statname,"a",f.size()/statcount,"#Hbond existance",header3,f,Acor);
				PrintAcor(acorname,f,Acor[0],kstep,astep);
				Data dmax,dmin,davg,dstd;
				getVectorStatistics(TotalDist[j],dmin,dmax,davg,dstd);
				Data dmax2,dmin2,davg2,dstd2;
				getAngleVectorStatistics(TotalAngle[j],dmin2,dmax2,davg2,dstd2);

				if(histflag)
				{
				string histnamed="hb_"+printAtom1(a)+"_"+printAtom1(b)+".histd";
				vector<HistStruct> st_dist;
				makeHist(TotalDist[j],st_dist,dmin,dmax,bindist);
				PrintHist(histnamed,st_dist);
				
				string histnamea="hb_"+printAtom1(a)+"_"+printAtom1(b)+".hista";
				vector<HistStruct> st_angle;
				makeHist(TotalAngle[j],st_angle,0.0,180.0,binangle);
				PrintAngleHist(histnamea,st_angle);
				
				string histname2="hb_"+printAtom1(a)+"_"+printAtom1(b)+".hist2";
				FILE *fp=fopen(histname2.c_str(),"w");
				if(!fp)  Error(WriteError(filename));
                                PrintLine(fp,"#"+printString("Dist",DISTANCE_WIDTH-1)+" "+
                                        printTorsionHeader("Angle")+" "+
                                        printFrameHeader("Common"));
                                for(int k=0;k<st_dist.size();k++)
                                {
                                        for(int m=0;m<st_angle.size();m++)
                                        {
                                        int icount=0;
                                        for(int l=0;l<f.size();l++)
                                        {
                                        if(TotalDist[j][l]>=st_dist[k].value-bindist/2.0
                                         && TotalDist[j][l]<=st_dist[k].value+bindist/2.0
                                         && TotalAngle[j][l]>=st_angle[m].value-binangle/2.0
                                         && TotalAngle[j][l]<=st_angle[m].value+binangle/2.0)
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
	for(int j=0;j<f.size();j++)
	myCounter.add(f[j],countAcc[j]?printAtom2(table[donor_id-1]):"",countAcc[j]);
	}
	myCounter.print("hb_count.dat");
	myCounter.printStat("hb_count.stat");
	if(smooth_flag)
		myCounter.printSmooth("hb_count.sda");
	if(histflag)
		myCounter.printHist("hb_count.hist");
	myTex.print();
	Log("\n#Statistics");
	Log("bb-bb "+printNumber(bb_bb_counter));
	Log("bb-sc "+printNumber(bb_sc_counter));
	Log("sc-sc "+printNumber(sc_sc_counter));
	int itotal=bb_bb_counter+bb_sc_counter+sc_sc_counter;
	Log("total "+printNumber(itotal));
}
