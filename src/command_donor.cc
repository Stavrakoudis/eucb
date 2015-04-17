# include <command_donor.h>
# include <global.h>
# include <donors.h>
# include <math.h>
# include <pdb.h>
# include <tex.h>
# include <counter.h>

CommandDonor::CommandDonor()
	:Command("DONOR")
{
	donorTeam = NULL;
	acceptorTeam = NULL;
}

void	CommandDonor::setTeams(Team *a,Team *b)
{
	donorTeam = a;
	acceptorTeam = b;
}

void	CommandDonor::setFlag(int Flag)
{
	donor_or_acceptor_first = Flag;
}

void	CommandDonor::acceptorFinder()
{
	vector<IntVector> CopyDonorCon;
	vector<IntVector> CopyAcceptorCon;
	vector<IntVector> BestAcceptorConDonor;
	vector<IntVector> BestAcceptorConHydro;
	IntVector ResAcceptor;
	vector<IntVector> ResDonor;
	vector<int> f;
	vector<int> CopyDonor;
	vector<int> CopyAcceptor;
	vector<IntVector> Flag;
	Counter myCounter;

	Tex	myTex("acceptor.tex");
	vector<string> s;
	s.resize(2);
	s[0]="Acceptor";
	s[1]="Donor";
	myTex.setHeader(s);
	myTex.setCaption("Hydrogen bonds");
	SelectDonorAcceptor(CopyDonor,CopyDonorCon,CopyAcceptor,CopyAcceptorCon,donorTeam,acceptorTeam);

	string donor_chain,acceptor_chain;
	dcd->setFrameStep(smart_skip);
	fflush(stdout);
	int Total=0;
	dcd->hbondsCheck(smart_distance,backbone_critical_angle,backbone_critical_percent,
		CopyDonor,CopyDonorCon,CopyAcceptor,Flag);
	for(int i=0;i<Flag.size();i++)
	{
		for(int j=0;j<Flag[i].size();j++)
		{
			if(Flag[i][j])
			{
				int pos=isin(ResAcceptor,j); //bazoume tin thesi tou donor kai oxi ayton.
				if(pos==-1)
				{
				int s = ResAcceptor.size();
				ResAcceptor.resize(s+1);		
				ResAcceptor[s]=j;

				ResDonor.resize(s+1);
				ResDonor[s].resize(1);
				ResDonor[s][0]=i;
				}
				else
				ResDonor[pos].push_back(i);
			}
		}
	}
	

	int atom_pos1,atom_pos2;
	vector<AcceptorStruct> st;
	dcd->setFrameStep(step);
	st.resize(ResAcceptor.size());
	for(int i=0;i<st.size();i++)
	{
		st[i].acceptor_id = CopyAcceptor[ResAcceptor[i]];
		st[i].Donor = ResDonor[i];	
		st[i].Hydro.resize(ResDonor[i].size());
		
		for(int j=0;j<ResDonor[i].size();j++) st[i].Hydro[j] = CopyDonorCon[ResDonor[i][j]];
	}
	dcd->acceptors(f,st,CopyDonor);

       int bb_bb_counter=0,bb_sc_counter=0,sc_sc_counter=0;


	for(int i=0;i<st.size();i++)
	{
		int null_counter=0;
		int acc_id = st[i].acceptor_id;
		int icount=0;
		Data pt;
		for(int j=0;j<f.size();j++)
		{
			for(int k=0;k<st[i].Donor.size();k++)
			{
				if(st[i].dist[k][j]<=backbone_critical_distance &&
				   st[i].angle[k][j]>=backbone_critical_angle)
				{
					icount++;
					break;
				}
			}
			pt=icount*1.0/f.size();
	//		if(pt>=backbone_critical_percent) break;
		}
		pt=icount*1.0/f.size();
		if(pt>=backbone_critical_percent)
		{
			atom a;
			atom_pos1=acc_id-1;
			a = table[acc_id-1];
			string filename="hba_"+printAtom1(a)+".dat";
			string statname="hba_"+printAtom1(a)+".stat";
			string histname="hba_"+printAtom1(a)+".hist";
			string acorname="hba_"+printAtom1(a)+".cor";
			FILE *fp=fopen(filename.c_str(),"w");
			if(!fp) Error(WriteError(filename));
			IntVector CountDonor;
			CountDonor.resize(st[i].Donor.size());
			for(int j=0;j<CountDonor.size();j++) CountDonor[j]=0;
			null_counter=0;
			int protein_counter=0;
			int water_counter=0;
			int sidechain_counter=0;
			int backbone_counter=0;
			DoubleVector MinDistance;
			MinDistance.resize(f.size());
			IntVector2 Acor;
			Acor.resize(5);
			Acor[0].resize(f.size());
			Acor[1].resize(f.size());
			Acor[2].resize(f.size());
			Acor[3].resize(f.size());
			Acor[4].resize(f.size());
			for(int j=0;j<f.size();j++)
			{
				MinDistance[j]=st[i].dist[0][j];
				Print(fp,printFrame(f[j])+" ");
				int donor_count=0;
				IntVector donor_found;
				donor_found.resize(0);
				int ifound=0;
				int water_found=0;
				int protein_found=0;
				int sidechain_found=0;
				int backbone_found=0;
				for(int k=0;k<st[i].Donor.size();k++)
				{
					if(st[i].dist[k][j]<=backbone_critical_distance &&
				  		 st[i].angle[k][j]>=backbone_critical_angle)
					{
						
						if(st[i].dist[k][j]<MinDistance[j])
							MinDistance[j]=st[i].dist[k][j];
						ifound=1;
						donor_count++;
						donor_found.push_back(k);
						CountDonor[k]++;
					int pos=st[i].Donor[k];
					int acc_id = CopyDonor[pos];
				
					if(!protein_found && isProteinDonor(acc_id,table))
					{
						protein_found=1;
						protein_counter++;
					}
					if(!water_found && isWaterDonor(acc_id,table))
					{
						water_found=1;
						water_counter++;
					}
					if(!sidechain_found && isSideChainDonor(acc_id,CopyDonorCon[pos],table))
					{
						sidechain_found=1;
						sidechain_counter++;
					}
					if(!backbone_found && isBackBoneDonor(acc_id,CopyDonorCon[pos],table))
					{
						backbone_found=1;
						backbone_counter++;
					}
					}
				}
				if(!ifound) {null_counter++; Acor[0][j]=0;} else Acor[0][j]=1;
				if(protein_found) Acor[1][j]=1; else Acor[1][j]=0;
				if(backbone_found) Acor[2][j]=1;else Acor[2][j]=0;
				if(sidechain_found) Acor[3][j]=1; else Acor[3][j]=0;
				if(water_found) Acor[4][j]=1; else Acor[4][j]=0;


				myCounter.add(f[j],donor_count?printAtom2(table[acc_id-1]):"",donor_count);
				Print(fp,printFrame(donor_count)+" ");
				for(int k=0;k<donor_count;k++)
				{
					int donor_id = CopyDonor[st[i].Donor[donor_found[k]]];
					atom b = table[donor_id-1];
					/**/
					int s1=printAtom4(a).size();
                                int s2=printAtom4(b).size();
                                string ss=printAtom4(a);
                                int kk;
                                for(kk=0;kk<18-s1;kk++)
                                        ss+=" ";
                                ss+=printAtom4(b);
                                for( kk=0;kk<18-s2;kk++)
                                        ss+=" ";
				Data p1=0.0;
				for(int l1=0;l1<f.size();l1++)
				if(st[i].dist[k][l1]<=backbone_critical_distance &&
				  		 st[i].angle[k][l1]>=backbone_critical_angle)
					p1=p1+1;
				p1=p1/f.size();
                                Log(ss+" "+printPercent(p1));
                                if(isBackbone(a) && isBackbone(b))
                                {
                                        bb_bb_counter++;
                                }
                                if(isBackbone(a) && isSidechain(b))
                                {
                                        bb_sc_counter++;
                                }
				  if(isBackbone(b) && isSidechain(a))
                                {
                                        bb_sc_counter++;
                                }
                                if(isSidechain(a) && isSidechain(b))
                                {
                                        sc_sc_counter++;
                                }

					/**/
					Print(fp,printString(b.chain_id,5)+" "+printNumber(b.res_id,5)+" "+
						printString(b.atom_name,5)+" ");
					Print(fp,printDistance(st[i].dist[donor_found[k]][j])+""+
						printAngle(st[i].angle[donor_found[k]][j])+" ");
				}
				PrintLine(fp,"");
			}
			fclose(fp);
                        PrintAcor(acorname,f,Acor[0],kstep,astep);
			string header1="#Existence statistics";
			string header2="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Total")+" "+printFrameHeader("Protein")+" "+
				printFrameHeader("Backbone")+" "+printFrameHeader("Sidechain")+" "+
				printFrameHeader("Water");
			PrintIntStatFile(statname,"w",f.size()/statcount,header1,header2,f,Acor);


			fp=fopen(statname.c_str(),"a");
			if(!fp) Error(WriteError(statname));
			PrintLine(fp,"#Atom statistics");
			for(int j=0;j<st[i].Donor.size();j++)
			{
				if(CountDonor[j])
				{
					int donor_id = CopyDonor[st[i].Donor[j]];
					atom_pos2=donor_id-1;
					atom b = table[donor_id-1];
					Print(fp,printString(b.chain_id,5)+" "+printNumber(b.res_id,5)+" "+
						printString(b.atom_name,5)+" ");
					PrintLine(fp,printFrame(CountDonor[j]));

					Data p=CountDonor[j]*1.0/f.size();
					IntVector atom_pos;
					atom_pos.resize(2);
					atom_pos[0]=atom_pos1+1;
					atom_pos[1]=atom_pos2+1;
					myTex.addAtomList(atom_pos,p);
				}
			}
			fclose(fp);
			if(histflag)
			{
				vector<HistStruct> st_dist;
				Data dmin,dmax,davg,dstd;
				
				getVectorStatistics(MinDistance,dmin,dmax,davg,dstd);
				makeHist(MinDistance,st_dist,dmin,dmax,bindist);
				FILE *fout = fopen(histname.c_str(),"w");
				if(!fout) Error(WriteError(histname));
				for(int l=0;l<st_dist.size();l++)
					PrintLine(fout,printNumber(st_dist[l].value,8,2)+" "+
						printPercent(st_dist[l].percent)+" "+
						printNumber(st_dist[l].count,5));
				fclose(fout);
			}
			//OMADOPOIHSH
			for(int j=0;j<st[i].Donor.size();j++)
			{
				int donor_id = CopyDonor[st[i].Donor[j]];
				atom donor_atom=table[donor_id-1];
				if(CountDonor[j] && donor_atom.res_name=="ARG" && donor_atom.atom_name=="NH1")
				{
					for(int k=0;k<st[i].Donor.size();k++)
					{
						int donor2_id = CopyDonor[st[i].Donor[k]];
						atom donor2_atom=table[donor_id-1];
						if(CountDonor[k] && donor_atom.res_name=="ARG" && donor_atom.atom_name=="NH2")
						{
							int common=0;
							for(int l=0;l<f.size();l++)
							{
								if((st[i].dist[j][l]<=backbone_critical_distance && st[i].angle[j][l]>=backbone_critical_angle) || (st[i].dist[k][l]<=backbone_critical_distance && st[i].angle[k][l]>=backbone_critical_angle)) common++;
							}
							break;
						}
					}
				} 
			}

			//END OMADOPOIHSH
		}
	}
	myTex.print();
	Log("\n#Statistics");
        Log("bb-bb "+printNumber(bb_bb_counter));
        Log("bb-sc "+printNumber(bb_sc_counter));
        Log("sc-sc "+printNumber(sc_sc_counter));
        int itotal=bb_bb_counter+bb_sc_counter+sc_sc_counter;
        Log("total "+printNumber(itotal));

	myCounter.print("hba_count.dat");
	myCounter.printStat("hba_count.stat");
	if(smooth_flag)
		myCounter.printSmooth("hba_count.sda");
	if(histflag)
		myCounter.printHist("hba_count.hist");
}

void	CommandDonor::donorFinder()
{
	vector<IntVector> CopyDonorCon;
	vector<IntVector> CopyAcceptorCon;
	vector<int> f;
	vector<int> CopyDonor;
	vector<int> CopyAcceptor;
	string donor_chain,acceptor_chain;
	Counter	myCounter;

	SelectDonorAcceptor(CopyDonor,CopyDonorCon,CopyAcceptor,
		CopyAcceptorCon,donorTeam,
		acceptorTeam);
	dcd->setFrameStep(smart_skip);
	IntVector X1;
	IntVector X2;
	IntVector HydroFlag;
	IntVector	  ResDonor;
	vector<IntVector> ResAcceptor;
	DoubleVector2 	TotalDist,TotalAngle;
	vector<DoubleVector> TotalDist2;
	vector<IntVector> Flag;
	int atom_pos1,atom_pos2;

	dcd->hbondsCheck(smart_distance,backbone_critical_angle,backbone_critical_percent,
		CopyDonor,CopyDonorCon,CopyAcceptor,Flag);
	vector<HbondStruct> st;
	
	Tex	myTex("donor.tex");
	vector<string> s;
	s.resize(2);
	s[0]="Donor";
	s[1]="Acceptor";
	myTex.setHeader(s);
	myTex.setCaption("Hydrogen bonds");
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
		st[i].Acceptor=ResAcceptor[i];
		st[i].Hydro=CopyDonorCon[ResDonor[i]];
	}
	dcd->hbonds(f,st);
	DoubleVector MinDistance;

	
        int bb_bb_counter=0,bb_sc_counter=0,sc_sc_counter=0;
        int bb_bb_frame_counter=0,bb_sc_frame_counter=0,sc_sc_frame_counter=0;

	Log("Critical distance = "+printDistance(backbone_critical_distance));
	Log("Critical angle  = "+printAngle(backbone_critical_angle));
	for(int i=0;i<ResDonor.size();i++)
	{
		int donor_id =st[i].donor_id;
		X1=CopyDonorCon[ResDonor[i]];
		X2=ResAcceptor[i];
		TotalDist=st[i].dist;
		IntVector X2Counter;
		X2Counter.resize(X2.size());
		for(int j=0;j<X2Counter.size();j++) X2Counter[j]=0;

		MinDistance.resize(f.size());
		for(int j=0;j<MinDistance.size();j++) MinDistance[j]=1e+100;

		TotalAngle=st[i].angle;
		TotalDist2=st[i].HydroAcc;
		HydroFlag = st[i].HydroFlag;
		int count1=0;
		int nullcounter=0;
		int total=0;
		int protein_counter=0;
		int water_counter=0;
		int sidechain_counter=0;
		int backbone_counter=0;
		Data pt;
	
		for(int j=0;j<f.size();j++)
		{
			for(int k=0;k<X2.size();k++)
			{
				if(TotalDist[k][j]<=backbone_critical_distance && 
			           TotalAngle[k][j]>=backbone_critical_angle)
				{
					count1++;
					break;
				}
			}
			pt=count1*1.0/f.size();
			if(pt>=backbone_critical_percent) break;
		}
		pt=count1*1.0/f.size();
		if(pt>=backbone_critical_percent)
		{
			atom a;
			a = table[donor_id-1];
			atom_pos1=donor_id-1;
			string filename="hbd_"+printAtom1(a)+".dat";
			string statname="hbd_"+printAtom1(a)+".stat";
			string histname="hbd_"+printAtom1(a)+".hist";
			string acorname="hbd_"+printAtom1(a)+".cor";
			FILE *fp=fopen(filename.c_str(),"w");
			if(!fp) Error(WriteError(filename));
			IntVector2 Acor;
			Acor.resize(5);
			Acor[0].resize(f.size());
			Acor[1].resize(f.size());
			Acor[2].resize(f.size());
			Acor[3].resize(f.size());
			Acor[4].resize(f.size());
		
			for(int j=0;j<f.size();j++)
			{
				vector<int> acc;
				acc.resize(0);
				int count_acc = 0;
				Print(fp,printFrame(f[j])+" ");
				int found=0;
				int water_found=0;
				int protein_found=0;
				int sidechain_found=0;
				int backbone_found=0;
				MinDistance[j]=TotalDist[0][j];
				for(int k=0;k<X2.size();k++)
				{
					if(TotalDist[k][j]<=backbone_critical_distance && 
					   TotalAngle[k][j]>=backbone_critical_angle)
					{
						count_acc++;
						acc.push_back(k);
						X2Counter[k]++;
						found=1;
						int acc_id = X2[k];
						int pos = isin(CopyAcceptor,acc_id);
						if(!protein_found && isProteinAcceptor(acc_id,table))
						{
							protein_found=1;
							protein_counter++;
						}
						if(!water_found && isWaterAcceptor(acc_id,table))
						{
							water_counter++;
							water_found=1;
						}
						if(!sidechain_found && isSideChainAcceptor(acc_id,CopyAcceptorCon[pos],table))
						{
							sidechain_found=1;
							sidechain_counter++;
						}
						if(!backbone_found && isBackBoneAcceptor(acc_id,CopyAcceptorCon[pos],table))
						{
							backbone_found=1;
							backbone_counter++;
						}
					}
					if(TotalDist[k][j]<=MinDistance[j])
						MinDistance[j]=TotalDist[k][j];
				}
				if(!found) {nullcounter++; Acor[0][j]=0;} else Acor[0][j]=1;
				if(protein_found) Acor[1][j]=1; else Acor[1][j]=0;
				if(backbone_found) Acor[2][j]=1;else Acor[2][j]=0;
				if(sidechain_found) Acor[3][j]=1; else Acor[3][j]=0;
				if(water_found) Acor[4][j]=1; else Acor[4][j]=0;

				myCounter.add(f[j],count_acc?printAtom2(table[donor_id-1]):"",count_acc);

				Print(fp,printNumber(count_acc,5)+" ");
				for(int k=0;k<count_acc;k++)
				{
					atom b = table[X2[acc[k]]-1];
					Print(fp,printString(b.chain_id,5)+" "+printNumber(b.res_id,5)+" "+
						printString(b.atom_name,5)+" ");
					Print(fp,printDistance(TotalDist[acc[k]][j])+" "+
						printAngle(TotalAngle[acc[k]][j])+" ");
					int hydro_id =st[i].HydroFlag[j];
					Data dmin=st[i].HydroAcc[acc[k]][j];
					Print(fp,printDistance(dmin)+" "+table[hydro_id-1].atom_name);
				}
				PrintLine(fp,"");
			}
			fclose(fp);
                        PrintAcor(acorname,f,Acor[0],kstep,astep);
			string header1="#Existence statistics";
			string header2="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Total")+" "+printFrameHeader("Protein")+" "+
				printFrameHeader("Backbone")+" "+printFrameHeader("Sidechain")+" "+
				printFrameHeader("Water");
			PrintIntStatFile(statname,"w",f.size()/statcount,header1,header2,f,Acor);

		
			fp=fopen(statname.c_str(),"a");
			if(!fp) Error(WriteError(statname));
			PrintLine(fp,"#Atom statistics");
			for(int l=0;l<X2Counter.size();l++)
			{
				if(X2Counter[l])
				{
					int acc_id = X2[l];
					atom_pos2=X2[l]-1;
					atom b = table[X2[l]-1];
					/**/
					  int s1=printAtom4(a).size();
                                int s2=printAtom4(b).size();
                                string ss=printAtom4(a);
                                int kk;
                                for(kk=0;kk<18-s1;kk++)
                                        ss+=" ";
                                ss+=printAtom4(b);
                                for( kk=0;kk<18-s2;kk++)
                                        ss+=" ";
				Data p1=0.0;
				for(int l1=0;l1<f.size();l1++)
				{
					if(TotalDist[l][l1]<=backbone_critical_distance && 
					   TotalAngle[l][l1]>=backbone_critical_angle)
						p1=p1+1;
				}
				p1=p1/f.size();
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

					/**/
					Print(fp,printString(b.chain_id,5)+" "+printNumber(b.res_id,5)+" "+
						printString(b.atom_name,5)+" ");
					PrintLine(fp,printFrame(X2Counter[l]));
					Data dist=0.0;	
					Data p=X2Counter[l]*1.0/f.size();
					IntVector atom_pos;
					atom_pos.resize(2);
					atom_pos[0]=atom_pos1+1;
					atom_pos[1]=atom_pos2+1;
					myTex.addAtomList(atom_pos,p);
				}
			}

			if(histflag)
			{
			vector<HistStruct> st_dist;
			Data dmin,dmax,davg,dstd;
			getVectorStatistics(MinDistance,dmin,dmax,davg,dstd);
			makeHist(MinDistance,st_dist,dmin,dmax,bindist);
			FILE *fout = fopen(histname.c_str(),"w");
			if(!fout) Error(WriteError(histname));
			for(int l=0;l<st_dist.size();l++)
				PrintLine(fout,printNumber(st_dist[l].value,8,2)+" "+
					printPercent(st_dist[l].percent)+" "+
					printNumber(st_dist[l].count,5));
			fclose(fout);
			
			}

			//ENOPOIHSH
			vector<TexGroup> Group;	
			for(int l=0;l<X2Counter.size();l++)
			{
				int acc_id = X2[l];
				atom acceptor1=table[acc_id-1];
				if(X2Counter[l] && acceptor1.res_name == "GLU"  && 
					acceptor1.atom_name=="OE1")
				{
					for(int m=0;m<X2Counter.size();m++)
					{
						int acc2_id=X2[m];
						atom acceptor2=table[acc2_id-1];
						if(X2Counter[m] && acceptor2.res_name=="GLU" && 
							acceptor2.atom_name=="OE2" && 
							acceptor1.res_id == acceptor2.res_id &&
							acceptor1.chain_id== acceptor2.chain_id)
						{
							TexGroup p;
							p.first=acc_id;
							p.second=acc2_id;
							int common=0;
							for(int ll=0;ll<f.size();ll++)
							{
								if((TotalDist[l][ll]<=backbone_critical_distance && 
					   			TotalAngle[l][ll]>=backbone_critical_angle) ||
								(TotalDist[m][ll]<=backbone_critical_distance && 
					   			TotalAngle[m][ll]>=backbone_critical_angle))
									common++;
							}
							p.percent=common*1.0/f.size();
							Group.push_back(p);
							break;
						}
					}
				}

				if(X2Counter[l] && acceptor1.res_name == "ASP"  && 
					acceptor1.atom_name=="OD1")
				{
					for(int m=0;m<X2Counter.size();m++)
					{
						int acc2_id=X2[m];
						atom acceptor2=table[acc2_id-1];
						if(X2Counter[m] && acceptor2.res_name=="ASP" && 
							acceptor2.atom_name=="OD2" && 
							acceptor1.res_id == acceptor2.res_id &&
							acceptor1.chain_id== acceptor2.chain_id)
						{
							TexGroup p;
							p.first=acc_id;
							p.second=acc2_id;
							int common=0;
							for(int ll=0;ll<f.size();ll++)
							{
								if((TotalDist[l][ll]<=backbone_critical_distance && 
					   			TotalAngle[l][ll]>=backbone_critical_angle) ||
								(TotalDist[m][ll]<=backbone_critical_distance && 
					   			TotalAngle[m][ll]>=backbone_critical_angle))
									common++;
							}
							p.percent=common*1.0/f.size();
							Group.push_back(p);
							break;
						}
					}
				}
				
			}
			
			//END OF ENOPOIHSH
			myTex.addGroupList(Group);
		}
	}
	myTex.print();

	 Log("\n#Statistics");
        Log("bb-bb "+printNumber(bb_bb_counter));
        Log("bb-sc "+printNumber(bb_sc_counter));
        Log("sc-sc "+printNumber(sc_sc_counter));
        int itotal=bb_bb_counter+bb_sc_counter+sc_sc_counter;
        Log("total "+printNumber(itotal));

	myCounter.print("hbd_count.dat");
	myCounter.printStat("hbd_count.stat");
	if(smooth_flag)
		myCounter.printSmooth("hbd_count.sda");
	if(histflag)
		myCounter.printHist("hbd_count.hist");
}


void	CommandDonor::Run()
{
	if(donor_or_acceptor_first==1)
	{
		donorFinder();
	}
	else 
	{
		acceptorFinder();
		fclose(fplog);
		rename("DONOR.log","ACCEPTOR.log");
		fplog=fopen("ACCEPTOR.log","a");
	}
}
