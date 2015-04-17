# include <global.h>
# include <command_contact.h>
# include <pdb.h>
# include <donors.h>

CommandContact::CommandContact()
	:Command("CONTACT")
{
	selection1 = "all";
	selection2 = "all";
}

void	CommandContact::setSelection1(string s)
{
	selection1 = s;
}

void	CommandContact::setSelection2(string s)
{	
	selection2 = s;
}

typedef struct
{
	string chain_id1;
	string chain_id2;
	int    res_id1;
	int    res_id2;
	string  res_name1;
	string  res_name2;
	int 	vdw_count;
	int     hb_count;
	int     salt_count;
}ContactPair;

typedef struct 
{
	string chain_id;
	int    res_id;
	string res_name;
	int 	pdb_vdw_count;
	int     pdb_hb_count;
	int     pdb_salt_count;
	int 	dcd_vdw_count;
	int     dcd_hb_count;
	int     dcd_salt_count;
}SynopsisStruct;

void CommandContact::Run()
{
	IntVector atom1;
	IntVector atom2;
	if(pdb_name=="") 
	{
		Log("You must specify the pdb file with the -pdb option");
		psf_error("You must specify the pdb file with the -pdb option");
	}
	Team *contact_team1=team1;
	Team *contact_team2=team2;
	if(contact_team1 == NULL || contact_team2==NULL)  Error(SeqError("contact"));
	for(int i=0;i<table.size();i++)
	{
		if(table[i].atom_name[0]=='H') continue;
		string chain;
		if(contact_team1==NULL || contact_team1->find(table[i],chain))
		{
			if(selection1=="all") atom1.push_back(i);
			else
			if(selection1=="ca" && table[i].atom_name=="CA") atom1.push_back(i);
			else
			if(selection1=="backbone3" && isBackbone3(table[i])) atom1.push_back(i);
			else
			if(selection1=="backbone4" && isBackbone4(table[i])) atom1.push_back(i);
			else
			if(selection1=="backbone" && isBackbone(table[i])) atom1.push_back(i);
			else
			if(selection1=="sidechain" && isBackbone(table[i])) atom1.push_back(i);
			else
			{
				Log("Wrong selection for sidechain "+selection1);
				psf_error("Wrong selection for sidechain "+selection1);
			}
		}
		if(contact_team2==NULL || contact_team2->find(table[i],chain))
		{
			if(selection2=="all") atom2.push_back(i);
			else
			if(selection2=="ca" && table[i].atom_name=="CA") atom2.push_back(i);
			else
			if(selection2=="backbone3" && isBackbone3(table[i])) atom2.push_back(i);
			else
			if(selection2=="backbone4" && isBackbone4(table[i])) atom2.push_back(i);
			else
			if(selection2=="backbone" && isBackbone(table[i])) atom2.push_back(i);
			else
			if(selection2=="sidechain" && isBackbone(table[i])) atom2.push_back(i);
			else
			{
				Log("Wrong selection for sidechain "+selection2);
				psf_error("Wrong selection for sidechain "+selection2);
			}
		}
	}
	Log("#Critical percent  = "+printPercent(contact_percent));
	Log("#Critical distance = "+printDistance(contact_distance));
	IntVector F;
	vector<NoeStruct> noe;
	for(int i=0;i<atom1.size();i++)
	{
		for(int j=0;j<atom2.size();j++)
		{
			if(atom1[i]==atom2[j]) continue;
			if(table[atom1[i]].chain_id == table[atom2[j]].chain_id &&
			 table[atom1[i]].res_id == table[atom2[j]].res_id) continue;
			if(table[atom1[i]].chain_id == table[atom2[j]].chain_id &&
				abs(table[atom1[i]].res_id - table[atom2[j]].res_id)<=res_diff)
			{
				continue;
			}
			Data pdb_dist =getPointDistance(pdbpos[atom1[i]],pdbpos[atom2[j]]);
			if(pdb_dist>=2.0*contact_distance) continue;
			NoeStruct n;
			n.atom1=atom1[i]+1;
			n.atom2=atom2[j]+1;
			noe.push_back(n);
		}
	}
	dcd->getNoe(F,noe);
	IntVector valid;
	for(int i=0;i<noe.size();i++)
	{
		Data dist_pdb=getPointDistance(pdbpos[noe[i].atom1-1],pdbpos[noe[i].atom2-1]);
		if(dist_pdb<=contact_distance) 
		{
			valid.push_back(i);
			continue;
		}
		else
		{
			Data p=getBelowPercent(noe[i].D,contact_distance);
			if(p>=contact_percent) valid.push_back(i);
			continue;
		}
		noe[i].D.resize(0);	
	}
	Log("#Total valid contact pairs = "+printNumber((int)valid.size()));
	IntVector DonorAcceptor;
	IntVector AcceptorDonor;
	for(int i=0;i<valid.size();i++)
	{
		int isdonor_atom1=0;
		int isdonor_atom2=0;
		int isacceptor_atom1=0;
		int isacceptor_atom2=0;
		for(int j=0;j<Donor.size();j++)
		{
			if(Donor[j]==noe[valid[i]].atom1) isdonor_atom1=1;	
			else
			if(Donor[j]==noe[valid[i]].atom2) isdonor_atom2=1;
		}
		for(int j=0;j<Acceptor.size();j++)
		{
			if(Acceptor[j]==noe[valid[i]].atom1) isacceptor_atom1=1;
			else
			if(Acceptor[j]==noe[valid[i]].atom2) isacceptor_atom2=1;
		}
		if(isdonor_atom1 && isacceptor_atom2) DonorAcceptor.push_back(valid[i]);
		if(isacceptor_atom1 && isdonor_atom2) AcceptorDonor.push_back(valid[i]);
	}
	Log("#Donor_Acceptor pairs = "+printNumber((int)DonorAcceptor.size()));
	Log("#Acceptor_Donor pairs = "+printNumber((int)AcceptorDonor.size()));
	for(int i=0;i<valid.size();i++)
	{
		IntVector Acor;
		Acor.resize(noe[valid[i]].D.size());
		for(int j=0;j<noe[valid[i]].D.size();j++)
		{
			if(noe[valid[i]].D[j]<=contact_distance) Acor[j]=1; else Acor[j]=0;
		}
		Data p=getBelowPercent(noe[valid[i]].D,contact_distance);
		atom a,b;
		a=table[noe[valid[i]].atom1-1];
		b=table[noe[valid[i]].atom2-1];
		string filename="contact_"+printAtom1(a)+"_"+printAtom1(b)+".dat";
		string statname="contact_"+printAtom1(a)+"_"+printAtom1(b)+".stat";
		string smoothname="contact_"+printAtom1(a)+"_"+printAtom1(b)+".sda";
		string histname="contact_"+printAtom1(a)+"_"+printAtom1(b)+".hist";
		string acorname="contact_"+printAtom1(a)+"_"+printAtom1(b)+".cor";
                PrintAcor(acorname,F,Acor,kstep,astep);

		Data dist_pdb=getPointDistance(pdbpos[noe[valid[i]].atom1-1],
			pdbpos[noe[valid[i]].atom2-1]);
		PrintDistances(filename,F,noe[valid[i]].D);
		FILE *fp;
		fp=fopen(statname.c_str(),"w");
		if(!fp)  Error(WriteError(statname));
		Data dmin,dmax,davg,dstd;
		getVectorStatistics(noe[valid[i]].D,dmin,dmax,davg,dstd);
		string header1="#Distance statistics";
		PrintStatFile(statname,"w",F.size()/statcount,header1,F,noe[valid[i]].D);
		if(smooth_flag)
		{
			IntVector sf;
			DoubleVector sd;
			makeSmooth(F,noe[valid[i]].D,sf,sd,smooth_start,smooth_step);
			PrintDistances(smoothname,sf,sd);
		}
		if(histflag)
		{
			vector<HistStruct> st;
			makeHist(noe[valid[i]].D,st,dmin,dmax,bindist);
			PrintHist(histname,st);
		}
	}
	vector<HbondStruct> st;
	st.resize(DonorAcceptor.size()+AcceptorDonor.size());
	for(int i=0;i<DonorAcceptor.size();i++)
	{
		int donor_pos;
		IntVector A;
		A.resize(1);A[0]=noe[DonorAcceptor[i]].atom2;
		st[i].donor_id=noe[DonorAcceptor[i]].atom1;
		for(int j=0;j<Donor.size();j++) if(Donor[j]==st[i].donor_id) {donor_pos=j;break;}
		st[i].Acceptor=A;
		st[i].Hydro=DonorCon[donor_pos];
	}
	for(int i=0;i<AcceptorDonor.size();i++)
	{
		/*
		int donor_pos;
		IntVector A;
		A.resize(1);A[0]=noe[DonorAcceptor[i]].atom1;
		st[DonorAcceptor.size()+i].donor_id=noe[AcceptorDonor[i]].atom2;
		for(int j=0;j<Donor.size();j++) 
		if(Donor[j]==st[DonorAcceptor.size()+i].donor_id) {donor_pos=j;break;}
		st[DonorAcceptor.size()+i].Acceptor=A;
		st[DonorAcceptor.size()+i].Hydro=DonorCon[donor_pos];
		*/
		int donor_pos;
		IntVector A;
		A.resize(1);A[0]=noe[AcceptorDonor[i]].atom1;
		st[DonorAcceptor.size()+i].donor_id=noe[AcceptorDonor[i]].atom2;
		for(int j=0;j<Donor.size();j++) 
		if(Donor[j]==st[DonorAcceptor.size()+i].donor_id) {donor_pos=j;break;}
		st[DonorAcceptor.size()+i].Acceptor=A;
		st[DonorAcceptor.size()+i].Hydro=DonorCon[donor_pos];
	}
	dcd->hbonds(F,st);
	IntVector DistanceType;
	//Distance types
	//0-VDW
	//1-SALT
	//2-HB
	DistanceType.resize(valid.size());
	for(int i=0;i<valid.size();i++) 
		DistanceType[i]=0;
	for(int i=0;i<st.size();i++)
	{
		DoubleVector2 TotalDist;
		DoubleVector2 TotalAngle;
		TotalDist=st[i].HydroAcc;
		TotalAngle=st[i].angle;
		int count=0;
		for(int j=0;j<F.size();j++)
		{
			if(TotalDist[0][j]>=backbone_critical_distance && 
			TotalAngle[0][j]>=backbone_critical_angle) count++;
		}
		if(count>=backbone_critical_percent)
		{
			for(int j=0;j<valid.size();j++)
			{
				if(noe[valid[j]].atom1==st[i].donor_id && 
					noe[valid[j]].atom2==st[i].Acceptor[0])
				{
					DistanceType[j]=2;
					break;
				}
				if(noe[valid[j]].atom2==st[i].donor_id &&
					noe[valid[j]].atom1==st[i].Acceptor[0])
				{
					DistanceType[j]=2;
					break;
				}
			}
		}
	}

	Log(printString("#Atom1",-15)+" "+printString("Atom2",-15)+" "+
			printDistanceHeader("PDBdist")+" "+
			printString("Type1",-5)+" "+printString("Type2",-5)+" "+
			printPercentHeader("%")+" "+
			printDistanceHeader("Min")+" "+printDistanceHeader("Max")+" "+
			printDistanceHeader("Avg")+" "+printDistanceHeader("Std"));
	vector<ContactPair> pdb_pair;
	vector<ContactPair> dcd_pair;
	for(int i=0;i<valid.size();i++)
	{
		atom a,b;
		a=table[noe[valid[i]].atom1-1];
		b=table[noe[valid[i]].atom2-1];
		Data dmin,dmax,davg,dstd,p;
		p=getBelowPercent(noe[valid[i]].D,contact_distance);
		Data dist_pdb=getPointDistance(pdbpos[noe[valid[i]].atom1-1],
			pdbpos[noe[valid[i]].atom2-1]);
		string type1="VDW";
		if(dist_pdb>contact_distance) type1="-";
		string type2="VDW";
		if(p<contact_percent) type2="-";
		int pos_atom1=0;
		int pos_atom2=0;
		int neg_atom1=0;
		int neg_atom2=0;
		if(DistanceType[i]==2) type2="HB";
		else
		if(type2!="-")
		{
		  if(isPositive(a)) pos_atom1=1;
		  if(isPositive(b)) pos_atom2=1;
		  if(isNegative(a)) neg_atom1=1;
		  if(isNegative(b)) neg_atom2=1;
		  if(pos_atom1 && neg_atom2) type2="SALT";
		  if(pos_atom2 && neg_atom1) type2="SALT";
		}
		int pos1=isin(DonorAcceptor,valid[i]);
		int pos2=isin(AcceptorDonor,valid[i]);
		if((pos1!=-1 || pos2!=-1) && dist_pdb<=backbone_critical_distance) type1="HB";
		else
		if(type1!="-")
		{
		if(pos_atom1 && neg_atom2) type1="SALT";
		if(pos_atom2 && neg_atom1) type1="SALT";
		}
		
		getVectorStatistics(noe[valid[i]].D,dmin,dmax,davg,dstd);
		Log(printAtom4(a)+" "+printAtom4(b)+" "+printDistance(dist_pdb)+" "+
			printString(type1,-5)+" "+printString(type2,-5)+" "+printPercent(p)+" "+
		printDistance(dmin)+" "+printDistance(dmax)+" "+printDistance(davg)+" "+
		printDistance(dstd));

		int ipos=-1;
		for(int j=0;j<pdb_pair.size();j++)
		{
			if(pdb_pair[j].res_name1==a.res_name && pdb_pair[j].chain_id1==a.chain_id
			  && pdb_pair[j].res_id1 == a.res_id && pdb_pair[j].res_name2==b.res_name
			  && pdb_pair[j].res_id2 == b.res_id && pdb_pair[j].chain_id2==b.chain_id)
			{
				ipos=j;
				break;
			}
			if(pdb_pair[j].res_name1==b.res_name && pdb_pair[j].chain_id1==b.chain_id
			  && pdb_pair[j].res_id1 == b.res_id && pdb_pair[j].res_name2==a.res_name
			  && pdb_pair[j].res_id2 == a.res_id && pdb_pair[j].chain_id2==a.chain_id)
			{
				ipos=j;
				break;
			}
		}
		if(ipos==-1)
		{
			int s=pdb_pair.size();
			pdb_pair.resize(s+1);
			pdb_pair[s].res_name1=a.res_name;
			pdb_pair[s].chain_id1=a.chain_id;
			pdb_pair[s].res_id1=a.res_id;	
			pdb_pair[s].res_name2=b.res_name;
			pdb_pair[s].chain_id2=b.chain_id;
			pdb_pair[s].res_id2=b.res_id;	
			pdb_pair[s].vdw_count=pdb_pair[s].hb_count=pdb_pair[s].salt_count=0;
			ipos=s;
		}
		if(type1=="VDW") pdb_pair[ipos].vdw_count++;
		if(type1=="HB") pdb_pair[ipos].hb_count++;
		if(type1=="SALT") pdb_pair[ipos].salt_count++;

		ipos=-1;
		for(int j=0;j<dcd_pair.size();j++)
		{
			if(dcd_pair[j].res_name1==a.res_name && dcd_pair[j].chain_id1==a.chain_id
			  && dcd_pair[j].res_id1 == a.res_id && dcd_pair[j].res_name2==b.res_name
			  && dcd_pair[j].res_id2 == b.res_id && dcd_pair[j].chain_id2==b.chain_id)
			{
				ipos=j;
				break;
			}
			if(dcd_pair[j].res_name1==b.res_name && dcd_pair[j].chain_id1==b.chain_id
			  && dcd_pair[j].res_id1 == b.res_id && dcd_pair[j].res_name2==a.res_name
			  && dcd_pair[j].res_id2 == a.res_id && dcd_pair[j].chain_id2==a.chain_id)
			{
				ipos=j;
				break;
			}
		}
		if(ipos==-1)
		{
			int s=dcd_pair.size();
			dcd_pair.resize(s+1);
			dcd_pair[s].res_name1=a.res_name;
			dcd_pair[s].chain_id1=a.chain_id;
			dcd_pair[s].res_id1=a.res_id;	
			dcd_pair[s].res_name2=b.res_name;
			dcd_pair[s].chain_id2=b.chain_id;
			dcd_pair[s].res_id2=b.res_id;	
			dcd_pair[s].vdw_count=dcd_pair[s].hb_count=dcd_pair[s].salt_count=0;
			ipos=s;
		}
		if(type2=="VDW") dcd_pair[ipos].vdw_count++;
		if(type2=="HB") dcd_pair[ipos].hb_count++;
		if(type2=="SALT") dcd_pair[ipos].salt_count++;

	}
	Log("#STATISTICS");
	Log(printString("#Res1",-15)+" "+printString("Res2",-15)+" "+
		printString("VDW1",5)+" "+
		printString("HB1",5)+" "+printString("SALT1",5)+" "+
		printString("VDW2",5)+" "+
		printString("HB2",5)+" "+printString("SALT2",5));
	vector<SynopsisStruct> syn1;
	vector<SynopsisStruct> syn2;
	for(int i=0;i<pdb_pair.size();i++)
	{
		char s1[100];
		char s2[100];
		sprintf(s1,"%s:%d:%s",pdb_pair[i].chain_id1.c_str(),
			pdb_pair[i].res_id1,pdb_pair[i].res_name1.c_str());
		sprintf(s2,"%s:%d:%s",pdb_pair[i].chain_id2.c_str(),
			pdb_pair[i].res_id2,pdb_pair[i].res_name2.c_str());
		Log(printString(s1,-15)+" "+printString(s2,-15)+" "+
		printNumber(pdb_pair[i].vdw_count,5)+" "+
		printNumber(pdb_pair[i].hb_count,5)+" "+
		printNumber(pdb_pair[i].salt_count,5)+" "+
		printNumber(dcd_pair[i].vdw_count,5)+" "+
		printNumber(dcd_pair[i].hb_count,5)+" "+
		printNumber(dcd_pair[i].salt_count,5));

		int ipos=-1;
		for(int j=0;j<syn1.size();j++)
		{
			if(syn1[j].chain_id == pdb_pair[i].chain_id1 && 
			syn1[j].res_id  == pdb_pair[i].res_id1 &&
			syn1[j].res_name == pdb_pair[i].res_name1) 
			{
				ipos=j;
				break;
			}
		}
		if(ipos==-1)
		{
			int s=syn1.size();
			syn1.resize(s+1);
			syn1[s].chain_id=pdb_pair[i].chain_id1;
			syn1[s].res_id=pdb_pair[i].res_id1;
			syn1[s].res_name=pdb_pair[i].res_name1;
			syn1[s].pdb_vdw_count=syn1[s].pdb_hb_count=syn1[s].pdb_salt_count=0;
			syn1[s].dcd_vdw_count=syn1[s].dcd_hb_count=syn1[s].dcd_salt_count=0;
			ipos=s;
		}
		syn1[ipos].pdb_vdw_count+=pdb_pair[i].vdw_count;
		syn1[ipos].pdb_hb_count+=pdb_pair[i].hb_count;
		syn1[ipos].pdb_salt_count+=pdb_pair[i].salt_count;
		syn1[ipos].dcd_vdw_count+=dcd_pair[i].vdw_count;
		syn1[ipos].dcd_hb_count+=dcd_pair[i].hb_count;
		syn1[ipos].dcd_salt_count+=dcd_pair[i].salt_count;

		ipos=-1;
		for(int j=0;j<syn2.size();j++)
		{
			if(syn2[j].chain_id == pdb_pair[i].chain_id2 && 
			syn2[j].res_id  == pdb_pair[i].res_id2 &&
			syn2[j].res_name == pdb_pair[i].res_name2) 
			{
				ipos=j;
				break;
			}
		}
		if(ipos==-1)
		{
			int s=syn2.size();
			syn2.resize(s+1);
			syn2[s].chain_id=pdb_pair[i].chain_id2;
			syn2[s].res_id=pdb_pair[i].res_id2;
			syn2[s].res_name=pdb_pair[i].res_name2;
			syn2[s].pdb_vdw_count=syn2[s].pdb_hb_count=syn2[s].pdb_salt_count=0;
			syn2[s].dcd_vdw_count=syn2[s].dcd_hb_count=syn2[s].dcd_salt_count=0;
			ipos=s;
		}
		syn2[ipos].pdb_vdw_count+=pdb_pair[i].vdw_count;
		syn2[ipos].pdb_hb_count+=pdb_pair[i].hb_count;
		syn2[ipos].pdb_salt_count+=pdb_pair[i].salt_count;
		syn2[ipos].dcd_vdw_count+=dcd_pair[i].vdw_count;
		syn2[ipos].dcd_hb_count+=dcd_pair[i].hb_count;
		syn2[ipos].dcd_salt_count+=dcd_pair[i].salt_count;
	}
	Log("#SYNOPSIS");
	Log(printString("#Res",-15)+" "+
		printString("VDW1",5)+" "+
		printString("HB1",5)+" "+printString("SALT1",5)+" "+
		printString("VDW2",5)+" "+
		printString("HB2",5)+" "+printString("SALT2",5));
	for(int i=0;i<syn1.size();i++)
	{
		char s1[100];
		sprintf(s1,"%s:%d:%s",syn1[i].chain_id.c_str(),
			syn1[i].res_id,syn1[i].res_name.c_str());
		Log(printString(s1,-15)+" "+printNumber(syn1[i].pdb_vdw_count,5)+" "+
		printNumber(syn1[i].pdb_hb_count,5)+" "+printNumber(syn1[i].pdb_salt_count,5)+" "+
		printNumber(syn1[i].dcd_vdw_count,5)+" "+
		printNumber(syn1[i].dcd_hb_count,5)+" "+printNumber(syn1[i].dcd_salt_count,5));
	}
	for(int i=0;i<syn2.size();i++)
	{
		char s1[100];
		sprintf(s1,"%s:%d:%s",syn2[i].chain_id.c_str(),
			syn2[i].res_id,syn2[i].res_name.c_str());
		Log(printString(s1,-15)+" "+printNumber(syn2[i].pdb_vdw_count,5)+" "+
		printNumber(syn2[i].pdb_hb_count,5)+" "+printNumber(syn2[i].pdb_salt_count,5)+" "+
		printNumber(syn2[i].dcd_vdw_count,5)+" "+
		printNumber(syn2[i].dcd_hb_count,5)+" "+printNumber(syn2[i].dcd_salt_count,5));
	}
}


CommandContact::~CommandContact()
{
}
