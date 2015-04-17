# include <global.h>
# include <math.h>
# include <command_pdbangle.h>
# include <pdb.h>

CommandPdbAngle::CommandPdbAngle()
	:Command("PDBTORSION")
{
	phi_flag=psi_flag=omega_flag=chi1_flag=0;
}

void	CommandPdbAngle::setAngleList(vector<string> &s)
{
	for(int i=0;i<s.size();i++) 
	{
		if(s[i]=="phi") phi_flag=1;
		if(s[i]=="psi") psi_flag=1;
		if(s[i]=="ome") omega_flag=1;
		if(s[i]=="chi1")  chi1_flag=1;
		if(s[i]=="all") phi_flag=psi_flag=omega_flag=chi1_flag=1;
		if(s[i]=="backbone") phi_flag = psi_flag=omega_flag=1;
		if(s[i]=="sidechain") chi1_flag=1;
	}
}

typedef struct
{
	int ipos;
	int phi;
	int psi;
	int omega;
	int chi1;
}TorInfo;

/*	Produce angle stat files.
 * */
void	CommandPdbAngle::Run()
{
	Team *angle_team = team;
	if(angle_team==NULL) Error(SeqError("pdbtors"));
	if(pdbpos.size()==0) 
	{
		Log("You must specify the pdb file with the -pdb option ");
		psf_error("You must specify the pdb file with the -pdb option ");
	}
	FILE *logfile=fopen("phij.log","w");
	if(!logfile) Error(WriteError("phij.log"));
	string chain;
	int firstAminoAcid = table[0].res_id;
	int lastAminoAcid = table[table.size()-1].res_id;
	vector<int> foundAminoAcid;
	vector<TorInfo> info;
	
	vector<int> frame;
	vector<Data> angle_phi;
	vector<Data> angle_psi;
	vector<Data> angle_omega;
	vector<Data> angle_chi1;
	string c1;
	vector<int> posAtom;
	angle_team->enumerateAtoms(table,posAtom);

	NoeDihedralStruct noed;
	vector<NoeDihedralStruct> NoeDihedral;

	DoubleVector PdbAngle;

	//pass1 
	for(int i=0;i<posAtom.size();i++)
	{
		int ipos=posAtom[i];
		chain = table[ipos].chain_id;
		if(isin(foundAminoAcid,table[ipos].res_id)==-1)
		{
			firstAminoAcid=-1;
			lastAminoAcid=-1;
			int  pos_N=-1,pos_C=-1,pos_CA=-1;
			for(int j=0;j<table.size();j++)
			{
				if(table[j].chain_id==chain)
				{
					if(firstAminoAcid==-1) 
					{
						firstAminoAcid=table[j].res_id;
						if(table[j].atom_name=="N") pos_N=j;
						if(table[j].atom_name=="CA") pos_CA=j;
					}
					if(table[j].res_id==firstAminoAcid)
					{
						if(table[j].atom_name=="N") pos_N=j;
						if(table[j].atom_name=="CA") pos_CA=j;
					}
					lastAminoAcid=table[j].res_id;
					if(table[j].atom_name=="C") pos_C=j;
				}
			}
			int cycleFlag = isConnected(table,pos_N+1,pos_C+1);
			int k=foundAminoAcid.size();
			foundAminoAcid.resize(k+1);
			foundAminoAcid[k]=table[ipos].res_id;
			info.resize(k+1);
			info[k].ipos = ipos;
			info[k].phi=info[k].psi=info[k].omega=info[k].chi1=0;
			if(phi_flag)
			{
			int x1=-1,x2=-2,x3=-1,x4=-1;
			//special case
			if(table[ipos].res_id==firstAminoAcid)
			{
				for(int l=0;l<posAtom.size();l++)
				{
					int lpos=posAtom[l];
					if(table[lpos].res_id==table[ipos].res_id)
					{ 
					if(table[lpos].atom_name=="CY") x1=table[lpos].atom_id;
					if(table[lpos].atom_name=="N")  x2=table[lpos].atom_id;
					if(table[lpos].atom_name=="CA") x3=table[lpos].atom_id;
					if(table[lpos].atom_name=="C")  x4=table[lpos].atom_id;
					}
					if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
				}
			}
			else
			for(int l=0;l<table.size();l++)
			{
				int lpos=l;
				if(  table[lpos].res_id==table[ipos].res_id-1 
					&& table[lpos].atom_name=="C")
				
				{
					x1=table[lpos].atom_id;
					break;
				}
			}

			if(x1==-1 && cycleFlag && table[ipos].res_id==firstAminoAcid) 
					x1=table[pos_C].atom_id;
				for(int l=0;l<posAtom.size();l++)
				{
					int lpos=posAtom[l];
					if(table[lpos].res_id==table[ipos].res_id)
					{
					  if(table[lpos].atom_name=="N") x2=table[lpos].atom_id;
					  if(table[lpos].atom_name=="CA")x3=table[lpos].atom_id;
					  if(table[lpos].atom_name=="C") x4=table[lpos].atom_id;
					}
					if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
				}
				if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1)
				{
					noed.atom1=x1;
					noed.atom2=x2;
					noed.atom3=x3;
					noed.atom4=x4;
					NoeDihedral.push_back(noed);
					PdbAngle.push_back(
						torsion(
						pdbpos[x1].x,pdbpos[x1].y,pdbpos[x1].z,
						pdbpos[x2].x,pdbpos[x2].y,pdbpos[x2].z,
						pdbpos[x3].x,pdbpos[x3].y,pdbpos[x3].z,
						pdbpos[x4].x,pdbpos[x4].y,pdbpos[x4].z
						)
					);
					Log("Phi "+table[x1].chain_id+" "+printNumber(table[x1].res_id));
					int s=info.size();
					info[s-1].phi=1;
				}
			}
			if(psi_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(table[ipos].res_id==lastAminoAcid)
				{
					for(int l=0;l<posAtom.size();l++)
					{
						int lpos=posAtom[l];
						if(table[lpos].res_id==table[ipos].res_id)
						{
							if(table[lpos].atom_name=="N") x1=table[lpos].atom_id;
							if(table[lpos].atom_name=="CA")  x2=table[lpos].atom_id;
							if(table[lpos].atom_name=="C") x3=table[lpos].atom_id;
							if(table[lpos].atom_name=="NT")  x4=table[lpos].atom_id;
						}
						if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
					}
				}
				else
				{
					for(int l=0;l<posAtom.size();l++)
					{
					int lpos=posAtom[l];
					if(table[lpos].res_id==table[ipos].res_id)
					{
					 if(table[lpos].atom_name=="N") x1=table[lpos].atom_id;
					 if(table[lpos].atom_name=="CA")x2=table[lpos].atom_id;
					 if(table[lpos].atom_name=="C") x3=table[lpos].atom_id;
					}
					}
					for(int l=0;l<table.size();l++)
					{
					int lpos=l;
					if(table[lpos].res_id==table[ipos].res_id+1 
							&& table[lpos].atom_name=="N")
					{
						x4=table[lpos].atom_id;
						break;
					}
					if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
					}
				}
				if(x4==-1 && cycleFlag && table[ipos].res_id==lastAminoAcid) 
					x4=table[pos_N].atom_id;
				if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1)
				{
					noed.atom1=x1;
					noed.atom2=x2;
					noed.atom3=x3;
					noed.atom4=x4;
					NoeDihedral.push_back(noed);
					PdbAngle.push_back(
						torsion(
						pdbpos[x1].x,pdbpos[x1].y,pdbpos[x1].z,
						pdbpos[x2].x,pdbpos[x2].y,pdbpos[x2].z,
						pdbpos[x3].x,pdbpos[x3].y,pdbpos[x3].z,
						pdbpos[x4].x,pdbpos[x4].y,pdbpos[x4].z
						)
					);
					Log("Psi "+table[x1].chain_id+" "+printNumber(table[x1].res_id));
					int s=info.size();
					info[s-1].psi=1;
				}
			}
			if(omega_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(table[ipos].res_id==firstAminoAcid)
				{
					for(int l=0;l<posAtom.size();l++)
					{
						int lpos=posAtom[l];
						if(table[lpos].res_id==table[ipos].res_id)
						{
							if(table[lpos].atom_name=="CAY") x1=table[lpos].atom_id;
							if(table[lpos].atom_name=="CY")  x2=table[lpos].atom_id;
							if(table[lpos].atom_name=="N") x3=table[lpos].atom_id;
							if(table[lpos].atom_name=="CA")  x4=table[lpos].atom_id;
						}
						if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
					}
				}
				else
				if(table[ipos].res_id==lastAminoAcid)
				{
					for(int l=0;l<posAtom.size();l++)
					{
						int lpos=posAtom[l];
						if(table[lpos].res_id==table[ipos].res_id)
						{
							if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
							if(table[lpos].atom_name=="C")  x2=table[lpos].atom_id;
							if(table[lpos].atom_name=="N") x3=table[lpos].atom_id;
							if(table[lpos].atom_name=="HT1")  x4=table[lpos].atom_id;
						}
						if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
					}
				}
				else
				{
					for(int l=0;l<posAtom.size();l++)
					{
					int lpos=posAtom[l];
					if(table[lpos].res_id==table[ipos].res_id 
					&& table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
					if(table[lpos].res_id==table[ipos].res_id 
					&& table[lpos].atom_name=="C") x2=table[lpos].atom_id;
					}
					for(int l=0;l<table.size();l++)
					{
					int lpos=l;
					if(table[lpos].res_id==table[ipos].res_id+1 
					&& table[lpos].atom_name=="N") x3=table[lpos].atom_id;
					if(table[lpos].res_id==table[ipos].res_id+1 
					&& table[lpos].atom_name=="CA") x4=table[lpos].atom_id;
					if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
					}
				}
				if(x3==-1 && x4==-1 && cycleFlag && table[ipos].res_id==lastAminoAcid)
				{
					x3=table[pos_N].atom_id;
					x4=table[pos_CA].atom_id;
				}
				if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1)
				{
					noed.atom1=x1;
					noed.atom2=x2;
					noed.atom3=x3;
					noed.atom4=x4;
					NoeDihedral.push_back(noed);

					PdbAngle.push_back(
						torsion(
						pdbpos[x1].x,pdbpos[x1].y,pdbpos[x1].z,
						pdbpos[x2].x,pdbpos[x2].y,pdbpos[x2].z,
						pdbpos[x3].x,pdbpos[x3].y,pdbpos[x3].z,
						pdbpos[x4].x,pdbpos[x4].y,pdbpos[x4].z
						)
					);
					Log("Omega "+table[x1].chain_id+" "+printNumber(table[x1].res_id));
					int s=info.size();
					info[s-1].omega=1;
				}
			}
			if(chi1_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				for(int l=0;l<posAtom.size();l++)
				{
					int lpos=posAtom[l];	
					if(table[lpos].res_id==table[ipos].res_id)
					{
					 if(table[lpos].res_name == "GLY" || table[lpos].res_name == "ALA") continue;
                                         
					 if(table[lpos].atom_name=="N") x1=table[lpos].atom_id;
					 if(table[lpos].atom_name=="CA")x2=table[lpos].atom_id;
					 if(table[lpos].atom_name=="CB") x3=table[lpos].atom_id;
					 if(table[lpos].atom_name=="CG") x4=table[lpos].atom_id;
					//exceptions
					 if(table[lpos].res_name == "CYS" && table[lpos].atom_name=="SG") x4=table[lpos].atom_id;
					 if(table[lpos].res_name == "SER" && table[lpos].atom_name=="OG") x4=table[lpos].atom_id;
					 if(table[lpos].res_name == "VAL" && table[lpos].atom_name=="CG1") x4=table[lpos].atom_id;
					 if(table[lpos].res_name == "ILE" && table[lpos].atom_name=="CG1") x4=table[lpos].atom_id;
					 if(table[lpos].res_name == "THR" && table[lpos].atom_name=="CG2") x4=table[lpos].atom_id;
					}
					if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
				}
				if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1)
				{
					noed.atom1=x1;
					noed.atom2=x2;
					noed.atom3=x3;
					noed.atom4=x4;
					NoeDihedral.push_back(noed);
					PdbAngle.push_back(
						torsion(
						pdbpos[x1].x,pdbpos[x1].y,pdbpos[x1].z,
						pdbpos[x2].x,pdbpos[x2].y,pdbpos[x2].z,
						pdbpos[x3].x,pdbpos[x3].y,pdbpos[x3].z,
						pdbpos[x4].x,pdbpos[x4].y,pdbpos[x4].z
						)
					);
					Log("Chi1 "+table[x1].chain_id+" "+printNumber(table[x1].res_id));
					int s=info.size();
					info[s-1].chi1=1;
				}
			}
		}
	}
	//pass2 
	dcd->getNoeDihedral(frame,NoeDihedral);
	int icount_dihedral=0;
	for(int i=0;i<foundAminoAcid.size();i++)
	{
		int ipos = info[i].ipos;
		if(info[i].phi==1) 
		{
			angle_phi=NoeDihedral[icount_dihedral].D;
			for(int j=0;j<angle_phi.size();j++)
			{
				Data a=fabs(angle_phi[j]-PdbAngle[icount_dihedral]);
				if(a>=180.0) a=360-a;
				angle_phi[j]=fabs(a);
			}
			icount_dihedral++;
		}
		if(info[i].psi==1) 
		{	
			angle_psi=NoeDihedral[icount_dihedral].D;
			for(int j=0;j<angle_psi.size();j++)
			{
				Data a=fabs(angle_psi[j]-PdbAngle[icount_dihedral]);
				if(a>=180.0) a=360-a;
				angle_psi[j]=fabs(a);
			}
			icount_dihedral++;
		}
		if(info[i].omega==1) 
		{
			angle_omega=NoeDihedral[icount_dihedral].D;
			for(int j=0;j<angle_omega.size();j++)
			{
				Data a=fabs(angle_omega[j]-PdbAngle[icount_dihedral]);
				if(a>=180.0) a=360-a;
				angle_omega[j]=fabs(a);
			}
			icount_dihedral++;
		}
		if(info[i].chi1==1) 
		{
			angle_chi1=NoeDihedral[icount_dihedral].D;
			for(int j=0;j<angle_chi1.size();j++)
			{
				Data a=fabs(angle_chi1[j]-PdbAngle[icount_dihedral]);
				if(a>=180.0) a=360-a;
				angle_chi1[j]=fabs(a);
			}
			icount_dihedral++;
		}
		FILE *fp;
		string filename,histname,statname,smoothname,histname2;
		filename = "pdb_tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".dat";
		statname = "pdb_tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".stat";
		histname = "pdb_tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".hist";
		smoothname = "pdb_tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".sda";
		histname2 = "pdb_tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".hist2";
		Data delta=bindihe;

		if(smooth_flag)
		{
			IntVector sf;
			DoubleVector sphi,spsi,somega,schi1;
			if(phi_flag && angle_phi.size()) makeAngleSmooth(frame,angle_phi,sf,sphi,smooth_start,smooth_step);
			if(psi_flag && angle_psi.size()) makeAngleSmooth(frame,angle_psi,sf,spsi,smooth_start,smooth_step);
			if(omega_flag && angle_omega.size()) makeAngleSmooth(frame,angle_omega,sf,somega,smooth_start,smooth_step);
			if(phi_flag && angle_chi1.size()) makeAngleSmooth(frame,angle_chi1,sf,schi1,smooth_start,smooth_step);
			fp=fopen(smoothname.c_str(),"w");
			if(!fp) Error(WriteError(smoothname));
			for(int l=0;l<sf.size();l++)
			{
				Print(fp,printFrame(sf[l])+" ");
				if(phi_flag && angle_phi.size()) 
					Print(fp,printAngle(sphi[l])+" ");
				if(psi_flag && angle_psi.size()) 
					Print(fp,printAngle(spsi[l])+" ");
				if(omega_flag && angle_omega.size()) 
					Print(fp,printAngle(somega[l])+" ");
				if(chi1_flag && angle_chi1.size()) 
					Print(fp,printAngle(schi1[l])+" ");
				PrintLine(fp,"");
			}
			fclose(fp);
		}
	
		fp=fopen(filename.c_str(),"w");
		FILE *stat=fopen(statname.c_str(),"w");
		if(!stat) Error(WriteError(statname));
		fclose(stat);

		int cols=0;
		Data avg_JHNHA=0.0;

		vector<HistStruct> st_phi;
		vector<HistStruct> st_psi;
		vector<HistStruct> st_omega;
		vector<HistStruct> st_chi1;

		Print(fp,"#");
		Print(fp,printString("frame ",8));
		Data amin,amax,aavg,astd;
		if(phi_flag && angle_phi.size()) 
		{
			getAngleVectorStatistics(angle_phi,amin,amax,aavg,astd);
			Print(fp,printString("phi ",8));
			if(histflag) makeHist(angle_phi,st_phi,0,180.0,bindihe);
			cols++;
			string header1="#Phi statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_phi);
		}
		if(psi_flag && angle_psi.size()) 
		{
			Print(fp,printString("psi ",8));
			getAngleVectorStatistics(angle_psi,amin,amax,aavg,astd);
			if(histflag) makeHist(angle_psi,st_psi,0,180.0,bindihe);
			cols++;
			string header1="#Psi statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_psi);
		}
		if(omega_flag && angle_omega.size()) 
		{
			Print(fp,printString("ome ",8));
			getAngleVectorStatistics(angle_omega,amin,amax,aavg,astd);
			if(histflag) makeHist(angle_omega,st_omega,0,180.0,bindihe);
			cols++;
			string header1="#Omega statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_omega);
		}
		if(chi1_flag && angle_chi1.size()) 
		{
			Print(fp,printString("ch1 ",8));
			getAngleVectorStatistics(angle_chi1,amin,amax,aavg,astd);
			if(histflag) makeHist(angle_chi1,st_chi1,0,180.0,bindihe);
			cols++;
			string header1="#Chi1 statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_chi1);
		}
		stat=fopen(statname.c_str(),"a");
		PrintLine(stat,"");
		PrintLine(fp,"");
		for(int l=0;l<frame.size();l++)
		{
			Print(fp,printFrame(frame[l])+" ");
			if(phi_flag && angle_phi.size())
			{
				Print(fp,printAngle(angle_phi[l])+" ");
			}
			if(psi_flag && angle_psi.size())
				Print(fp,printAngle(angle_psi[l])+" ");
			if(omega_flag && angle_omega.size())
				Print(fp,printAngle(angle_omega[l])+" ");
			if(chi1_flag && angle_chi1.size())
				Print(fp,printAngle(angle_chi1[l]));
			PrintLine(fp,"");
		}
		fclose(stat);
		PrintLine(logfile,printString(table[ipos].chain_id,3)+" "+printNumber(table[ipos].res_id,4)+" "+printString(table[ipos].res_name,3));
		fclose(fp);
	
		if(histflag)
		{
			fp=fopen(histname.c_str(),"w");
			if(!fp) Error(WriteError(histname));
			int histsize=0;
			int hist_use=0;
			Print(fp,"#");
			Print(fp,printString("center",6)+" ");
			Print(fp,printString(" ",3));
			if(phi_flag && angle_phi.size()) 
			{
				Print(fp,printString("phi(%)",8)+" "+printString("phi(#)",8)+" ");
				histsize=st_phi.size();
				hist_use=1;
			}
			if(psi_flag && angle_psi.size()) 
			{
				Print(fp,printString("psi(%)",8)+" "+printString("psi(#)",8)+" ");
				histsize=st_psi.size();
				hist_use=2;
			}
			if(omega_flag && angle_omega.size()) 
			{
				Print(fp,printString("ome(%)",8)+" "+printString("ome(#)",8)+" ");
				histsize=st_omega.size();
				hist_use=3;
			}
			if(chi1_flag && angle_chi1.size()) 
			{
				Print(fp,printString("ch1(%)",8)+" "+printString("chi1(#)",8)+" ");
				histsize=st_chi1.size();
				hist_use=4;
			}
			PrintLine(fp,"");
			for(int l=0;l<histsize;l++)
			{
				if(l==histsize-1) continue;
				if(hist_use==1) Print(fp,printAngle(st_phi[l].value)+" ");
				else
				if(hist_use==2) Print(fp,printAngle(st_psi[l].value)+" ");
				else
				if(hist_use==3) Print(fp,printAngle(st_omega[l].value)+" ");
				else		Print(fp,printAngle(st_chi1[l].value)+" ");
				Print(fp,printString(" ",5));
				if(phi_flag && angle_phi.size())
					Print(fp,printPercent(st_phi[l].percent)+" "+printFrame(st_phi[l].count)+"   ");
				if(psi_flag && angle_psi.size())
					Print(fp,printPercent(st_psi[l].percent)+" "+printFrame(st_psi[l].count)+"   ");
				if(omega_flag && angle_omega.size())
					Print(fp,printPercent(st_omega[l].percent)+" "+printFrame(st_omega[l].count)+"   ");
				if(chi1_flag && angle_chi1.size())
					Print(fp,printPercent(st_chi1[l].percent)+" "+printFrame(st_chi1[l].count)+"   ");
				PrintLine(fp,"");
			}
			fclose(fp);
		}
		if(angle_phi.size() && angle_psi.size() && histflag)
		{
			fp=fopen(histname2.c_str(),"w");
			if(!fp) Error(WriteError(histname2));
			for(int k=0;k<st_phi.size();k++)
			{
				for(int m=0;m<st_psi.size();m++)
				{
					int icount=0;
					for(int l=0;l<angle_psi.size();l++)
					{
						if(angle_phi[l]>=st_phi[k].value-bindihe/2.0 &&
						 angle_phi[l]<=st_phi[k].value+bindihe/2.0 &&
						 angle_psi[l]>=st_psi[m].value-bindihe/2.0 &&
						 angle_psi[l]<=st_psi[m].value+bindihe/2.0)
							icount++;
					}
					PrintLine(fp,printAngle(st_phi[k].value)+" "+
					 printAngle(st_psi[m].value)+" "+printFrame(icount));
				}
				PrintLine(fp,"");
			}
			fclose(fp);
		}
	}
	fclose(logfile);
}
