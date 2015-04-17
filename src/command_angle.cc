# include <global.h>
# include <math.h>
# include <command_angle.h>

CommandAngle::CommandAngle()
	:Command("TORSION")
{
	phi_flag=psi_flag=omega_flag=chi1_flag=0;
	JHNHA_a=6.4;
	JHNHA_b=-1.4;
	JHNHA_c=1.9;
	zero_flag = 0;
	deg_flag=0;
}

void	CommandAngle::newDihedral(int x1,int x2,int x3,int x4,string angle)
{
	NoeDihedralStruct noed;
	noed.atom1=x1;
	noed.atom2=x2;
	noed.atom3=x3;
	noed.atom4=x4;
	NoeDihedral.push_back(noed);
	Log("Found "+angle+" torsion "+table[x1].chain_id+" "+
		printNumber(table[x1].res_id));
	Log("Torsion atoms are \n"+
		   printNumber(table[x1-1].atom_id)+" "+
		  printAtom3(table[x1-1])+
          	"\n"+printNumber(table[x2-1].atom_id)+" "+
		printAtom3(table[x2-1])+" "+
	  	"\n"+printNumber(table[x3-1].atom_id)+" "+
		printAtom3(table[x3-1])+" "+
	 	"\n"+printNumber(table[x4-1].atom_id)+" "+
		printAtom3(table[x4-1]));
}

int	CommandAngle::getAtomsPhi(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	if(table[ipos].res_id==firstAminoAcid)
	{
		for(int l=0;l<posAtom.size();l++)
		{
			int lpos=posAtom[l];
			if(table[lpos].res_id==table[ipos].res_id
		 		&& table[lpos].chain_id==table[ipos].chain_id)
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
		&& table[lpos].chain_id==table[ipos].chain_id
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
		if(table[lpos].res_id==table[ipos].res_id
		&& table[lpos].chain_id==table[ipos].chain_id)
		{
			  if(table[lpos].atom_name=="N") x2=table[lpos].atom_id;
			  if(table[lpos].atom_name=="CA")x3=table[lpos].atom_id;
			  if(table[lpos].atom_name=="C") x4=table[lpos].atom_id;
		}
		if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
	}
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}

int	CommandAngle::getAtomsPsi(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	if(table[ipos].res_id==lastAminoAcid)
	{
		for(int l=0;l<posAtom.size();l++)
		{
			int lpos=posAtom[l];
			if(table[lpos].res_id==table[ipos].res_id
			&& table[lpos].chain_id==table[ipos].chain_id)
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
			if(table[lpos].res_id==table[ipos].res_id
			&& table[lpos].chain_id==table[ipos].chain_id)
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
		&& table[lpos].chain_id==table[ipos].chain_id
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
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}


int	CommandAngle::getAtomsOmega(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	x1=-1;x2=-1;x3=-1;x4=-1;
	if(table[ipos].res_id==firstAminoAcid)
	{
		for(int l=0;l<posAtom.size();l++)
		{
			int lpos=posAtom[l];
			if(table[lpos].res_id==table[ipos].res_id
			&& table[lpos].chain_id==table[ipos].chain_id)
			{
				if(table[lpos].atom_name=="CAY") 	
					x1=table[lpos].atom_id;
				if(table[lpos].atom_name=="CY")  
					x2=table[lpos].atom_id;
				if(table[lpos].atom_name=="N") 
					x3=table[lpos].atom_id;
				if(table[lpos].atom_name=="CA")  
					x4=table[lpos].atom_id;
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
			if(table[lpos].res_id==table[ipos].res_id
			&& table[lpos].chain_id==table[ipos].chain_id)
			{
				if(table[lpos].atom_name=="CA") 
					x1=table[lpos].atom_id;
				if(table[lpos].atom_name=="C")  
					x2=table[lpos].atom_id;
				if(table[lpos].atom_name=="N") 
					x3=table[lpos].atom_id;
				if(table[lpos].atom_name=="HT1")  
					x4=table[lpos].atom_id;
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
				&& table[lpos].chain_id==table[ipos].chain_id
				&& table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			if(table[lpos].res_id==table[ipos].res_id 
				&& table[lpos].chain_id==table[ipos].chain_id
			&& table[lpos].atom_name=="C") x2=table[lpos].atom_id;
		}
		for(int l=0;l<table.size();l++)
		{
			int lpos=l;
			if(table[lpos].res_id==table[ipos].res_id+1 
			&& table[lpos].chain_id==table[ipos].chain_id
			&& table[lpos].atom_name=="N") x3=table[lpos].atom_id;
			if(table[lpos].res_id==table[ipos].res_id+1 
			&& table[lpos].chain_id==table[ipos].chain_id
			&& table[lpos].atom_name=="CA") x4=table[lpos].atom_id;
			if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
		}
	}
	if(x3==-1 && x4==-1 && cycleFlag && table[ipos].res_id==lastAminoAcid)
	{
		x3=table[pos_N].atom_id;
		x4=table[pos_CA].atom_id;
	}
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}

int	CommandAngle::getAtomsChi1(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	x1=-1;x2=-1;x3=-1;x4=-1;
	for(int l=0;l<posAtom.size();l++)
	{
		int lpos=posAtom[l];	
		if(table[lpos].res_id==table[ipos].res_id
		&& table[lpos].chain_id==table[ipos].chain_id)
		{
		 if(table[lpos].res_name == "GLY" || table[lpos].res_name == "ALA") continue;
                                       
		 if(table[lpos].atom_name=="N") x1=table[lpos].atom_id;
		 if(table[lpos].atom_name=="CA")x2=table[lpos].atom_id;
		 if(table[lpos].atom_name=="CB") x3=table[lpos].atom_id;
		 if(table[lpos].atom_name=="CG") 
			x4=table[lpos].atom_id;
		 if(table[lpos].res_name == "CYS" && table[lpos].atom_name=="SG") 
			x4=table[lpos].atom_id;
		 if(table[lpos].res_name == "SER" && table[lpos].atom_name=="OG") 
			x4=table[lpos].atom_id;
		 if(table[lpos].res_name == "VAL" && table[lpos].atom_name=="CG1") 
			x4=table[lpos].atom_id;
		 if(table[lpos].res_name == "ILE" && table[lpos].atom_name=="CG1") 
			x4=table[lpos].atom_id;
		 if(table[lpos].res_name == "THR" && table[lpos].atom_name=="CG2") 
			x4=table[lpos].atom_id;
		}
		if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
	}
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}

int	CommandAngle::getAtomsChi3(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	x1=-1,x2=-1,x3=-1,x4=-1;
	for(int l=0;l<posAtom.size();l++)
	{
		int lpos=posAtom[l];	
		if(table[lpos].res_id!=table[ipos].res_id
		|| table[lpos].chain_id!=table[ipos].chain_id) continue;
		 if(table[lpos].res_name == "ARG")
		{
			 if(table[lpos].atom_name=="CB") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="NE") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="LYS")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CE") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="GLU")
		{
			 if(table[lpos].atom_name=="CB") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="OE1") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="PRO")
		{
			 if(table[lpos].atom_name=="CB") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="N") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="PRO")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="MET")
		{
			 if(table[lpos].atom_name=="CB") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="SD") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CE") x4=table[lpos].atom_id;
		}
		if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
	}
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}

int	CommandAngle::getAtomsChi4(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	x1=-1;x2=-1;x3=-1;x4=-1;
	for(int l=0;l<posAtom.size();l++)
	{
		int lpos=posAtom[l];	
		if(table[lpos].res_id!=table[ipos].res_id
		|| table[lpos].chain_id!=table[ipos].chain_id) continue;
		 if(table[lpos].res_name == "ARG")
		{
			 if(table[lpos].atom_name=="CG") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="NE") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CZ") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="LYS")
		{
			 if(table[lpos].atom_name=="CG") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CE") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="NZ") x4=table[lpos].atom_id;
		}
		if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
	}
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}

int	CommandAngle::getAtomsChi2(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4)
{
	x1=-1;x2=-1;x3=-1;x4=-1;
	for(int l=0;l<posAtom.size();l++)
	{
		int lpos=posAtom[l];	
		if(table[lpos].res_id!=table[ipos].res_id
			|| table[lpos].chain_id!=table[ipos].chain_id) continue;
		atom a=table[lpos];
		 if(a.res_name == "ARG" || a.res_name=="CYS")
		 {
		 	if(a.atom_name=="CA") x1=a.atom_id;
		 	if(a.atom_name=="CB") x2=a.atom_id;
			if(a.atom_name=="CG") x3=a.atom_id;
		 	if(a.atom_name=="CD") x4=a.atom_id;
		}
		if(table[lpos].res_name=="ASP")
		{
		 	if(a.atom_name=="CA") x1=a.atom_id;
		 	if(a.atom_name=="CB") x2=a.atom_id;
		 	if(a.atom_name=="CG") x3=a.atom_id;
		 	if(a.atom_name=="OD1")x4=a.atom_id;
		}
		if(table[lpos].res_name=="GLU")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="ASN")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="ND2") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="PRO")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="ILE")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG1") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="LEU")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD2") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="MET")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="SD") x4=table[lpos].atom_id;
		}
		if(table[lpos].res_name=="PHE" || table[lpos].res_name=="TYR")
		{
			 if(table[lpos].atom_name=="CA") x1=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CB")x2=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CG") x3=table[lpos].atom_id;
			 if(table[lpos].atom_name=="CD1") x4=table[lpos].atom_id;
		}
		if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) break;
	}
	if(x1!=-1 && x2!=-1 && x3!=-1 && x4!=-1) return 1;
	return 0;
}

void	CommandAngle::setJHNHA(Data a,Data b,Data c)
{
	JHNHA_a=a;
	JHNHA_b=b;
	JHNHA_c=c;
}

void	CommandAngle::setAngleList(vector<string> &s)
{
	for(int i=0;i<s.size();i++) 
	{
		if(s[i]=="phi") phi_flag=1;
		else
		if(s[i]=="psi") psi_flag=1;
		else
		if(s[i]=="ome") omega_flag=1;
		else
		if(s[i]=="chi1")  chi1_flag=1;
		else
		if(s[i]=="chi2") chi2_flag=1;
		else
		if(s[i]=="chi3") chi3_flag=1;
		else
		if(s[i]=="chi4") chi4_flag=1;
		else
		if(s[i]=="all") phi_flag=psi_flag=omega_flag=chi1_flag=chi2_flag=chi3_flag=chi4_flag=1;
		else
		if(s[i]=="backbone") phi_flag = psi_flag=omega_flag=1;
		else
		if(s[i]=="sidechain") chi1_flag=chi2_flag=chi3_flag=chi4_flag=1;
		else
		if(s[i]=="0") zero_flag =1;
		else
		if(s[i]=="pi") deg_flag=1;
		else
		{
			Log("Unrecognized torsion option "+s[i]);
			psf_error("Unrecognized torsion option "+s[i]);
		}
	}
}

typedef struct
{
	int ipos;
	int phi;
	int psi;
	int omega;
	int chi1;
	int chi2;
	int chi3;
	int chi4;
}TorInfo;

/*	Produce angle stat files.
 * */
void	CommandAngle::Run()
{
	Team *angle_team = team;
	if(angle_team==NULL)  	Error(SeqError("tors"));
	FILE *logfile=fopen("JHNHA.log","w");
	if(!logfile)  		Error(WriteError("JHNHA.log"));
	string chain;
	firstAminoAcid = table[0].res_id;
	lastAminoAcid = table[table.size()-1].res_id;
	vector<int> foundAminoAcid;
	vector<TorInfo> info;
	
	IntVector frame;
	DoubleVector angle_phi;
	DoubleVector angle_psi;
	DoubleVector angle_omega;
	DoubleVector angle_chi1;
	DoubleVector angle_chi2;
	DoubleVector angle_chi3;
	DoubleVector angle_chi4;
	string c1;
	IntVector posAtom;
	angle_team->enumerateAtoms(table,posAtom);

	NoeDihedralStruct noed;

	//pass1 
	for(int i=0;i<posAtom.size();i++)
	{
		int ipos=posAtom[i];
		chain = table[ipos].chain_id;
		if(isin(foundAminoAcid,table[ipos].res_id)==-1)
		{
			firstAminoAcid=-1;
			lastAminoAcid=-1;
			pos_N=-1,pos_C=-1,pos_CA=-1;
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
			cycleFlag = isConnected(table,pos_N+1,pos_C+1);
			int k=foundAminoAcid.size();
			foundAminoAcid.resize(k+1);
			foundAminoAcid[k]=table[ipos].res_id;
			info.resize(k+1);
			info[k].ipos = ipos;
			info[k].phi=info[k].psi=info[k].omega=info[k].chi1=info[k].chi2=info[k].chi3=info[k].chi4=0;
			if(phi_flag)
			{
				int x1=-1,x2=-2,x3=-1,x4=-1;
				if(getAtomsPhi(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"phi");
					int s=info.size();
					info[s-1].phi=1;
				}
			}
			if(psi_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(getAtomsPsi(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"psi");
					int s=info.size();
					info[s-1].psi=1;
				}
			}
			if(omega_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(getAtomsOmega(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"omega");
					int s=info.size();
					info[s-1].omega=1;
				}
			}
			if(chi2_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(getAtomsChi2(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"chi2");
					int s=info.size();
					info[s-1].chi2=1;
				}
			}
			if(chi3_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(getAtomsChi3(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"chi3");
					int s=info.size();
					info[s-1].chi3=1;
				}
			}
			if(chi4_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(getAtomsChi4(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"chi4");
					int s=info.size();
					info[s-1].chi4=1;
				}
			}
			if(chi1_flag)
			{
				int x1=-1,x2=-1,x3=-1,x4=-1;
				if(getAtomsChi1(ipos,posAtom,x1,x2,x3,x4))
				{
					newDihedral(x1,x2,x3,x4,"chi1");
					int s=info.size();
					info[s-1].chi1=1;
				}
			}
		}
	}
	//pass2 
	dcd->getNoeDihedral(frame,NoeDihedral);
	//transform angles
	if(zero_flag) zero_angle_flag = 1;
	if(deg_flag)  degrees_angle_flag=1;

	int icount_dihedral=0;
	for(int i=0;i<foundAminoAcid.size();i++)
	{
		int ipos = info[i].ipos;
		if(info[i].phi==1) {angle_phi=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].psi==1) {angle_psi=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].omega==1) {angle_omega=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].chi1==1) {angle_chi1=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].chi2==1) {angle_chi2=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].chi3==1) {angle_chi3=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].chi4==1) {angle_chi4=NoeDihedral[icount_dihedral].D;icount_dihedral++;}
		if(info[i].phi==0 && info[i].psi==0 && info[i].omega==0 && 
		info[i].chi1==0 && info[i].chi2==0 && info[i].chi3==0 && info[i].chi4==0) continue;
		FILE *fp;
		string filename,histname,statname,smoothname,histname2;
		filename="tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".dat";
		statname="tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".stat";
		histname="tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".hist";
		smoothname="tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".sda";
		histname2="tors_"+table[ipos].chain_id+"_"+printNumber(table[ipos].res_id)+".hist2";
		Data minvalue=-180.0;
		Data maxvalue= 180.0;
		Data delta=bindihe;

		if(smooth_flag)
		{
			IntVector sf;
			DoubleVector sphi,spsi,somega,schi1,schi2,schi3,schi4;
			if(info[i].phi) 
			 makeAngleSmooth(frame,angle_phi,sf,sphi,smooth_start,smooth_step);
			if(info[i].psi) 
			 makeAngleSmooth(frame,angle_psi,sf,spsi,smooth_start,smooth_step);
			if(info[i].omega) 
			 makeAngleSmooth(frame,angle_omega,sf,somega,smooth_start,smooth_step);
			if(info[i].chi1) 
			 makeAngleSmooth(frame,angle_chi1,sf,schi1,smooth_start,smooth_step);
			if(info[i].chi2) 
			 makeAngleSmooth(frame,angle_chi2,sf,schi2,smooth_start,smooth_step);
			if(info[i].chi3) 
			 makeAngleSmooth(frame,angle_chi3,sf,schi3,smooth_start,smooth_step);
			if(info[i].chi4) 
			 makeAngleSmooth(frame,angle_chi4,sf,schi4,smooth_start,smooth_step);
			fp=fopen(smoothname.c_str(),"w");
			if(!fp) Error(WriteError(smoothname));
			for(int l=0;l<sf.size();l++)
			{
				Print(fp,printFrame(sf[l])+" ");
				if(info[i].phi)   Print(fp,printAngle(sphi[l])+" ");
				if(info[i].psi)   Print(fp,printAngle(spsi[l])+" ");
				if(info[i].omega) Print(fp,printAngle(somega[l])+" ");
				if(info[i].chi1)  Print(fp,printAngle(schi1[l])+" ");
				if(info[i].chi2)  Print(fp,printAngle(schi2[l])+" ");
				if(info[i].chi3)  Print(fp,printAngle(schi3[l])+" ");
				if(info[i].chi4)  Print(fp,printAngle(schi4[l])+" ");
				PrintLine(fp,"");
			}
			fclose(fp);
		}
	
		fp=fopen(filename.c_str(),"w");
		if(!fp) Error(WriteError(filename));
		FILE *stat=fopen(statname.c_str(),"w");
		if(!stat) Error(WriteError(statname));
		fclose(stat);

		Data avg_JHNHA=0.0;

		vector<HistStruct> st_phi;
		vector<HistStruct> st_psi;
		vector<HistStruct> st_omega;
		vector<HistStruct> st_chi1;
		vector<HistStruct> st_chi2;
		vector<HistStruct> st_chi3;
		vector<HistStruct> st_chi4;

		Print(fp,"#");
		Print(fp,printString("Frame",FRAME_WIDTH-1));
		Data amin,amax,aavg,astd;
		Data angle_min=1e+10;
		Data angle_max=-1e+10;
		if(info[i].phi)
		{
			getAngleVectorStatistics(angle_phi,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		if(info[i].psi)
		{
			getAngleVectorStatistics(angle_psi,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		if(info[i].omega)
		{
			getAngleVectorStatistics(angle_omega,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		if(info[i].chi1)
		{
			getAngleVectorStatistics(angle_chi1,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		if(info[i].chi2)
		{
			getAngleVectorStatistics(angle_chi2,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		if(info[i].chi3)
		{
			getAngleVectorStatistics(angle_chi3,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		if(info[i].chi4)
		{
			getAngleVectorStatistics(angle_chi4,amin,amax,aavg,astd);
			if(amin<angle_min) angle_min=amin;
			if(amax>angle_max) angle_max=amax;
		}
		
		
		if(info[i].phi)
		{
			Print(fp,printTorsionHeader("phi")+" ");
			if(histflag) makeHist(angle_phi,st_phi,angle_min,angle_max,bindihe);
			string header1="#Phi statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_phi);
		}
		if(info[i].psi)
		{
			Print(fp,printTorsionHeader("psi")+" ");
			if(histflag) makeHist(angle_psi,st_psi,angle_min,angle_max,bindihe);
			string header1="#Psi statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_psi);
		}
		if(info[i].omega)
		{
			Print(fp,printTorsionHeader("ome")+" ");
			if(histflag) makeHist(angle_omega,st_omega,angle_min,angle_max,bindihe);
			string header1="#Omega statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_omega);
		}
		if(info[i].chi1)
		{
			Print(fp,printTorsionHeader("chi1")+" ");
			if(histflag) makeHist(angle_chi1,st_chi1,angle_min,angle_max,bindihe);
			string header1="#Chi1 statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_chi1);
		}
		if(info[i].chi2)
		{
			Print(fp,printTorsionHeader("chi2")+" ");
			if(histflag) makeHist(angle_chi2,st_chi2,angle_min,angle_max,bindihe);
			string header1="#Chi2 statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_chi2);
		}
		if(info[i].chi3)
		{
			Print(fp,printTorsionHeader("chi3")+" ");
			if(histflag) makeHist(angle_chi3,st_chi3,angle_min,angle_max,bindihe);
			string header1="#Chi3 statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_chi3);
		}
		if(info[i].chi4)
		{
			Print(fp,printTorsionHeader("chi4")+" ");
			if(histflag) makeHist(angle_chi4,st_chi4,angle_min,angle_max,bindihe);
			string header1="#Chi4 statistics";
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header1,frame,angle_chi4);
		}
		stat=fopen(statname.c_str(),"a");
		PrintLine(stat,"");
		PrintLine(fp,"");
		IntVector res;
		IntVector start;
		IntVector last;
		DoubleVector jn;
		SplitTable(frame.size()/statcount,start,last,frame,frame,res);
		jn.resize(start.size());
		for(int l=0;l<jn.size();l++) jn[l]=0.0;
		DoubleVector jn_table;
		jn_table.resize(frame.size());
		

		for(int l=0;l<frame.size();l++)
		{
			int lpos=-1;
			int low,upper;
			for(int k1=0;k1<start.size();k1++) 
			{
				low=start[k1];
				upper=last[k1];
				if(frame[l]>=low && frame[l]<=upper)	
				{
					lpos = k1;
					break;
				}
			}	
			
			Print(fp,printFrame(frame[l])+" ");
			if(info[i].phi)
			{
				Data jhnha=JHNHA_a * 
					pow(cos(deg2rad(angle_phi[l])-M_PI/3.0),2.0)+
					JHNHA_b*cos(deg2rad(angle_phi[l])-M_PI/3.0)+JHNHA_c;
				avg_JHNHA+=jhnha;
				jn[lpos]+=jhnha/(upper-low+1);
				jn_table[l]=jhnha;
				Print(fp,printAngle(angle_phi[l])+" ");
			}
			if(info[i].psi)   Print(fp,printAngle(angle_psi[l])+" ");
			if(info[i].omega) Print(fp,printAngle(angle_omega[l])+" ");
			if(info[i].chi1)  Print(fp,printAngle(angle_chi1[l]));
			if(info[i].chi2)  Print(fp,printAngle(angle_chi2[l]));
			if(info[i].chi3)  Print(fp,printAngle(angle_chi3[l]));
			if(info[i].chi4)  Print(fp,printAngle(angle_chi4[l]));
			PrintLine(fp,"");
		}
		avg_JHNHA/=frame.size();
		fclose(stat);
		if(info[i].phi) 
		{
			PrintStatFile(statname,"a",frame.size()/statcount,"#JHNHA Statistics",
				frame,jn_table);
		}
		PrintLine(logfile,printAtom4(table[ipos])+" "+printNumber(avg_JHNHA,8,2));
		fclose(fp);
	
		if(histflag)
		{
			fp=fopen(histname.c_str(),"w");
			if(!fp) Error(WriteError(histname));
			int histsize=0;
			int hist_use=0;
			Print(fp,printTorsionHeader("#center")+" ");
			if(info[i].phi)
			{
				Print(fp,printPercentHeader("phi(%)")+" "+
					 printFrameHeader("phi(#)")+" ");
				histsize=st_phi.size();
				hist_use=1;
			}
			if(info[i].psi)
			{
				Print(fp,printPercentHeader("psi(%)")+" "+
					 printFrameHeader("psi(#)")+" ");
				histsize=st_psi.size();
				hist_use=2;
			}
			if(info[i].omega)
			{
				Print(fp,printPercentHeader("ome(%)")+" "+
					 printFrameHeader("ome(#)")+" ");
				histsize=st_omega.size();
				hist_use=3;
			}
			if(info[i].chi1)
			{
				Print(fp,printPercentHeader("ch1(%)")+" "+
					 printFrameHeader("chi1(#)")+" ");
				histsize=st_chi1.size();
				hist_use=4;
			}
			if(info[i].chi2)
			{
				Print(fp,printPercentHeader("ch2(%)")+" "+
					 printFrameHeader("chi2(#)")+" ");
				histsize=st_chi2.size();
				hist_use=5;
			}
			if(info[i].chi3)
			{
				Print(fp,printPercentHeader("ch3(%)")+" "+
					 printFrameHeader("chi3(#)")+" ");
				histsize=st_chi3.size();
				hist_use=6;
			}
			if(info[i].chi4)
			{
				Print(fp,printPercentHeader("ch4(%)")+" "+
					 printFrameHeader("chi4(#)")+" ");
				histsize=st_chi4.size();
				hist_use=7;
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
				else
				if(hist_use==4)	Print(fp,printAngle(st_chi1[l].value)+" ");
				else
				if(hist_use==5)	Print(fp,printAngle(st_chi2[l].value)+" ");
				else
				if(hist_use==6) Print(fp,printAngle(st_chi3[l].value)+" ");
				else
				if(hist_use==7) Print(fp,printAngle(st_chi4[l].value)+" ");
				if(info[i].phi)
					Print(fp,printPercent(st_phi[l].percent)+" "+
					printFrame(st_phi[l].count)+" ");
				if(info[i].psi)
					Print(fp,printPercent(st_psi[l].percent)+" "+
					printFrame(st_psi[l].count)+" ");
				if(info[i].omega)
					Print(fp,printPercent(st_omega[l].percent)+" "+
					printFrame(st_omega[l].count)+" ");
				if(info[i].chi1)
					Print(fp,printPercent(st_chi1[l].percent)+" "+
					printFrame(st_chi1[l].count)+" ");
				if(info[i].chi2)
					Print(fp,printPercent(st_chi2[l].percent)+" "+
					printFrame(st_chi2[l].count)+" ");
				if(info[i].chi3)
					Print(fp,printPercent(st_chi3[l].percent)+" "+
					printFrame(st_chi3[l].count)+" ");
				if(info[i].chi4)
					Print(fp,printPercent(st_chi4[l].percent)+" "+
					printFrame(st_chi4[l].count)+" ");
				PrintLine(fp,"");
			}
			fclose(fp);
		}
		if(info[i].phi && info[i].psi && histflag)
		{
			fp=fopen(histname2.c_str(),"w");
			if(!fp)  Error(WriteError(histname2));
			Print(fp,"#");
			PrintLine(fp,printString("Phi",ANGLE_WIDTH-1)+" "+
				     printTorsionHeader("Psi")+" "+
				     printFrameHeader("Common"));
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
	if(zero_flag) zero_angle_flag = 0;
	if(deg_flag)  degrees_angle_flag=0;
}
