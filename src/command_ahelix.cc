# include <command_ahelix.h>
# include <global.h>

CommandHelix::CommandHelix(Data p)
	:Command("HELIX")
{
	percent = p;
}

static int	helix(Data p)
{
	if(p>=(-57-30) && p<=(-57+30)) return 1;
	return 0;
}

void	CommandHelix::Run()
{
	if(team == NULL) psf_error("You must specify the sequence with the -seq option ");
	vector<string> chain;
	team->enumerateChains(chain);	
	vector<NoeDihedralStruct> phi;
	vector<NoeDihedralStruct> psi;
	for(int i=0;i<chain.size();i++)
	{
		//int i=0; // here
		vector<int> res_id;
		team->enumerateResId(chain[i],res_id);
		
		
		for(int j=0;j<res_id.size();j++)
		{
			NoeDihedralStruct st;
			int id1=res_id[j];
			if(j!=0) 
			{
				st.atom1=getAtomPos(res_id[j-1],"C");
				st.atom2=getAtomPos(id1,"N");
				st.atom3=getAtomPos(id1,"CA");
				st.atom4=getAtomPos(id1,"C");
				phi.push_back(st);
			}
			if(j!=res_id.size()-1)
			{
				st.atom1=getAtomPos(id1,"N");
				st.atom2=getAtomPos(id1,"CA");
				st.atom3=getAtomPos(id1,"C");
				st.atom4=getAtomPos(res_id[j+1],"N");
				psi.push_back(st);
			}
		}
	}
	IntVector f;
	dcd->getNoeDihedral(f,phi);
	dcd->getNoeDihedral(f,psi);
	IntVector isHelix;
	isHelix.resize(phi.size());
	for(int j=0;j<phi.size();j++) isHelix[j]=0;
	for(int j=0;j<phi.size()-3;j++)
	{
		int icount=0;
		for(int i=0;i<f.size();i++)
		{
			double f1=phi[j].D[i];
			double f2=phi[j+1].D[i];
			double f3=phi[j+2].D[i];
			double p1=psi[j].D[i];
			double p2=psi[j+1].D[i];
			double p3=psi[j+2].D[i];
	 		if(helix(f1) && helix(f2) && helix(f3) && 	
				helix(p1) && helix(p2) && helix(p3)) icount++;
		}
		if(icount*1.0/f.size() >= percent)
		{
			isHelix[j]=1;
			isHelix[j+1]=1;
			isHelix[j+2]=1;
		}
	}
	for(int i=0;i<phi.size();i++)
	{
		IntVector Helix;
		int icount=i;
		while(icount <phi.size() && isHelix[icount])
		{
			Helix.push_back(icount);
			icount++;
		}
		i=icount;
		if(Helix.size()<3) continue;
		atom firstAtom = table[phi[Helix[0]].atom2-1];
		atom lastAtom  = table[phi[Helix[Helix.size()-1]].atom2-1];
		string filename="ahelix_"+firstAtom.chain_id+"_"+printNumber(firstAtom.res_id)+"_"+
			printNumber(lastAtom.res_id)+".dat";
		string statname="ahelix_"+firstAtom.chain_id+"_"+printNumber(firstAtom.res_id)+"_"+
			printNumber(lastAtom.res_id)+".stat";
		FILE *fp=fopen(filename.c_str(),"w");
		Print(fp,printString("#Frame",-FRAME_WIDTH)+" ");
		for(int j=0;j<Helix.size();j++)
		{
			int res=table[phi[Helix[j]].atom2-1].res_id; 
			Print(fp,printFrame(res)+" ");
		}
		PrintLine(fp,"");
		DoubleVector Acor;
		Acor.resize(f.size());	
		for(int j=0;j<f.size();j++)
		{
			Print(fp,printFrame(f[j])+" ");
			Acor[j]=0;
			for(int k=0;k<Helix.size();k++)
			{
				double d1=phi[Helix[k]].D[j];
				double d2=psi[Helix[k]].D[j];
				int value=0;
				if(helix(d1) && helix(d2)) value=1;
				Print(fp,printFrame(value)+" ");
				Acor[j]+=value;
			}	
			PrintLine(fp,"");
		}
		fclose(fp);
		PrintStatFile(statname,"w",f.size()/statcount,"#Count statistics",f,Acor);
	}
}
