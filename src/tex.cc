# include <tex.h>		
# include <global.h>
Tex::Tex(string name)
{
	filename = name;
	header.resize(0);
	Atom.resize(0);
	percent.resize(0);
	caption="";
	enable_atoms=1;
}

void	Tex::enableAtoms(int f)
{
	enable_atoms=f;
}

void	Tex::setCaption(string s)
{
	caption=s;
}

void	Tex::setHeader(vector<string> s)
{
	header.resize(s.size());
	header = s;
}

void	Tex::addAtomList(IntVector a,Data p)
{
	int s=Atom.size();
	Atom.resize(s+1);
	Atom[s]=a;
	percent.resize(s+1);
	percent[s]=p;
	isgroup.push_back(0);
}

void	Tex::addGroupList(vector<TexGroup> &g)
{
		//removed for now 
	return ;
	for(int i=0;i<g.size();i++)
	{
		int atom_id1=g[i].first;
		int atom_id2=g[i].second;
		for(int j=0;j<Atom.size();j++)
		{	
			int atom_id3=Atom[j][1];
			if(atom_id3==atom_id1)
			{
				percent[j]=g[j].percent;
				isgroup[j]=1;
			}
			if(atom_id3==atom_id2)
			{
				//mark as invalid
				Atom[j][0]=-1;
			}
		}
	}
	for(int i=0;i<Atom.size();i++)
	{
		if(Atom[i][0]==-1)
		{
			for(int j=i;j<Atom.size()-1;j++) 
			{
				Atom[j]=Atom[j+1];
				isgroup[j]=isgroup[j+1];
			}
			int s=Atom.size();
			Atom.resize(s-1);
			isgroup.resize(s-1);
		}
	}
}

void	Tex::sort()
{
	for(int i=0;i<Atom.size();i++)
	{
		for(int j=0;j<Atom.size()-1;j++)
		{
		
			if(table[Atom[j+1][0]-1].res_id<table[Atom[j][0]-1].res_id)
			{
				IntVector t;
				t=Atom[j];
				Atom[j]=Atom[j+1];
				Atom[j+1]=t;
				int p;
				p=isgroup[j];
				isgroup[j]=isgroup[j+1];
				isgroup[j+1]=p;
			}
		}
	}
	for(int i=0;i<Atom.size();i++)
	{
		for(int j=0;j<Atom.size()-1;j++)
		{
			if(table[Atom[j+1][0]-1].res_id==table[Atom[j][0]-1].res_id && 
				table[Atom[j+1][1]-1].res_id<table[Atom[j][1]-1].res_id)
			{
				IntVector t;
				t=Atom[j];
				Atom[j]=Atom[j+1];
				Atom[j+1]=t;
				int p;
				p=isgroup[j];
				isgroup[j]=isgroup[j+1];
				isgroup[j+1]=p;
			}
		}
	}
}

void	Tex::print()
{
	sort();
	FILE *fptex=fopen(filename.c_str(),"w");
	if(!fptex) psf_error("file "+filename+" not ready for opening");
	PrintLine(fptex,"\\documentclass{article}");
        PrintLine(fptex,"\\begin{document}");
	PrintLine(fptex,"\\noindent");
//	PrintLine(fptex,"\\textbf{Table X.}");
	PrintLine(fptex,"");
        PrintLine(fptex,"\\begin{table}");
        PrintLine(fptex,"\\caption{"+caption+"}");
	Print(fptex,"\\begin{tabular}{");
	for(int i=0;i<header.size();i++)
		Print(fptex," c ");
	if(pdbfile!="") 
	{
		if(header.size()==2)
			Print(fptex," r ");
		else	
			Print(fptex," r  r ");
	}
	PrintLine(fptex,"r }\\hline");
	for(int i=0;i<header.size();i++)
		Print(fptex,"\\textbf{"+header[i]+"} & ");
	if(pdbfile!="") 
	{
		if(header.size()==2)
			Print(fptex,"\\textbf{PDB} & ");
		else
			Print(fptex,"\\textbf{PDB1} & \\textbf{PDB2} & ");
	}
	PrintLine(fptex,"\\textbf{MD}\\\\ \\hline");
	for(int i=0;i<Atom.size();i++)
	{
		if(Atom[i][0]==-1) continue;
		for(int j=0;j<Atom[i].size();j++)
		{
			atom a=table[Atom[i][j]-1];
			if(isgroup[i] && j==1) 
			{
				Print(fptex,printAtomLatexGroup(a,enable_atoms)+"&");
			}
			else
				Print(fptex,printAtomLatex(a,enable_atoms)+"&");				
		}
		if(pdbfile!="")
		{
			PdbAtom pa,pb,pc;
			pa=pdb_table[Atom[i][0]-1];
			pb=pdb_table[Atom[i][1]-1];
			if(header.size()==2)
			{
				Data dist=getPdbDistance(pa,pb);
				Print(fptex,printNumber(dist,5,2)+" & ");
			}
			else
			{
				pc=pdb_table[Atom[i][2]-1];
				Data dist1=getPdbDistance(pa,pb);
				Data dist2=getPdbDistance(pb,pc);
				Print(fptex,printNumber(dist1,5,2)+" & "+printNumber(dist2,5,2)+" & ");
			}
		}
		PrintLine(fptex,printNumber(100.0*percent[i],5,1)+"\\\\");
	}
	PrintLine(fptex,"\\hline");
	PrintLine(fptex,"\\end{tabular}");
        PrintLine(fptex,"\\end{table}");
        PrintLine(fptex,"\\end{document}");
	fclose(fptex);
}
