# include <command_watpdb.h>
# include <pdb.h>
# include <global.h>
# include <math.h>
# include <donors.h>

CommandWatpdb::CommandWatpdb(Data d)
	:Command("WATPDB")
{
	distance = d;
}


static Data	getDistance(DoubleVector a,DoubleVector b)
{
	Data s=0.0;
	for(int i=0;i<a.size();i++) s=s+(a[i]-b[i])*(a[i]-b[i]);
	return sqrt(s);
}

void CommandWatpdb::Run()
{
	if(pdbpos.size() == 0) 
	{
		Log("You must load the pdb file with the -pdb command");
		psf_error("You must load the pdb file with the -pdb command");
	}
	if(team == NULL)
	{
		psf_error("You should provide a sequence with the -seq option ");
	}
	string basename="";
	int index;
	for(index=pdb_name.size()-1;pdb_name[index]!='.';index--);
	index--;
	while(index>=0 && pdb_name[index]!='/') 
	{
		basename=pdb_name[index]+basename;
		index--;
	}
	DoubleVector2  D;
	IntVector atom_list;
	IntVector water;
	if(team) team->enumerateAtoms(table,atom_list);
	else	
	{
		atom_list.resize(table.size());
		for(int i=0;i<atom_list.size();i++) atom_list[i]=i;
	}
	IntVector flagwater;
	for(int i=0;i<table.size();i++)
	{
		if((isWater(table[i]) && table[i].atom_name=="OH2")  ||
			table[i].atom_name == "CLA" || table[i].atom_name=="SOD")
		{
			water.push_back(i);
			flagwater.push_back(0);
		}
	}
	Log("Start from "+printNumber(first)+" to "+printNumber(last)+" frame ");
	for(int i=first;i<=last;i+=step)
	{
		for(int k=0;k<water.size();k++) flagwater[k]=0;
		dcd->fetchAtoms(i,D);
		string filename = basename+printNumber0(i,5)+".pdb";
		FILE *fp=fopen(filename.c_str(),"w");
		if(!fp) Error(WriteError(filename));
		for(int j=0;j<atom_list.size();j++)
		{
			PrintAtomInPdb(fp,D,atom_list[j]);
		}
		for(int j=0;j<atom_list.size();j++)
		{
			int dpos=isin(Donor,atom_list[j]);
			int apos=isin(Acceptor,atom_list[j]);
			if(dpos==-1 && apos==-1) continue;
			for(int k=0;k<water.size();k++)
			{
				double d=getDistance(D[atom_list[j]],D[water[k]]);
				if(d<=distance && !flagwater[k])
				{
					PrintAtomInPdb(fp,D,water[k]);
					flagwater[k]=1;
					atom w=table[water[k]];
					for(int m=0;m<w.bond.size();m++)
						PrintAtomInPdb(fp,D,w.bond[m]-1);
				}
			}
		}	
		PrintLine(fp,"END");
		fclose(fp);
	}
}

void	CommandWatpdb::PrintAtomInPdb(FILE *fp,DoubleVector2 &D,int id)
{
	Print(fp,printString("ATOM",-5)+" ");
	Print(fp,printNumber(id+1,5)+"  ");
	PdbAtom a=pdb_table[id];
	string s=a.atom_name;
	for(int k=s.size();k<=4;k++) s=s+" ";
	Print(fp,printString(s,-4));
	Print(fp,printString(a.res_name,-3)+" ");
	Print(fp,printString(a.pdb_chain,-1)+" ");
	Print(fp,printNumber(a.res_id,4)+"   ");
	Print(fp,printNumber(D[id][0],8,3)+""+
	printNumber(D[id][1],8,3)+""+
	printNumber(D[id][2],8,3)+" "+
	printNumber(a.val1,5,2)+" "+
	printNumber(a.val2,5,2)+" "+

	printString(a.chain_id,6)+" "+
	printString(a.pdb_letter,4));
	PrintLine(fp,"");
}

