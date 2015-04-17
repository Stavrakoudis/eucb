# include <command_pdbwrite.h>
# include <pdb.h>
# include <global.h>

CommandPdbWrite::CommandPdbWrite()
	:Command("PDBWRITE")
{
}


void	writePdb(string filename,int frame,IntVector &atom_list)
{
	int i=frame;
	vector<DoubleVector>  D;
		dcd->fetchAtoms(i,D);
		FILE *fp=fopen(filename.c_str(),"w");
		if(!fp) psf_error("Unable to write "+filename);
		for(int j=0;j<atom_list.size();j++)
		{
			Print(fp,printString("ATOM",-5));
			Print(fp,printNumber(j+1,5)+"  ");
			string s=pdb_table[atom_list[j]].atom_name;
			for(int k=s.size();k<=4;k++) s=s+" ";
			Print(fp,printString(s,-4));
			Print(fp,printString(pdb_table[atom_list[j]].res_name,-3)+" ");
			Print(fp,printString(pdb_table[atom_list[j]].pdb_chain,-1)+" ");
			Print(fp,printNumber(pdb_table[atom_list[j]].res_id,4)+"   ");
			Print(fp,printNumber(D[atom_list[j]][0],8,3)+""+
				printNumber(D[atom_list[j]][1],8,3)+""+
				printNumber(D[atom_list[j]][2],8,3)+" "+
				printNumber(pdb_table[atom_list[j]].val1,5,2)+" "+
				printNumber(pdb_table[atom_list[j]].val2,5,2)+" "+

				printString(pdb_table[atom_list[j]].chain_id,6)+" "+
				printString(pdb_table[atom_list[j]].pdb_letter,4));
			PrintLine(fp,"");
		}	
		PrintLine(fp,"END");
		fclose(fp);
}

void CommandPdbWrite::Run()
{
	if(pdbpos.size() == 0) 
	{
		Log("You must load the pdb file with the -pdb command");
		psf_error("You must load the pdb file with the -pdb command");
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
	IntVector atom_list;
	if(team) team->enumerateAtoms(table,atom_list);
	else	
	{
		atom_list.resize(table.size());
		for(int i=0;i<atom_list.size();i++) atom_list[i]=i;
	}
	Log("Start from "+printNumber(first)+" to "+printNumber(last)+" frame ");
	for(int i=first;i<=last;i+=step)
	{
		string filename = basename+printNumber0(i,5)+".pdb";
		writePdb(filename,i,atom_list);
	
	}
}
