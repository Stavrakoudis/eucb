# include <global.h>
# include <psf.h>
# include <command_rmsf.h>

CommandRmsf::CommandRmsf()
	:Command("RMSF")
{
}

void 	CommandRmsf::printRMSF(string chain)
{
	vector<int> list;
	vector<Data>  rmsf;
	for(int i=0;i<table.size();i++)
	{
		string c;
		if(team && !team->find(table[i],c)) continue;
		if(team && team->getAtomNameSize()) 
		{
			if(table[i].chain_id==chain) list.push_back(i+1);
		}
		else
		if(table[i].atom_name=="CA" && table[i].chain_id==chain) list.push_back(i+1);
	}
	if(list.size()==0) return;
	rmsf.resize(list.size());
	dcd->printRMSF(list,rmsf);
	char outputname[1024];
	sprintf(outputname,"rmsf_%s.dat",chain.c_str());
	FILE *fp=fopen(outputname,"w");
	if(!fp) Error(WriteError(outputname));
	PrintLine(fp,"#"+printString("ResNumber",10)+"\t"+
		     printString("ResName",10)+"\t"+printString("RMSF",6));
	for(int i=0;i<list.size();i++)
	{
		PrintLine(fp,printNumber(table[list[i]-1].res_id,10)+"\t"+
			     printString(table[list[i]-1].res_name,10)+"\t"+
			     printNumber(rmsf[i],8,4));
	}
	fclose(fp);
}

void	CommandRmsf::Run()
{
	Team *rmsf_team = team;
	string chain;
	vector<string> all_chain;
	all_chain.resize(0);
	if(rmsf_team)           rmsf_team->enumerateChains(all_chain);
	if(all_chain.size()==0) enumerateChains(table,all_chain);
	for(int i=0;i<all_chain.size();i++) printRMSF(all_chain[i]);
}
