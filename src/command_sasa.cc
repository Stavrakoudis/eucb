# include <command_sasa.h>
# include <global.h>
CommandSasa::CommandSasa(Data s)
	:Command("SASA")
{
	srad = s;
}

void	CommandSasa::Run()
{
	IntVector posAtom;
	IntVector F;
	DoubleVector D;
	for(int i=0;i<table.size();i++)
	{
		if(table[i].atom_name[0]=='H') continue;
		string chain;
		if(!team || team->find(table[i],chain)) posAtom.push_back(i);
	}
	dcd->getSasa(F,posAtom,srad,D);
	string filename="sasa.dat";
	string statname="sasa.stat";
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp) psf_error("Can not open "+filename+" for writing");
	PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
		printDistanceHeader("Sasa"));
	for(int i=0;i<F.size();i++)
	{
		PrintLine(fp,printFrame(F[i])+" "+printDistance(D[i]));
	}
	fclose(fp);
	PrintStatFile(statname,"w",F.size()/statcount,"#Sasa statistics",F,D);
}

CommandSasa::~CommandSasa()
{
}
