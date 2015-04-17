# include <command_energy.h>
# include <global.h>

CommandEnergy::CommandEnergy(vector<string> &s)
	:Command("ENERGY")
{
	list.resize(s.size());
	for(int i=0;i<s.size();i++) list[i]=s[i];
}


void CommandEnergy::Run()
{
	if(parfile=="") psf_error("You must specify the param file with the -par option");
	IntVector atomPos;
	if(team) team->enumerateAtoms(table,atomPos);
	else
	for(int i=0;i<table.size();i++) atomPos.push_back(i);
	IntVector F;
	DoubleVector2 E;
	dcd->getEnergy(F,E,atomPos,list);	
	FILE *fp=fopen("energy.dat","w");
	Print(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" ");	
	for(int i=0;i<list.size();i++)
		Print(fp,printString(list[i],12)+" ");
	PrintLine(fp,"");
	for(int i=0;i<F.size();i++)
	{
		Print(fp,printFrame(F[i])+" ");
		for(int j=0;j<list.size();j++)
			Print(fp,printNumber(E[j][i],12,4)+" ");
		PrintLine(fp,"");
	}
	fclose(fp);
}


CommandEnergy::~CommandEnergy()
{
}

