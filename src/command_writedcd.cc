# include <command_writedcd.h>
# include <pdb.h>
# include <global.h>
		
extern void	writePdb(string filename,int frame,IntVector &atom_list);
static string basename(string s)
{
	int dotpos=-1;
	for(int i=0;i<s.size();i++) 
		if(s[i]=='.') {dotpos=i;break;}
	if(dotpos==-1 || dotpos==0) return "";
	string ret="";
	for(int i=0;i<dotpos;i++)
	{
		ret=ret+s[i];
	}
	return ret;
}

CommandWriteDcd::CommandWriteDcd(string s,int f)
	:Command("WRITEDCD")
{
	fname = s;
	fit =f ;
}

void	CommandWriteDcd::Run()
{
	if(fname  == dcdfile) 
		psf_error("The filename must be different from "+dcdfile);
	IntVector posAtom;
	if(team) team->enumerateAtoms(table,posAtom);
	else    for(int i=0;i<table.size();i++) posAtom.push_back(i);
	dcd->writeDcd(fname,posAtom,fit);
	if(!team) return;
	string psfname=basename(fname)+".psf";
	if(psfname==psffile || psfname==".psf")
		psf_error("The psf output file must be different than current psf file");
	write_psf(psfname,table,posAtom);
	string newpdbname=basename(fname)+".pdb";
	writePdb(newpdbname,1,posAtom);
}

CommandWriteDcd::~CommandWriteDcd()
{
}
