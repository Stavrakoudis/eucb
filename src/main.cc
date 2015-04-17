# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include <global.h>
# include <getoptions.h>
# include <donors.h>
# include <command.h>
# include <helpfiles.h>
# include <pdb.h>

char **Argv;
int  Argc;

/*	Perform the initialization of the system:
 *		1. Open the psf file
 *		2. Open the dcd file
 *		3. If the donor - acceptor or backbone facility is activated
 *		   open and parse the donorsrc file.
 * */
extern int commandPos(string );
void	init()
{
	parse_psf((char*)psffile.c_str(),table);
	if(pdbfile!="")
	{
		int code = readPdb(pdbfile);
		if(!code)
		{
			FILE *fp=fopen(pdbfile.c_str(),"r");
			if(fp) 
			{
				fclose(fp);
				psf_warning("Incompatible psf and pdb files .");
			}
			else pdbfile="";	
		}
	}
	if(donor_act || backbone_flag || watbridge_flag || commandPos("CONTACT")!=-1 || commandPos("WATPDB")!=-1 || commandPos("IHBONDS")!=-1)
	{
		parse_donorsrc(table);
	}
	dcd = new Dcd((char*)dcdfile.c_str());
	if(dcd->getnatoms()!=table.size()) 
	{
		delete dcd;
		return psf_error("Incompatible number of atoms!");
	}
	dcd->setFrameStart(first);
	dcd->setFrameEnd(last);
	dcd->setFrameStep(step);
}

/*	Perform the system searches. Check what flags are active 
 *	and perform the correspondig searches.
 * */
void	run()
{
	for(int i=0;i<command_list.size();i++) command_list[i]->Run();
}	

/*	Terminate the application:
 *		1. Close the psf and dcd files.
 *		2. Delete the previously allocated teams.
 * */
void 	done()
{
	for(int i=0;i<command_list.size();i++) delete command_list[i];
	delete dcd;
	if(acceptorTeam) 	delete acceptorTeam;
	if(donorTeam)    	delete donorTeam;
}

/*	The main function of the program:
 *		1. Parse the command line.
 *		2. Initialize the system.
 *		3. Perform the necessary searches.
 *		4. Terminate the system.
 * */
int main(int argc,char **argv)
{
	makeHelpTables();
	Argv = argv;
	Argc = argc;
	parse_cmd_line(argc,argv);
	init();
	run();
	done();
	return 0;
}
