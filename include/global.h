# ifndef __GLOBAL__H
# define __GLOBAL__H
# include <vector>
# include <psf.h>
# include <dcd.h>
# include <team.h>
# include <string>

/*
 * 	Donor acceptor definitions.
 * 	donor_act: Set to 1 if the -donor is activated.
 * 	acceptor_act: Set to 1 if the -acceptor is activated.
 * 	donor_or_acceptor_first: It determines if the -donor or 
 * 	the -acceptor is first located in the command line.
 * 	donorTeam: The sequence of donors.
 * 	acceptorTeam: The sequence of acceptors.
 * */
extern int	donor_act;
extern int 	acceptor_act;
extern	int 	donor_or_acceptor_first;
extern	Team 	*donorTeam;
extern  Team	*acceptorTeam;

/*
 * 	All the atoms in the psf file.
 * */
extern vector<atom> table;

/*
 * 	psffile: The name of the psf file.
 * 	dcdfile: The name of the dcd file.
 * 	pdbfile: The name of the pdb file.
 * */
extern string   psffile;
extern string   dcdfile;
extern string	pdbfile;

/*
 * 	critical_distance: The critical distance in A^0 (default 5.5).
 * 	critical_percent:  The critical percent in [0,1] (default 0.05).
 * 	critical_angle:    The critical angle in degrees (default 120.0).
 * */
extern Data	critical_distance,critical_percent,critical_angle;

/*	
 *	backbone_critical_distance: The critical distance for hbonds (default 3.2).
 *	backbone_critical_percent:  The critical percent  for hbonds (default 0.05).
 *	backbone_critical_angle:    The critical angle    for hbonds (default 120.0).
 * */
extern Data	backbone_critical_distance,backbone_critical_percent,backbone_critical_angle;

/*	
 *	smart_skip: The frames to be skipped by the hbonds commands (default 10).
 *	smart_distance: The distance used in smart mode (default 2 * 3.2).
 *
 * */
extern int	smart_skip;
extern Data smart_distance;

/*
 * 	dcd: A pointer to the dcd file.
 * 	first: The first frame to be read from the dcd file.
 * 	last:  The last frame to read from the dcd file.
 * 	step:  The number of frames to be skipped in the dcd file.
 * */
extern int	first,last,step;
extern Dcd	*dcd;

extern int	backbone_flag,isbackbone;
extern string    donor_selection;
extern string    acceptor_selection;

/*	Definitions for the -watbridge parameter
 * */
extern  int	watbridge_flag;
extern  int	watbridge_level;

/*	Definitions for the -noe command
 * */
extern Data	noe_critical_distance,noe_critical_percent;


/*	Declarations for the hist command.
 * */
extern int histflag;
extern Data bindist;
extern Data binangle;
extern Data bindihe;

/*	Declarations for the smooth command
 * */
extern int smooth_flag;
extern int smooth_start;
extern int smooth_step;

/*	The list of commands activated by the user.
 * */
class Command;
extern	vector<Command*> command_list;


/*	Definitions for weights
 * */
extern int weight_flag;

extern int contact_flag;
extern Data contact_distance;
extern Data contact_percent;
extern int 	res_diff;


/*	Definitions for acorr
 * */
extern int astart;
extern int kstep;
extern int astep;


/*	Definitions for stat files.
*/
extern int statcount;
extern		int 	getAtomPos(int res_id,string atom_name);

/*
	par file
*/
extern	string	parfile;

extern	Data	diel;
extern  int	zero_angle_flag;
extern  int	degrees_angle_flag;

/*
	used in the -ref option.
*/
extern int ref;
# endif
