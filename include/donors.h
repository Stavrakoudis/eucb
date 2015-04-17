# ifndef __DONORS__H
# define __DONORS__H
# include <vector>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <util.h>
# include <psf.h>
# include <team.h>
using namespace std;

# define ISDONOR	0
# define ISHYDRO	1
# define ISACCEPTOR	2

typedef vector<int> IntVector;

typedef struct
{
	int  flag;
	int id;
	string res_name;
	string atom_name;
	string   atom_type;
}DonorT;

extern vector<int> Donor;
extern vector<int> Hydro;
extern vector<int> Acceptor;
extern vector<IntVector> DonorCon;
extern vector<IntVector> AcceptorCon;

/**/
int	isWaterDonor(int id,vector<atom> &table);
int	isProteinDonor(int id,vector<atom> &table);
int	isBackBoneDonor(int id,vector<int> &con,vector<atom> &table);
int	isSideChainDonor(int id,vector<int> &con,vector<atom> &table);
int	isWaterAcceptor(int id,vector<atom> &table);
int	isProteinAcceptor(int id,vector<atom> &table);
int	isBackBoneAcceptor(int id,vector<int> &con,vector<atom> &table);
int	isSideChainAcceptor(int id,vector<int> &con,vector<atom> &table);
/**/
void	parseDonorFile(vector<atom> &table);
void	parse_donorsrc(vector<atom> &table);
void SelectDonorAcceptor(vector<int> &CopyDonor,vector<IntVector> &CopyDonorCon,
	vector<int> &CopyAcceptor,vector<IntVector> &CopyAcceptorCon,Team *donorTeam,
	Team *acceptorTeam);

# endif
