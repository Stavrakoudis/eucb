# ifndef __PARFILE__H
# define __PARFILE__H
# include <global.h>
# include <command.h>

typedef struct
{
	string atomtype1;
	string atomtype2;
	Data	kb;
	Data 	b0;
}ParBondStruct;

typedef struct
{
	string atomtype1;
	string atomtype2;
	string atomtype3;
	Data	ktheta;
	Data	theta0;
	Data	kub;
	Data	s0;
}ParAngleStruct;

typedef struct
{
	string	atomtype1;
	string	atomtype2;
	string	atomtype3;
	string	atomtype4;
	Data	kchi;
	int	n;
	Data	delta;
}ParDihedralStruct;

typedef	struct
{
	string	atomtype1;
	string	atomtype2;
	string	atomtype3;
	string	atomtype4;
	Data	kpsi;
	Data	psi0;
}ParImproperStruct;

typedef struct
{
	string	atomtype1;
	Data	epsilon;
	Data	rmin;
}ParVdwStruct;

class ParFile
{
	private:
		vector<ParBondStruct> bondstruct;
		vector<ParAngleStruct> anglestruct;
		vector<ParDihedralStruct> dihedralstruct;
		vector<ParImproperStruct> improperstruct;
		vector<ParVdwStruct>	  vdwstruct;
	public:
		ParFile(string s);
		Data	getBondEnergy(string a,string b,Data B);
		Data	getAngleEnergy(string a,string b,string c,Data theta,Data d);
		Data	getDihedralEnergy(string a,string b,string c,string d,Data chi);
		Data	getImproperEnergy(string a,string b,string c,string d,Data psi);
		Data	getVdwEnergy(string a,string b,Data d);
		~ParFile();
};
# endif
