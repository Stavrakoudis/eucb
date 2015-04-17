# ifndef __COMMAND_ANGLE__H
# define __COMMAND_ANGLE__H
# include <command.h>

class	CommandAngle :public Command
{
	private:
		/*
 		Variables:
		=======================================================
		phi_flag: Enable or disable the phi angles.
		psi_flag: Enable or disable the psi angles.
		omega_flag: Enable or disable the omega angles.
		chi1_flag:  Enable or disable the chi1 angles.
		JHNHA_a,JHNHA_b,JHNHA_c: Used in the JHNHA statistic.
		Methods:
		=======================================================
		CommandAngle(): The constructor of the class.
		setAngleList(): Specify which angles do we want.
		setJHNHA(): Set the parameters for the JHNHA statistic.
		Run(): Run the torsion command.
 		* */
		int	phi_flag;
		int	psi_flag;
		int	omega_flag;
		int	chi1_flag;
		int	chi2_flag;
		int	chi3_flag;
		int	chi4_flag;
		int	zero_flag;
		int	deg_flag;
		int	firstAminoAcid;
		int	lastAminoAcid;
		int	cycleFlag,pos_N,pos_C,pos_CA;
		Data	JHNHA_a;
		Data  JHNHA_b;
		Data  JHNHA_c;
		vector<NoeDihedralStruct> NoeDihedral;
		int	getAtomsChi1(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		int	getAtomsChi2(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		int	getAtomsChi3(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		int	getAtomsChi4(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		int	getAtomsOmega(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		int	getAtomsPhi(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		int	getAtomsPsi(int ipos,IntVector &posAtom,int &x1,int &x2,int &x3,int &x4);
		void	newDihedral(int x1,int x2,int x3,int x4,string angle);
	public:
		CommandAngle();
		void	setAngleList(vector<string> &s);
		void	setJHNHA(Data a,Data b,Data c);
		virtual void Run();
};
# endif
