# ifndef __DCD__H
# define __DCD__H
# include <util.h>
# include <vector>
using namespace std;
# define REQUEST_DISTANCE		1
# define REQUEST_ANGLE			2
# define REQUEST_DANGLE			3
# define REQUEST_ANGLE_DISTANCE		4

#define	HFIT				1
#define BBFIT				2
#define CAFIT				3
typedef struct
{
	int x1;
	int x2;
	IntVector water;
	DoubleVector distance;
	DoubleVector angle;
}Pair;


typedef vector<DoubleVector> DoubleVector2;
typedef vector<DoubleVector2> DoubleVector3;

typedef vector<IntVector> IntVector2;

typedef struct
{
        Data dist1;
        Data dist2;
        int id1;
        int id2;
        int id3;
        int     pos_id1;
        int     neg_id;
        int     pos_id2;
}PosNegPosStruct;


typedef vector<PosNegPosStruct>    PosNegPosTable;

typedef struct
{
	string res_name;
	string atom_name;
}ComplexSaltStruct;

typedef struct
{
        Data  dist;
        int     id1;
        int     id2;
}PosNegInfo;
typedef vector<PosNegInfo> PosTable2;

typedef struct
{
	Data dist1;
	Data dist2;
	int id1;
	int id2;
	int id3;
	int	neg_id1;
	int	pos_id;
	int	neg_id2;
}NegPosNegStruct;
typedef vector<NegPosNegStruct>    NegPosNegTable;

typedef struct
{
	int donor_id;
	IntVector Acceptor;
	vector<DoubleVector> dist;
	vector<DoubleVector> angle;
	IntVector Hydro;
	IntVector HydroFlag;
	vector<DoubleVector> HydroAcc;
}HbondStruct;

typedef struct
{
	int acceptor_id;
	IntVector Donor;
	vector<DoubleVector> dist;
	vector<DoubleVector> angle;
	vector<IntVector> Hydro;
	IntVector HydroFlag;
	vector<DoubleVector2> HydroDonor;
}AcceptorStruct;



typedef struct
{
	int id1;
	int id2;
	int id3;
}stack_table;

typedef struct
{
	int atom1_id1;
	int atom1_id2;
	int atom1_id3;
	int atom2_id1;
	int atom2_id2;
	int atom2_id3;
	DoubleVector D;
	DoubleVector A;
	Data percent;
}StackAtom;

typedef struct
{
	int npos;
	int hpos;
}NHStruct;

typedef struct
{
	int npos;
	IntVector hpos;
}AroposStruct;

typedef struct
{
	int atom1;
	int atom2;
	DoubleVector D;
}NoeStruct;

typedef struct
{
	int atom1;
	int atom2;
	int atom3;
	DoubleVector D;
}NoeAngleStruct;

typedef struct
{
	int atom1;
	int atom2;
	int atom3;
	int atom4;
	DoubleVector D;
}NoeDihedralStruct;

typedef struct
{
	int dist_atom1;
	int dist_atom2;
	int dist_atom3;
	int dist_atom4;
	int dihe_atom1;
	int dihe_atom2;
	int dihe_atom3;
	int dihe_atom4;
	int dihe_atom5;
	int dihe_atom6;
	int dihe_atom7;
	int dihe_atom8;
	DoubleVector D1;
	DoubleVector D2;
	DoubleVector A1;
	DoubleVector A2;
}PdoStruct;

typedef struct
{
	int	atom1;
	int	atom2;
	IntVector	Frame;
	DoubleVector	Distance1;
	DoubleVector	Distance2;
	IntVector	Water;
}HydroSolveStruct;

typedef struct
{
	int		water;
	IntVector	list;
	IntVector	count;
	vector<IntVector>	atomList;
}IwcnStruct;

typedef struct
{
	int	atom1;
	int	atom2;
	int	count;
}CountFramesStruct;
class Dcd
{
	private:
		/*
 		Variables:
		========================================================
		posCAs: The index of atoms in the dcd file.
		CAs:    The position (x,y,z) of atoms in the dcd file.
		dcd:	The file descriptor of the dcd file.
		natoms: The number of atoms in the dcd file.
		headerframes: The number of frames in the dcd file.
		frameStart: The first frame to be read from the file.
		frameEnd: The last frame to be read from the file.
		frameStep: The frames to skipped in every read.
		Methods: 
		========================================================
		Dcd(): The constructor of the class. It opens the file.
		build_CAs(): Reads atoms from the dcd file at the 
		    current frame.
		skipFrame(): Skip frames from the dcd file.
		readCAs(): Read atoms' position from the file.
		dihedral(): Calculate and return a dihedral angle.
		getDistance(): Calculate and return an euclidan distance.
		getAngle(): Calculate and return an angle.
		distAndAngle(): Calculate and return both distances and 
		  angles.
		setFrameStart(): Set the first frame to be read.
		setFrameEnd(): Set the last frame to be read.
		setFrameStep(): Set the frames to be skipped.
		getFrameStart(): Return the first frame to be read.
		getFrameEnd(): Return the last frame to be read.
		getFrameStep(): Return the number of frames skipped.
		getnatoms(): Return the number of atoms in the file.
		getframes(): Return the number of frames in the file.
		printDistance(): Return the distance between two atoms.
		printAngle(): Return the angle between three atoms.
		printDAngle(): Return the dihedral of four atoms.
		hbondsCheck(): Checks for hbonds conectivity.
		hbonds(): Return the hydrogen bond.
		acceptors(): Return the acceptors.
		countAtoms(): Return the number of atoms for a property.
		printRMSF(): Calculate and return the rmsf.
		getStackAverages(): Calculate and return stack averages.
		getAronh(): Calculate and return the NH facility.
		getAropos(): Calculate and return the the 
		aropos facility.
		getNoe(): Calculate and return the distances between 
		pairs of atoms.
		getNoeAngle(): Calculate and return the angle between
		3-tuples of atoms.
		getNoeDihedral(): Calculate and return the dihedral 
		angle between 4-tuples of atoms.
		getNoeCenter(): Calculate and return the distance 
		between the centers of two atom's lists.
		getRmsd(): Calculate and return the rmsd for a list 
		of atoms.
		~Dcd(): Close the file and free the allocated memory.
 		* */

		Data cell_a,cell_gamma,cell_b,cell_beta,cell_alpha,cell_c;
		IntVector posCAs;
		DoubleVector2 CAs;
		float *dcd_frame,*wdcd_frame;
		char filename[1024];
		char dcd_head[92];
		char dcd_head2[16];
		char title[81];
		unsigned int auint,bytes_per_set,wbytes_per_set;
		int have_cell;
		int cell_offset;
		int dcd;
		unsigned int natoms;
		int headerframes;
		int	frameStart,frameEnd,frameStep;

		void	build_CAs(int offset);
		int	skipFrame(int &frame);
		int	readCAs();
		Data	dihedral(DoubleVector x1,DoubleVector x2,DoubleVector x3,DoubleVector x4);
		Data	dihedral(int x1,int x2,int x3,int x4);
		Data	dihedral(Point x1,Point x2,Point x3,Point x4);
		Data	getDistance(Point a,Point b);
		Data	getDistance(int  DISTFIRST,int DISTSECOND);
		Data	getAngle(int DISTFIRST,int DISTSECOND,int DISTTHIRD);
		Data	getAngle(Point p1,Point p2,Point p3);
		Point	makePoint(int index);
		Point	makeCenter(IntVector &list,int wflag);	
	public:
		Dcd(char *s);
		int	nframes();
		void	distAndAngle(int DISTFIRST,int DISTSECOND,int DISTTHIRD,Data &dist,Data &angle);
		void	setFrameStart(int s);
		void	setFrameEnd(int s);
		void	setFrameStep(int s);
		const 	int getFrameStart() const;
		const	int getFrameEnd() const;
		const	int getFrameStep() const;
		const	unsigned int getnatoms() const;
		const	int getframes() const;
		void	hbondsCheck(Data cdist,Data cangle,Data cpercent,
			IntVector &X1,vector<IntVector> &X2,
			IntVector &X3,vector<IntVector> &Flag);
		void	hbonds(IntVector &f,vector<HbondStruct> &st);

		void	acceptors(IntVector &f,vector<AcceptorStruct> &st,IntVector &Donor);

		Data	countAtoms(int start,Data distance,vector<int> &list,vector<Data> &vframe);
		void	printRMSF(vector<int> &res_number,vector<Data> &rmsf);
		void	getStackAverages(IntVector &F,vector<stack_table> &stable,
				vector<StackAtom> &stack_atom);
		void	getPosNeg(IntVector &f,IntVector &pos,IntVector &neg,
				vector<ComplexSaltStruct> &poscsalt,
				vector<ComplexSaltStruct> &negcsalt,
				vector<PosTable2> &PosTable);

		void	getPosNegPos(IntVector &f,IntVector &pos,IntVector &neg,
					vector<ComplexSaltStruct> &poscsalt,
					vector<ComplexSaltStruct> &negcsalt,
					vector<PosNegPosTable> &PosTable);

		void	getNegPosNeg(IntVector &f,IntVector &pos,IntVector &neg,
					vector<ComplexSaltStruct> &poscsalt,
					vector<ComplexSaltStruct> &negcsalt,
					vector<NegPosNegTable> &PosTable);

		void	getGroupDistance(IntVector2 &group1,IntVector2 &group2,
				IntVector &f,DoubleVector2 &d,
				IntVector2 &atom1,IntVector2 &atom2);

		void	getAronh(IntVector &f,vector<DoubleVector> &Dist,
				vector<DoubleVector> &Dist2,
				vector<DoubleVector> &Angle,
				vector<stack_table> &stable,vector<NHStruct> &nhpos,
			Data stack_distance,Data stack_angle);
		void	getAropos(IntVector &f,vector<DoubleVector> &Dist,vector<DoubleVector> &Angle,
				vector<stack_table> &stable,vector<AroposStruct> &nhpos,
			Data stack_distance,Data stack_angle);
		void	getNoe(IntVector &f,vector<NoeStruct> &noe);
		void	getNoeAngle(IntVector &f,vector<NoeAngleStruct> &noe);
		void	getNoeDihedral(IntVector &f,vector<NoeDihedralStruct> &noe);
		void	getNoeCenter(IntVector &atom1,IntVector &atom2,IntVector &F,DoubleVector &D);
		void	getNoeCenterAngle(IntVector &atom1,IntVector &atom2,
			IntVector &atom3,IntVector &F,DoubleVector &D);
		void	getNoeCenterDihedral(IntVector &atom1,IntVector &atom2,
			IntVector &atom3,IntVector &atom4,IntVector &F,DoubleVector &D);
		void	getRmsd(IntVector &AtomPos,PointVector &pdbpos,IntVector &F,DoubleVector &D);
		void	getRmsd(IntVector &AtomPos,IntVector &F,DoubleVector2 &D);
		void	fetchAtoms(int frame,vector<DoubleVector> &D);
		void	getPdo(IntVector &F,vector<PdoStruct> &Pdo);
		void	getHydroSolve(IntVector &waterlist,IntVector &F,
			vector<HydroSolveStruct> &st,Data distance,Data angle);
		void	getCenterNear(IntVector &F,IntVector &atom1,IntVector &atom2,
				IntVector &Index,DoubleVector &distance);
		void	getIwcn(double d,IntVector &F,vector<IwcnStruct> &iwcn);
		void	countFrames(double d,IntVector &F,vector<CountFramesStruct> &count);
		void	getNear(IntVector &F,IntVector2 &List,DoubleVector2 &Dist);
		void	getNear(IntVector &F,IntVector2 &List1,IntVector &List2,DoubleVector2 &Dist,
			IntVector &Vote);
		void	getEntropyData(IntVector &F,DoubleVector3 &D,
				IntVector &posAtom,int N,
				int fitflag,PointVector &pdbpos);
		void	getAvgPos(IntVector &F,IntVector &posAtom,
			DoubleVector2 &D,int N,
				int fitflag,PointVector &pdbpos);
		void	writeDcd(string fname,IntVector &posAtom,int fit);
		void	getEnergy(IntVector &F,DoubleVector2 &E,IntVector &posAtom,vector<string> &list);
		void	getFitted(int fitflag,IntVector &PosAtom,PointVector &pdbpos,
				DoubleVector2 &orig,DoubleVector2 &dest);
		void	getSasa(IntVector &F,IntVector &posAtom,Data srad,DoubleVector &D);
		void	getDisu(IntVector &F,IntVector &Index1,IntVector &Index2,IntVector &Index3,
				IntVector &Index4,DoubleVector2 &D,DoubleVector2 &A);
		~Dcd();
};

# endif
