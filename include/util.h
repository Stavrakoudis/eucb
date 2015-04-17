# ifndef __UTIL__H
# define __UTIL__H
# include <string>
# include <vector>
# include <stdio.h>
# include <stdlib.h>
# define FRAME_WIDTH		8
# define DISTANCE_WIDTH		8
# define DISTANCE_DECIMAL	3
# define ANGLE_WIDTH		7
# define ANGLE_DECIMAL		1
# define PERCENT_WIDTH		8
# define PERCENT_DECIMAL	4


using namespace std;

/*
 * 	Basic data types used in the eucb.
 * 	===============================================================
 * 	DoubleVector: A vector of Data values.
 * 	IntVector   : A vector of integer values.
 * 	Point       : A structure of the coordinates of a 3-d point.
 * 	PointVector : A vector of 3-d Points.
 * 	HistStruct  : A structure holding values used in the histograms.
 * */

typedef float  Data;

typedef vector<Data>   DoubleVector;
typedef vector<int>    IntVector;
typedef vector<IntVector>    IntVector2;
typedef struct
{
	Data x;
	Data y;
	Data z;
}Point;
typedef vector<Point> PointVector;
typedef struct
{
	Data	value;	
	Data	percent;
	int  	count;
}HistStruct;



extern	Data	getPAverage(DoubleVector &x,Data p);
/*	
 *	Return the inner product of two points (a,b).
 * */
extern Data	getDotProduct(Point a,Point b);

/*	Return the cross product of two points (a,b).
 * */
extern void	getCrossProduct(Point a,Point b,Point &c);


/*
 * 	Return the point n that is a normal to the plane defined
 * 	by p1,p2,p3.
 * */
extern void	getPlaneNormal(Point p1,Point p2,Point p3,Point &n);

/*	
 *	Return the intersection point of the vector (La,Lb) and the plane
 *	defined by the points A,B,C,D.	
 * */
extern void	getInterSectionPoint(Point La,Point Lb,Data A,Data B,Data C,Data D,Point &pt);

/*
 * 	Return the angle between the lines defined by the (La,Lb) and
 * 	(Lc,Ld).
 * */
extern Data 	getAngleBetweenLines(Point La,Point Lb,Point Lc,Point Ld);

/*	
 *	Return The four elements of a plane defined by three ponts
 *	a,b,c.
 * */
extern void     getPlaneElements(Point a,Point b,Point c,Data &A,
				  Data &B,Data &C,Data &D);

/*	
 *	Return the determinant of the matrix [x1 y1 z1].
 * */
extern Data 	getDet(Point x1,Point y1,Point z1);

/*	Return the euclidean distance of points a,b.
 * */
extern Data	getPointDistance(Point a,Point b);
extern Data getPlaneAngle(Point x1,Point y1,Point z1,
		            Point x2,Point y2,Point z2);
string getNextWord(string s,int &index,char delimiter);

/*	
 *	Return the next comma element from the current line s.
 * */
extern string	getCommaElement(string s,int &index);

/*	
 *	Return all the elements delimited by comma in the current line s.
 * */
extern void	getCommaElements(string s,vector<string> &str);

/*	
 *	Display an error message s and terminate.
 * */
extern void 	psf_error(string s);

/*	
 *	Display an error message s but do not terminate the program.
 * */
extern void 	psf_warning(string s);

/*	
 *	Return the next word (word) from the current line (line) starting
 *	from the position index.
 * */
extern int 	getword(char *line,char *word,int &index);

/*	
 *	Return in vector list all the words of the line (line).
 * */
extern void	getAllWords(char *line,vector<string> &list);

/*	
 *	Store all the words of string s to the list str.
 * */
extern void	BreakList(char *s,vector<string> &str);

/*	
 *	Parse a comma or - seperated integer list and store the 
 *	founded integers to the reslist.
 * */
extern void	ParseIntList(string str,vector<int> &reslist);

/*
 * 	Parse a sequence list.
 * */
extern void	BreakItem(string str,string &chain,vector<int> &reslist);

/*
 * 	Return the position of element y in the vector x.
 * */
extern int	isin(vector<int> x,int y);

/*	
 *	Return maximum between a,b.
 * */
extern Data	dmax(Data a,Data b);

/*
 * 	Return the minimum between a,b.
 * */
extern Data	dmin(Data a,Data b);

/*	
 *	Return minimum, maximum, average and standard deviation of
 *	the vector x.
 * */
extern void	getVectorStatistics(DoubleVector &x,Data &min,Data &max,Data &avg,Data &std);

/*
 *	Return minimum, maximum, average and standard deviation 
 *	of the angle vector x.
 * */
extern void	getAngleVectorStatistics(DoubleVector &x,Data &min,Data &max,Data &avg,Data &std);

/*	
 *	Return the percent of points in x that falls below the 
 *	value critical_value.
 * */
extern Data	getBelowPercent(DoubleVector &x,Data critical_value);

/*
 * 	Return the number of elements in x that falls below the
 * 	value critical_value.
 * */
extern int 	getBelowCount(DoubleVector &x,Data critical_value);

/*	
 *	Read the file filename and put the histogram of the data 
 *	to the file histname.
 * */
extern void	makeHist(string filename,string histname,Data minvalue,Data maxvalue,Data delta,int cols);
extern void	makeHist(DoubleVector x,string histname,Data minvalue,Data maxvalue,Data delta);
extern void	makeHist(string filename,string histname,Data minvalue,Data maxvalue,Data delta);
extern void	makeHist(DoubleVector x,vector<HistStruct> &st,Data minvalue,Data maxvalue,Data delta);

/*
 * 	Print the histogram stored in st to the file filename.
 * */
extern	void	PrintHist(string filename,vector<HistStruct> &st);

/*	
 *	Print the angle histogram stored in st to the file filename.
 * */
extern	void	PrintAngleHist(string filename,vector<HistStruct> &st);

/*	
 *	Make smooth data for the vector f and x.
 * */
extern void	makeSmooth(IntVector &f,DoubleVector &x,IntVector &sf,DoubleVector &sx,int start,int step);

/*	
 *	Make smooth data for the vector f and x (x is a vector of angles).
 * */
extern void	makeAngleSmooth(IntVector &f,DoubleVector &x,IntVector &sf,DoubleVector &sx,int start,int step);

/*	Calculate and return the dihedral angle of four atoms.
 * */
extern Data torsion( Data x1, Data y1, Data z1,
                Data x2, Data y2, Data z2,
                Data x3, Data y3, Data z3,
               Data x4, Data y4, Data z4);

/*	Print numbers with format.
 *	The final outcome is a string.
 * */
extern	 string	printNumber(int x);
extern	 string	printNumber0(int x,int space);
extern	 string printNumber(int x,int space);
extern	 string printNumber(Data x);
extern	 string printNumber(Data x,int space,int decimal);
extern	 string printPercent(Data x);
extern	 string	printString(string s,int tab);

/*	Print specific type of numbers 
 *	using the above formats.
 * */
extern	 string printFrame(int f);
extern   string printDistance(Data x);
extern	 string printAngle(Data a);

/*	Print strings to the current output
 *	or to the screen.
 * */
extern	 void	Print(string s);
extern	 void	PrintLine(string s);
extern	 void	Print(FILE *fp,string s);
extern	 void	PrintLine(FILE *fp,string s);
/*	Print	vectors to files
 * */
extern   void	PrintDistances(string filename,IntVector &frame,DoubleVector &distance);
extern	 void	PrintAngles(string filename,IntVector &frame,DoubleVector &angle);
extern   void	PrintDistAndAngle(string filename,IntVector &frame,DoubleVector &distance,DoubleVector &angle);

extern 	 string	printFrameHeader();
extern	 string printDistanceHeader();
extern	 string printAngleHeader();
extern	 string printTorsionHeader();

extern 	 string	printFrameHeader(string s);
extern	 string printDistanceHeader(string s);
extern	 string printAngleHeader(string s);
extern	 string printTorsionHeader(string s);
extern   string printPercentHeader(string s);


/*
	Print a message relative to writing.
*/
extern	 string	WriteError(string filename);
extern	 string SeqError(string command);


/*	Print the correlation matrix
 * */
extern	void	PrintAcor(string filename,IntVector &F,IntVector Acor,int kstep,int astep);


extern void	getVectorStatistics(int icount,IntVector &start,IntVector &last,IntVector &F,
		DoubleVector &x,DoubleVector &min,DoubleVector &max,DoubleVector &avg,DoubleVector &std);
extern void	getAngleVectorStatistics(int icount,IntVector &start,IntVector &last,IntVector &F,
		DoubleVector &x,DoubleVector &min,DoubleVector &max,DoubleVector &avg,DoubleVector &std);
extern  double	getMinimum(DoubleVector &x,int icount);
extern  double	getMaximum(DoubleVector &x,int icount);
extern	double	getAverage(DoubleVector &x,int icount);
extern	double	getStd(DoubleVector &x,int icount);

extern	void	PrintStatFile(string filename,string mode,int icount,string header,IntVector &F,DoubleVector &x);
extern	void	PrintAngleStatFile(string filename,string mode,int icount,string header,IntVector &F,DoubleVector &x);
extern	void	PrintIntStatFile(string filename,string mode,int icount,string startheader,
		string header,IntVector &F,IntVector2 &x);
extern	void	SplitTable(int icount,IntVector &start,IntVector &last,IntVector &F,IntVector &x,IntVector &res);
extern	void	printVersion();
extern Data deg2rad(Data x);
extern Data rad2deg(Data x);
# define __UTIL__H
# endif
