# include <util.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>
# include <sstream>
# include <iomanip>
# include <iostream>
using namespace std;
double	getMinimum(DoubleVector &x,int icount)
{
	Data mmin=x[0];
	for(int i=1;i<icount;i++) if(x[i]<mmin) mmin=x[i];
	return mmin;
}

double	getMaximum(DoubleVector &x,int icount)
{
	Data mmax=x[0];
	for(int i=1;i<icount;i++) if(x[i]>mmax) mmax=x[i];
	return mmax;
}

double	getAverage(DoubleVector &x,int icount)
{
	Data mavg=0.0;
	for(int i=0;i<icount;i++) mavg+=x[i];
	return mavg/icount;
}

double	getStd(DoubleVector &x,int icount)
{
	double mavg=getAverage(x,icount);
	double mstd=0.0;
	for(int i=0;i<icount;i++) mstd+=pow(x[i]-mavg,2.0);
	return mstd/(icount * (icount-1));
}

void	makeAngleSmooth(IntVector &f,DoubleVector &x,IntVector &sf,DoubleVector &sx,int start,int step)
{
	sf.resize(0);
	sx.resize(0);
	Data average=0.0;
	int istart=0;
	if(step>=0)
	{
		Data sinaverage=0.0;
		Data cosaverage=0.0;
		for(int i=0;i<start;i++)
		{
			Data angle = x[i] * M_PI/180.0;
			sinaverage+=sin(angle);
			cosaverage+=cos(angle);
		}
		average=atan2(sinaverage/start,cosaverage/start);
		average=average*180.0/M_PI;

		sx.push_back(average);

		sf.push_back(f[start-1]);	
		istart=step;
	}
	for(int frame_count=istart;frame_count<f.size();frame_count+=start)
	{
		average=0.0;
		int icount=0;
		int i;
		Data sinaverage=0.0;
		Data cosaverage=0.0;
		for(i=frame_count;i<f.size() && i<frame_count+start;i++) 
		{
			icount++;
			Data angle = x[i] * M_PI/180.0;
			sinaverage+=sin(angle);
			cosaverage+=cos(angle);
		}
		average=atan2(sinaverage/icount,cosaverage/icount);
		average=average*180.0/M_PI;
		sx.push_back(average);
		sf.push_back(f[i-1]);
	}
}

void	makeSmooth(IntVector &f,DoubleVector &x,IntVector &sf,DoubleVector &sx,int start,int step)
{
	sf.resize(0);
	sx.resize(0);
	Data average=0.0;
	int istart=0;
	if(step>=0)
	{
		for(int i=0;i<start;i++) average+=x[i];
		sx.push_back(average/start);
		sf.push_back(f[start-1]);	
		istart=step;
	}
	for(int frame_count=istart;frame_count<f.size();frame_count+=start)
	{
		average=0.0;
		int icount=0;
		int i;
		for(i=frame_count;i<f.size() && i<frame_count+start;i++) {icount++;average+=x[i];}
		sx.push_back(average/icount);
		sf.push_back(f[i-1]);
	}
}

Data	getPAverage(DoubleVector &x,Data p)
{
	Data sum=0.0;
	for(int i=0;i<x.size();i++) sum+=pow(x[i],-p);
	sum = sum /x.size();
	return pow((double)sum,(double)-1.0/p);
}


void	getPlaneElements(Point a,Point b,Point c,Data &A,
				  Data &B,Data &C,Data &D)
{
	Point p1,p2,p3;
	p1.x=a.x;
	p1.y=b.x;
	p1.z=c.x;

	p2.x=a.y;
	p2.y=b.y;
	p2.z=c.y;

	p3.x=a.z;
	p3.y=b.z;
	p3.z=c.z;
	Data Det=getDet(p1,p2,p3);
	Point n;
	getPlaneNormal(p1,p2,p3,n);
	Data d=getDotProduct(p1,n);
	
	Point p11;
	p11.x=p11.y=p11.z=1.0;
	A=d/Det*getDet(p11,p2,p3);
	B=d/Det*getDet(p1,p11,p3);
	C=d/Det*getDet(p1,p2,p11);	
	D=d;
}

Data	getDotProduct(Point a,Point b)
{
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

void	getCrossProduct(Point a,Point b,Point &c)
{
	c.x=a.y*b.z-a.z*b.y;
	c.y=a.z*b.x-a.x*b.z;
	c.z=a.x*b.y-a.y*b.x;
}

void	getPlaneNormal(Point p1,Point p2,Point p3,Point &n)
{
	Point pt1,pt2;
	pt1.x=p2.x-p1.x;
	pt1.y=p2.y-p1.y;
	pt1.z=p2.z-p1.z;
	
	pt2.x=p3.x-p1.x;
	pt2.y=p3.y-p1.y;
	pt2.z=p3.z-p1.z;
	getCrossProduct(pt1,pt2,n);
}

void	getInterSectionPoint(Point La,Point Lb,Data A,Data B,Data C,Data D,Point &pt)
{
	Data t;
	t=(D-A*La.x-B*La.y-C*La.z)/(A*(Lb.x-La.x)+B*(Lb.y-La.y)+C*(Lb.z-La.z));
	pt.x=La.x+(Lb.x-La.x)*t;
	pt.y=La.y+(Lb.y-La.y)*t;
	pt.z=La.z+(Lb.z-La.z)*t;
}

Data 	getAngleBetweenLines(Point La,Point Lb,Point Lc,Point Ld)
{
	Data sum=0.0;
	Data dx1=Lb.x-La.x;
	Data dy1=Lb.y-La.y;
	Data dz1=Lb.z-La.z;
	Data dx2=Ld.x-Lc.x;
	Data dy2=Ld.y-Lc.y;
	Data dz2=Ld.z-Lc.z;
	sum=dx1*dx2+dy1*dy2+dz1*dz2;
	Data a=asin(sum/(getPointDistance(La,Lb)*getPointDistance(Lc,Ld)));
	return 90.0-a*180.0/M_PI;
}

#define VOP_3D_COORDS_CROSS_PRODUCT(TX,TY,TZ,UX,UY,UZ,VX,VY,VZ) \
	  TX = (UY * VZ) - (UZ * VY); \
  TY = (UZ * VX) - (UX * VZ); \
  TZ = (UX * VY) - (UY * VX)

Data torsion( Data x1, Data y1, Data z1,
                Data x2, Data y2, Data z2,
                Data x3, Data y3, Data z3,
               Data x4, Data y4, Data z4)
{
	Data Lx, Ly, Lz, Lnorm;
	Data Rx, Ry, Rz, Rnorm;
	Data Sx, Sy, Sz;
	Data angle;
	VOP_3D_COORDS_CROSS_PRODUCT(     Lx,      Ly,      Lz,
		                       (x2-x1), (y2-y1), (z2-z1),
		                       (x3-x2), (y3-y2), (z3-z2));
 	VOP_3D_COORDS_CROSS_PRODUCT(     Rx,      Ry,      Rz,
                              (x4-x3), (y4-y3), (z4-z3),
                            (x2-x3), (y2-y3), (z2-z3));
	Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
	Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

	VOP_3D_COORDS_CROSS_PRODUCT(Sx, Sy, Sz,
                            Lx, Ly, Lz,
                            Rx, Ry, Rz);
	angle = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm);

	if ( angle > 1.0 ) angle = 1.0;
	if ( angle < -1.0 ) angle = -1.0;
	angle = acos( angle );
	angle = 180.0 * angle / M_PI;
	if ( (Sx * (x3-x2) + Sy * (y3-y2) + Sz * (z3-z2)) < 0 ) angle = -angle;
	return(angle);
}

string getNextWord(string s,int &index,char delimiter)
{
	string ret="";
	while(index<s.size() && s[index]!=delimiter)
	{
		ret = ret + s[index];
		index++;
	}
	index++;
	return ret;
}

Data getPlaneAngle(Point x1,Point y1,Point z1,
		            Point x2,Point y2,Point z2)
{
	Point o;o.x=o.y=o.z=1.0;
	Data d1=getDet(x1,y1,z1);
	Data d2=getDet(x2,y2,z2);
	Data a1=getDet(o,y1,z1)/d1;
	Data b1=getDet(x1,o,z1)/d1;
	Data c1=getDet(x1,y1,o)/d1;
	Data a2=getDet(o,y2,z2)/d2;
	Data b2=getDet(x2,o,z2)/d2;
	Data c2=getDet(x2,y2,o)/d2;
	Data meanx1=(x1.x+x1.y+x1.z)/3.0;
	Data meanx2=(x2.x+x2.y+x2.z)/3.0;
	Data meany1=(y1.x+y1.y+y1.z)/3.0;
	Data meany2=(y2.x+y2.y+y2.z)/3.0;
	Data meanz1=(z1.x+z1.y+z1.z)/3.0;
	Data meanz2=(z2.x+z2.y+z2.z)/3.0;
	Data d=sqrt((meanx1-meanx2)*(meanx1-meanx2)+(meany1-meany2)*(meany1-meany2)+
		(meanz1-meanz2)*(meanz1-meanz2));
	Data p=(a1*a2+b1*b2+c1*c2) / ( sqrt(a1*a1 + b1*b1 + c1*c1) * sqrt(a2*a2 + b2*b2 + c2*c2) );
	p=acos(p);
	p = fabs(p*180/M_PI);
        Data p1= dmin(p,180-p);
	return p1;
}

Data	getPointDistance(Point a,Point b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

/*	Return the determinat of the matrix
 *	[x1 y1 z1]
 * */
Data getDet(Point x1,Point y1,Point z1)
{
	return x1.x*(y1.y*z1.z-y1.z*z1.x)-y1.x*(x1.y*z1.z-x1.z*z1.y)+z1.x*(x1.y*y1.z-x1.z*y1.y);
}

string	getCommaElement(string s,int &index)
{
	if(index>=s.size()) return "";
	string s1="";
	if(s[index]==',') index++;
	while(index<s.size() && s[index]!=',')
	{
		s1=s1+s[index];
		index++;
	}
	return s1;
}

void	getCommaElements(string s,vector<string> &str)
{
	str.resize(0);
	string ss;
	int index=0;
	while((ss=getCommaElement(s,index))!="")
		str.push_back(ss);
}

/*	This file contains various utility functions.
 * */

/*	Return the position of the item y in vector x.
 * */
int	isin(vector<int> x,int y)
{
	for(int i=0;i<x.size();i++)
		if(x[i]==y) return i;
	return -1;
}

string	WriteError(string filename)
{
	return "Can not open "+filename+" for writing.";
}

string	SeqError(string command)
{
	return "You must specify sequnce for command "+command;
}

/*	Print an error message to the user and terminate the program.
 * */
void psf_error(string s)
{
	fprintf(stderr,"\033[31m\033[1mFATAL ERROR OCCURED..\033[0m\n");
	fprintf(stderr,"\033[31m\033[1m%s\033[0m\n",s.c_str());
	exit(EXIT_FAILURE);
}

void	psf_warning(string s)
{
	fprintf(stderr,"\033[31m\033[1m%s\033[0m\n",s.c_str());
}

/*	Break the line in words and return the next word.
 * */
int getword(char *line,char *word,int &index)
{
	if(index==(int)strlen(line)) return 0;
	while(line[index]==' ')index++;
	int icount=0;
	while(line[index]!=' ' && index<strlen(line))
	{
		word[icount]=line[index];
		index++;
		icount++;
	}
	word[icount]=0;
	return 1;
}

void	getAllWords(char *line,vector<string> &list)
{
	char word[1024];
	int index=0;
	while(getword(line,word,index))
	{
		list.push_back(word);
	}
}

/*	
 * */
void	BreakList(char *s,vector<string> &str)
{
	string ss="";
	int index=0;
	while(1)
	{
		if(index==strlen(s)) break;
		if(s[index]==',')
		{
			str.push_back(ss);
			ss="";
		}
		else
			ss=ss+s[index];
		index++;
	}
	if(ss!="")
	{
		str.push_back(ss);
	}
}

/*	Return a vector with all integers in a string. It is used in the -seq option.
 * */
void	ParseIntList(string str,vector<int> &reslist)
{
	int i=0;
	string s="";
	int lastint=-1;
	while(1)
	{
		if(i==str.size()) break;
		if(str[i]=='-') 
		{
			lastint=atoi(s.c_str());
			s="";
		}
		else
		{
			s=s+str[i];
		}
		i++;
	}
	if(lastint==-1)
	{
		reslist.resize(1);
		reslist[0]=atoi(s.c_str());
	}
	else
	{
		int newint=atoi(s.c_str());
		reslist.resize(newint-lastint+1);
		for(int i=lastint;i<=newint;i++)
		{
			reslist[i-lastint]=i;
		}
	}
}

void	BreakItem(string str,string &chain,vector<int> &reslist)
{
	chain="";
	reslist.resize(0);
	string rest="";
	if(str.find(':')!=-1)
	{
		int i;
		for(i=0;str[i]!=':';i++)
			chain=chain+str[i];
		i++;
		for(;i<str.size();i++)
			rest=rest+str[i];
		ParseIntList(rest,reslist);
		return;
	}
	if(isdigit(str[0]))
		ParseIntList(str,reslist);
	else
		chain = str;
}

/*
 * 	Return the maximum value between
 * 	two Data numbers.
**/
Data dmax(Data a,Data b)
{
	return a>b?a:b;
}

/*
 * 	Return the minimum value between
 * 	two Data numbers.
**/
Data dmin(Data a,Data b)
{
	return a<b?a:b;
}

int	getBelowCount(vector<Data> &x,Data critical_value)
{
	int icount = 0;
	for(int i=0;i<x.size();i++)
		if(x[i]<=critical_value) icount++;
	return icount;
}
/*	Return the percent of x - values (in 0. format) that falls below
 *	the critical_value.
 * */
Data	getBelowPercent(DoubleVector &x,Data critical_value)
{
	int icount = 0;
	for(int i=0;i<x.size();i++)
		if(x[i]<=critical_value) icount++;
	return icount*1.0/x.size();
}

void getVectorStatistics(DoubleVector &x,Data &min,Data &max,Data &avg,Data &std)
{
	min = x[0];
	max = x[0];
	avg = 0.0;
	std = 0.0;	
	Data x1=0.0;
	Data x2=0.0;
	for(int i=0;i<x.size();i++)
	{
		if(x[i]<min) min = x[i];
		if(x[i]>max) max = x[i];
		avg = avg+x[i];
		x1+=x[i];
		x2+=x[i]*x[i];
	}
	avg/=x.size();
	std=sqrt(x2/x.size()-x1/x.size()*x1/x.size());
}

/*	Create histogram for the values of filename. The value
*	of delta determines the gap. We assume that the input 
*	file has cols columns.
 * */
void	makeHist(string filename,string histname,Data minvalue,Data maxvalue,Data delta,int cols)
{
	minvalue = floor(minvalue);
	maxvalue = ceil(maxvalue);
	extern int histflag;
	if(histflag==0) return;
	int nsize = (int)((maxvalue-minvalue)/delta)+1;
	typedef vector<int> intv;
	vector<intv> hist;
	hist.resize(cols);
	for(int i=0;i<cols;i++) hist[i].resize(nsize);

	for(int j=0;j<cols;j++)
	for(int i=0;i<nsize;i++) hist[j][i]=0;

	FILE *fp=fopen(filename.c_str(),"r");
	int icount=0;
	vector<Data> value;
	value.resize(cols);

	while(1)
	{
		int frame;
		int k=fscanf(fp,"%d",&frame);
		if(k<=0) break;
		int break_flag = 0;
		for(int i=0;i<cols;i++)
		{
			int k=fscanf(fp,"%f",&value[i]);
			if(k<=0) {break_flag=1;break;}
		}
		if(break_flag) break;
		icount++;
		for(int i=0;i<nsize;i++)
		{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			for(int j=0;j<cols;j++)
			{
				if(i==nsize-1 && fabs(value[j])==180.0) hist[j][i]++;
				else
				if(value[j]>a && value[j]<=b) hist[j][i]++;
			}
		}
	}
	fclose(fp);
	fp=fopen(histname.c_str(),"w");
	for(int i=0;i<nsize;i++)
	{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(b>maxvalue)
			{
				fprintf(fp,"%15.1lf ",a);
				for(int j=0;j<cols;j++)
					fprintf(fp,"%6.4lf %7d",hist[j][i]*1.0/icount,hist[j][i]);
			}
			else
			{
				fprintf(fp,"%15.1lf ",a+(b-a)/2.0);
				for(int j=0;j<cols;j++)
					fprintf(fp,"%6.4lf %7d",hist[j][i]*1.0/icount,hist[j][i]);
			}
			fprintf(fp,"\n");
	}
	fclose(fp);
}

void	makeHist(DoubleVector x,vector<HistStruct> &st,Data minvalue,Data maxvalue,Data delta)
{
	minvalue = floor(minvalue);
	maxvalue = ceil(maxvalue);
	extern int histflag;
	if(histflag==0) return;
	int nsize = (int)((maxvalue-minvalue)/delta)+1;
	st.resize(nsize);
	for(int i=0;i<nsize;i++) 
	{
		st[i].count = 0;
		st[i].percent = 0.0;
		st[i].value = 0.0;
	}
	int icount=x.size();
	for(int k=0;k<icount;k++)
	{
		Data value;
		value = x[k];
		for(int i=0;i<nsize;i++)
		{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(i==nsize-1 && fabs(value)==180.0) st[i].count++;
			else
			if(value>a && value<=b) st[i].count++;
		}
	}
	for(int i=0;i<nsize;i++)
	{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(b>maxvalue)
			{
				st[i].value = a;
				st[i].percent=st[i].count*1.0/icount;
			}
			else
			{
				st[i].value = a+(b-a)/2.0;
				st[i].percent=st[i].count*1.0/icount;
			}
	}
}

void	makeHist(DoubleVector x,string histname,Data minvalue,Data maxvalue,Data delta)
{
	minvalue = floor(minvalue);
	maxvalue = ceil(maxvalue);
	extern int histflag;
	if(histflag==0) return;
	int nsize = (int)((maxvalue-minvalue)/delta)+1;
	vector<int> hist;
	hist.resize(nsize);
	for(int i=0;i<nsize;i++) hist[i]=0;
	int icount=x.size();
	for(int k=0;k<icount;k++)
	{
		Data value;
		value = x[k];
		for(int i=0;i<nsize;i++)
		{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(i==nsize-1 && fabs(value)==180.0) hist[i]++;
			else
			if(value>a && value<=b) hist[i]++;
		}
	}
	FILE *fp=fopen(histname.c_str(),"w");
	for(int i=0;i<nsize;i++)
	{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(b>maxvalue)
		fprintf(fp,"%15.5lf %15.5lf %6d\n",a,hist[i]*1.0/icount,hist[i]);
		else
		fprintf(fp,"%15.5lf %15.5lf %6d\n",a+(b-a)/2.0,hist[i]*1.0/icount,hist[i]);
	}
	fclose(fp);
}

/*	Create histogram for the values of filename. The value
 *	of delta determines the gap.
 * */
void	makeHist(string filename,string histname,Data minvalue,Data maxvalue,Data delta)
{
	minvalue = floor(minvalue);
	maxvalue = ceil(maxvalue);
	extern int histflag;
	if(histflag==0) return;
	int nsize = (int)((maxvalue-minvalue)/delta)+1;
	vector<int> hist;
	hist.resize(nsize);
	for(int i=0;i<nsize;i++) hist[i]=0;
	FILE *fp=fopen(filename.c_str(),"r");
	int icount=0;
	while(1)
	{
		int frame;
		Data value;
		int k=fscanf(fp,"%d %f",&frame,&value);
		if(k<=0) break;
		icount++;
		for(int i=0;i<nsize;i++)
		{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(i==nsize-1 && fabs(value)==180.0) hist[i]++;
			else
			if(value>a && value<=b) hist[i]++;
		}
	}
	fclose(fp);
	fp=fopen(histname.c_str(),"w");
	for(int i=0;i<nsize;i++)
	{
			Data a=minvalue+i*delta;
			Data b=minvalue+(i+1)*delta;
			if(b>maxvalue)
		fprintf(fp,"%15.5lf %15.5lf %6d\n",a,hist[i]*1.0/icount,hist[i]);
		else
		fprintf(fp,"%15.5lf %15.5lf %6d\n",a+(b-a)/2.0,hist[i]*1.0/icount,hist[i]);
	}
	fclose(fp);
}


/*	Print numbers with a format.
 * */
string printNumber(int x)
{
	stringstream st;
	st<<x;
	return st.str();
	/*
	char s1[100];
	sprintf(s1,"%d",x);
	return s1;
	*/
}

string printNumber0(int x,int space)
{
	char s1[100];
	string format="%0"+printNumber(space)+"d";
	sprintf(s1,format.c_str(),x);
	return s1;
}

string printNumber(int x,int space)
{
	stringstream st;
	st<<setw(space)<<x;
	return st.str();
	/*
	char s1[100];
	string format="%"+printNumber(space)+"d";
	sprintf(s1,format.c_str(),x);
	return s1;
	*/
}

string printNumber(Data x)
{
	stringstream st;
	st<<x;
	return st.str();
	/*
	char s1[100];
	sprintf(s1,"%lf",x);
	*/
}

string printNumber(Data x,int space,int decimal)
{
	
	stringstream st;
	st<<fixed<<setw(space)<<setprecision(decimal)<<x;
	return st.str();
	/*	
	string format="%"+printNumber(space)+"."+printNumber(decimal)+"lf";
	char s1[100];
	sprintf(s1,format.c_str(),x);
	return s1;
	*/
}

string printFrame(int f)
{
	return printNumber(f,FRAME_WIDTH);
}

string printDistance(Data x)
{
	return printNumber(x,DISTANCE_WIDTH,DISTANCE_DECIMAL);
}

string printAngle(Data a)
{
	extern int zero_angle_flag;
	extern int degrees_angle_flag;
	if(zero_angle_flag) 
		if(a<0) a+=360;
	if(degrees_angle_flag)
	{
		a=deg2rad(a);
		return printNumber(a,ANGLE_WIDTH,4*ANGLE_DECIMAL);
	}
	return printNumber(a,ANGLE_WIDTH,ANGLE_DECIMAL);
}

string	printPercent(Data x)
{
	return printNumber(x,PERCENT_WIDTH,PERCENT_DECIMAL);
}

string	printString(string s,int tab)
{
	stringstream st;
	st<<setw(tab)<<s;
	return st.str();
	/*
	char s1[1024];
	string format="%"+printNumber(tab)+"s";
	sprintf(s1,format.c_str(),s.c_str());
	return s1;
	*/
}

void	Print(string s)
{
	printf("%s",s.c_str());
}

void	Print(FILE *fp,string s)
{
	fprintf(fp,"%s",s.c_str());
}

void	PrintLine(string s)
{
	printf("%s\n",s.c_str());
}

void	PrintLine(FILE *fp,string s)
{
	fprintf(fp,"%s\n",s.c_str());
}
void	getAngleVectorStatistics(DoubleVector &x,Data &minvalue,Data &maxvalue,Data &avgvalue,Data &variance)
{
	minvalue=1e+100;
	maxvalue=-1e+100;
	avgvalue=0.0;
	Data x1=0.0;
	Data x2=0.0;
	int    icount=0;
	Data  sinaverage=0.0;
	Data  cosaverage=0.0;
	for(int i=0;i<x.size();i++)
	{
		int frame;
		Data value;
		value = x[i];
		icount++;
		if(value<minvalue) minvalue=value;
		if(value>maxvalue) maxvalue=value;
		Data angle=value * M_PI/180.0;
		sinaverage+=sin(angle);
		cosaverage+=cos(angle);
	}
	avgvalue=atan2(sinaverage/icount,cosaverage/icount);
	double epsilon = sqrt(1.0 -sinaverage/icount * sinaverage/icount-
			cosaverage/icount * cosaverage/icount);
	variance = asin(epsilon) * (1.0+(2.0/sqrt(3.0)-1.0)*epsilon*epsilon*epsilon);
	variance = variance * 180.0/M_PI;
	avgvalue=avgvalue*180.0/M_PI;
}

void	PrintDistances(string filename,IntVector &frame,DoubleVector &distance)
{
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+printDistanceHeader("Dist"));
	for(int i=0;i<frame.size();i++)
		PrintLine(fp,printFrame(frame[i])+" "+printDistance(distance[i]));
	fclose(fp);
}

void	PrintAngles(string filename,IntVector &frame,DoubleVector &angle)
{
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+printAngleHeader("Tors"));
	for(int i=0;i<frame.size();i++)
		PrintLine(fp,printFrame(frame[i])+" "+printAngle(angle[i]));
	fclose(fp);
}

void	PrintDistAndAngle(string filename,IntVector &frame,DoubleVector &distance,DoubleVector &angle)
{
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+printDistanceHeader("Dist")+
		" "+printTorsionHeader("Angle"));
	for(int i=0;i<frame.size();i++)
		PrintLine(fp,printFrame(frame[i])+" "+printDistance(distance[i])+"  "+
			printAngle(angle[i]));
	fclose(fp);
}

void	PrintHist(string filename,vector<HistStruct> &st)
{
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Dist",DISTANCE_WIDTH-1)+" "+printPercentHeader("frames(%)")+
		" "+printFrameHeader("frames(#)"));
	for(int i=0;i<st.size();i++)
		PrintLine(fp,printDistance(st[i].value)+" "+printPercent(st[i].percent)+
			" "+printFrame(st[i].count));
	fclose(fp);
}

void	PrintAngleHist(string filename,vector<HistStruct> &st)
{
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Angle",ANGLE_WIDTH-1)+" "+printPercentHeader("frames(%)")+
		" "+printFrameHeader("frames(#)"));
	for(int i=0;i<st.size();i++)
		PrintLine(fp,printAngle(st[i].value)+" "+printPercent(st[i].percent)+
			" "+printFrame(st[i].count));
	fclose(fp);
}

string printFrameHeader()
{
	return printFrameHeader("Frame");
}

string printDistanceHeader()
{
	return printFrameHeader("Dist");
}

string printAngleHeader()
{
	return printAngleHeader("Angle");
}

string	printTorsionHeader()
{
	return	printTorsionHeader("Torsion");
}

string	printFrameHeader(string s)
{
	return printString(s,FRAME_WIDTH);
}

string printDistanceHeader(string s)
{
	return printString(s,DISTANCE_WIDTH);
}

string printAngleHeader(string s)
{
	return printString(s,ANGLE_WIDTH);
}

string printTorsionHeader(string s)
{
	return printString(s,ANGLE_WIDTH);
}

string printPercentHeader(string s)
{
	return printString(s,PERCENT_WIDTH);
}

void	PrintAcor(string filename,IntVector &F,IntVector Acor,int k,int w)

{
	return ; // It will be implemented later 
	IntVector2 Acor2;
	int parts=F.size()/k;
	if(F.size() % k!=0) parts++;
	Acor2.resize(parts);
	int icount=1;
	for(int i=0;i<parts;i++)
	{
		Acor2[i].resize(k);
		for(int j=0;j<k;j++)
		{
			Acor2[i][j]=Acor[icount++];
		}
	}
	FILE *fp=fopen(filename.c_str(),"w");
        if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Lag",FRAME_WIDTH-1)+" "+printString("cor",10));
	int ii=0;
	Data ff=1.0;
	PrintLine(fp,printFrame(ii)+" "+printNumber(ff,10,2));
	Data hmean=0.0;
	for(int i=1;i<Acor.size();i++) hmean+=1.0*Acor[i];
	hmean/=(Acor.size()-1);
	for(int i=0;i<k;i+=w)
	{
		Data h=0.0;
		for(int j=0;j<parts;j++)
			h+=1.0*Acor2[j][i];
		h/=parts;
		h/=hmean;
		int ii=i+1;
		PrintLine(fp,printFrame(ii)+" "+printNumber(h,10,2));	
	}	
	fclose(fp);
/*
	Data c0=0.0;
	Data ymean=0.0;
	for(int i=0;i<Acor.size();i++) ymean+=Acor[i];
	ymean/=Acor.size();
	for(int i=0;i<Acor.size();i++) c0+=pow(Acor[i]-ymean,2.0);
	c0/=Acor.size();
	if(fabs(c0)<1e-8) return;

	int dd=Acor[0];
	FILE *fp=fopen(filename.c_str(),"w");
        if(!fp) psf_error("Can not open "+filename+" for writting ");
	PrintLine(fp,"#"+printString("Lag",FRAME_WIDTH-1)+" "+printString("cor",10));
	int xx=0;Data yy=0.0;
	PrintLine(fp,printFrame(xx)+" "+printNumber(yy,10,2));
	for(int h=1;h<=Acor.size()-1;h+=astep)
	{
		Data ch=0.0;
		for(int t=1;t<=Acor.size()-h;t++)
		{
			ch+=(Acor[t-1]-ymean)*(Acor[t+h-1]-ymean);
		}
		ch/=Acor.size();
		int xx=h;
		Data rh=ch/c0;
		PrintLine(fp,printFrame(xx)+" "+printNumber(rh,10,2));
	}	
	fclose(fp);
	Acor[0]=dd;
	*/
}

void	getVectorStatistics(int icount,IntVector &start,IntVector &last,IntVector &F,DoubleVector &x,
		DoubleVector &min,DoubleVector &max,DoubleVector &avg,DoubleVector &std)
{
	min.resize(icount+1);
	max.resize(icount+1);
	avg.resize(icount+1);
	std.resize(icount+1);
	start.resize(icount);
	last.resize(icount);
	if(x.size() % icount!=0)
	{
		min.resize(icount+2);
		max.resize(icount+2);
		avg.resize(icount+2);
		std.resize(icount+2);
		start.resize(icount+1);
		last.resize(icount+1);
	}
	DoubleVector xx;
	for(int i=0;i<icount;i++)
	{
		xx.resize(x.size()/icount);
		for(int j=0;j<x.size()/icount;j++) 
		{
			xx[j]=x[i*x.size()/icount+j];
		}
		getVectorStatistics(xx,min[i],max[i],avg[i],std[i]);
		start[i]=F[i*x.size()/icount+0];
		last[i]=F[i*x.size()/icount+x.size()/icount-1];
	}
	if(x.size() % icount!=0)
	{
		int remain = x.size()-x.size()/icount*icount;
		xx.resize(remain);
		int k=0;
		for(int j=x.size()-1;j>=x.size()-1-remain;j--)
			xx[k++]=x[j];
		start[icount]=F[x.size()-1-remain];
		last[icount]=F[F.size()-1];
		getVectorStatistics(xx,min[icount],max[icount],avg[icount],std[icount]);
		getVectorStatistics(x,min[icount+1],max[icount+1],avg[icount+1],std[icount+1]);
	}
	else
	getVectorStatistics(x,min[icount],max[icount],avg[icount],std[icount]);
}


void	getAngleVectorStatistics(int icount,IntVector &start,IntVector &last,IntVector &F,
	DoubleVector &x,DoubleVector &min,DoubleVector &max,DoubleVector &avg,DoubleVector &std)
{
	min.resize(icount+1);
	max.resize(icount+1);
	avg.resize(icount+1);
	std.resize(icount+1);
	start.resize(icount);
	last.resize(icount);
	if(x.size() % icount!=0)
	{
		min.resize(icount+2);
		max.resize(icount+2);
		avg.resize(icount+2);
		std.resize(icount+2);
		start.resize(icount+1);
		last.resize(icount+1);
	}
	DoubleVector xx;
	for(int i=0;i<icount;i++)
	{
		xx.resize(x.size()/icount);
		for(int j=0;j<x.size()/icount;j++) 
		{
			xx[j]=x[i*x.size()/icount+j];
		}
		getAngleVectorStatistics(xx,min[i],max[i],avg[i],std[i]);
		start[i]=F[i*x.size()/icount+0];
		last[i]=F[i*x.size()/icount+x.size()/icount-1];
	}
	if(x.size() % icount!=0)
	{
		int remain = x.size()-x.size()/icount*icount;
		xx.resize(remain);
		int k=0;
		for(int j=x.size()-1;j>=x.size()-1-remain;j--)
			xx[k++]=x[j];
		start[icount]=F[x.size()-1-remain];
		last[icount]=F[F.size()-1];
		getAngleVectorStatistics(xx,min[icount],max[icount],avg[icount],std[icount]);
		getAngleVectorStatistics(x,min[icount+1],max[icount+1],avg[icount+1],std[icount+1]);
	}
	else
	getAngleVectorStatistics(x,min[icount],max[icount],avg[icount],std[icount]);
}

void	PrintStatFile(string filename,string mode,int icount,string header,IntVector &F,DoubleVector &x)
{
	if(icount==0) icount=1;
	IntVector start;
	IntVector last;
	FILE *fp=fopen(filename.c_str(),mode.c_str());
	if(!fp) 
		psf_error("CAN NOT OPEN "+filename+" FOR WRITTING");
	PrintLine(fp,header);
	int d=icount;
	if(F.size() % icount!=0) d++;
	PrintLine(fp,printString("#Begin block",-15)+" "+printFrame(d));
	DoubleVector min,max,avg,std;
	getVectorStatistics(icount,start,last,F,x,min,max,avg,std);
	PrintLine(fp,"#"+printString("Block",FRAME_WIDTH-1)+" "+
		printString("Frame1",FRAME_WIDTH)+" "+printString("Frame2",FRAME_WIDTH)+" "+
		printString("Min",10)+" "+printString("Max",10)+" "+
		printString("Avg",10)+" "+printString("Std",10));
	for(int i=0;i<icount;i++)
	{
		int low=start[i];
		int upper=last[i];
		int ii=i+1;
		PrintLine(fp,printFrame(ii)+" "+printFrame(low)+" "+printFrame(upper)+" "+
			" "+printNumber(min[i],10,4)+" "+
			printNumber(max[i],10,4)+" "+printNumber(avg[i],10,4)+" "+printNumber(std[i],10,4));	
	}	
	if(F.size() % icount!=0) 
	{
		int low=start[icount];
		int upper=F[F.size()-1];
		int ii=icount+1;
		PrintLine(fp,printFrame(ii)+" "+printFrame(low)+" "+printFrame(upper)+" "+
			" "+printNumber(min[icount+1],10,4)+" "+
			printNumber(max[icount+1],10,4)+" "+printNumber(avg[icount+1],10,4)+" "+printNumber(std[icount+1],10,4));	
	}
	int low=F[0];
	int upper=F[F.size()-1];
	PrintLine(fp,printString("#End block",-15)+" "+printFrame(d));
	PrintLine(fp,"\n#Complete trajectory");
	PrintLine(fp,printFrameHeader(" ")+" "+printFrame(low)+" "+printFrame(upper)+"  "+
		printNumber(min[min.size()-1],10,4)+" "+printNumber(max[max.size()-1],10,4)+" "+
		printNumber(avg[avg.size()-1],10,4)+" "+printNumber(std[std.size()-1],10,4));
	PrintLine(fp,"\n");
	fclose(fp);
}

void	PrintAngleStatFile(string filename,string mode,int icount,string header,IntVector &F,DoubleVector &x)
{
	if(icount==0) icount=1;
	IntVector start;
	IntVector last;
	FILE *fp=fopen(filename.c_str(),mode.c_str());
	if(!fp) 
		psf_error("CAN NOT OPEN "+filename+" FOR WRITTING");
	PrintLine(fp,header);
	int d=icount;
	if(F.size() % icount!=0) d++;
	PrintLine(fp,printString("#Begin block",-15)+" "+printFrame(d));
	DoubleVector min,max,avg,std;
	getAngleVectorStatistics(icount,start,last,F,x,min,max,avg,std);
	PrintLine(fp,"#"+printString("Block",FRAME_WIDTH-1)+" "+
		printString("Frame1",FRAME_WIDTH)+" "+printString("Frame2",FRAME_WIDTH)+" "+
		printAngleHeader("Min")+" "+printAngleHeader("Max")+" "+
		printAngleHeader("Avg")+" "+printAngleHeader("Std"));
	for(int i=0;i<icount;i++)
	{
		int low=start[i];
		int upper=last[i];
		int ii=i+1;
		PrintLine(fp,printFrame(ii)+" "+printFrame(low)+" "+printFrame(upper)+" "+
			" "+printAngle(min[i])+" "+printAngle(max[i])+
			" "+printAngle(avg[i])+" "+printAngle(std[i]));	
	}	
	if(F.size() % icount!=0) 
	{
		int low=start[icount];
		int upper=F[F.size()-1];
		int ii=icount+1;
		PrintLine(fp,printFrame(ii)+" "+printFrame(low)+" "+printFrame(upper)+" "+
			" "+printAngle(min[icount+1])+" "+printAngle(max[icount+1])+
			" "+printAngle(avg[icount+1])+" "+printAngle(std[icount+1]));	
	}
	int low=F[0];
	int upper=F[F.size()-1];
	PrintLine(fp,printString("#End block",-15)+" "+printFrame(d));
	PrintLine(fp,"\n#Complete trajectory");
	PrintLine(fp,printFrameHeader(" ")+" "+printFrame(low)+" "+printFrame(upper)+"  "+
		printAngle(min[min.size()-1])+" "+printAngle(max[max.size()-1])+" "+
		printAngle(avg[avg.size()-1])+" "+printAngle(std[std.size()-1]));
	PrintLine(fp,"\n");
	fclose(fp);
}

void	SplitTable(int icount,IntVector &start,IntVector &last,IntVector &F,IntVector &x,IntVector &res)
{
	if(icount==0) icount=1;
	start.resize(icount);
	last.resize(icount);
	res.resize(icount);
	if(F.size() % icount !=0) 
	{
		start.resize(icount+1);
		last.resize(icount+1);
		res.resize(icount+1);
	}
	IntVector xx;
	for(int i=0;i<icount;i++)
	{
		xx.resize(x.size()/icount);
		res[i]=0;
		for(int j=0;j<x.size()/icount;j++) 
		{
			xx[j]=x[i*x.size()/icount+j];
			res[i]+=xx[j];
		}
		start[i]=F[i*x.size()/icount+0];
		last[i]=F[i*x.size()/icount+x.size()/icount-1];
	}
	if(x.size() % icount!=0)
	{
		int remain = x.size()-x.size()/icount*icount;
		xx.resize(remain);
		int k=0;
		res[icount]=0;
		for(int j=x.size()-1;j>=x.size()-1-remain;j--)
		{
			xx[k++]=x[j];
			res[icount]+=x[j];
		}
		start[icount]=F[x.size()-1-remain];
		last[icount]=F[x.size()-1];
	}
}

void	PrintIntStatFile(string statname,string mode,int icount,string startheader,
	string header,IntVector &F,IntVector2 &x)
{
	IntVector start;
	IntVector last;
	IntVector2 res;
	if(icount==0) icount=1;
	FILE *fp=fopen(statname.c_str(),mode.c_str());
	if(!fp) 
		psf_error("CAN NOT OPEN "+statname+" FOR WRITTING");
	int dd=icount;
	if(F.size() % icount !=0) dd++;
	res.resize(x.size());
	IntVector xx;
	for(int i=0;i<x.size();i++)
	{
		SplitTable(icount,start,last,F,x[i],xx);
		res[i]=xx;
	}
	PrintLine(fp,startheader);
	PrintLine(fp,printString("#Begin block",-15)+" "+printFrame(dd));
	PrintLine(fp,header);
	IntVector sum;
	sum.resize(x.size());
	for(int i=0;i<x.size();i++) sum[i]=0;
	for(int i=0;i<start.size();i++)
	{
		int low = start[i];
		int ii=i+1;
		int upper = last[i];
		Print(fp,printFrame(ii)+" "+printFrame(low)+" "+printFrame(upper)+" ");
		for(int j=0;j<x.size();j++)
		{
			Print(fp,printFrame(res[j][i])+" ");
			sum[j]+=res[j][i];
		}
		PrintLine(fp,"");
	}
	PrintLine(fp,printString("#End block",-15)+" "+printFrame(dd));
	PrintLine(fp,"\n#Complete trajectory");
	Print(fp,printFrameHeader(" ")+" "+printFrame(F[0])+" "+printFrame(F[F.size()-1])+" ");
	for(int i=0;i<x.size();i++) 
		Print(fp,printFrame(sum[i])+" ");
	PrintLine(fp,"");
	fclose(fp);
}

Data deg2rad(Data x)
{
	return x* M_PI/180.0;
}

Data rad2deg(Data x)
{
	return x*180.0/M_PI;
}

void	printVersion()
{
	int	major_version = 1;
	int	minor_version = 0;
	int	sub_version   = 0;
	char str[1024];
	sprintf(str,"%s %s",__DATE__,__TIME__);
	string  datetime=str;
	PrintLine("Euclidean computational biology - version "+
		printNumber(major_version)+"."+printNumber(minor_version)+"."+
		printNumber(sub_version));
	PrintLine("Please visit http://stavrakoudis.econ.uoi.gr/eucb");
	PrintLine("\tfor release news and updates");
	PrintLine("Ioannis G. Tsoulos - itsoulos@gmail.com");
	PrintLine("Athanassios Stavrakoudis - astavrak@cc.uoi.gr");

}
