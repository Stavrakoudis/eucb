# include <counter.h>
# include <math.h>
Counter::Counter()
{
	count.resize(0);
	diff_res.resize(0);
}

void	Counter::add(int Frame,string Residue,int Count)
{
	int found=0;
	for(int i=0;i<diff_res.size();i++)
		if(diff_res[i]==Residue) {found=1;break;}
	if(Residue!="" && !found) diff_res.push_back(Residue);

	found=-1;
	for(int i=0;i<count.size();i++)
		if(count[i].frame==Frame) {found=i;break;}
	if(found==-1)
	{
		int ss=count.size();
		count.resize(ss+1);
		found=ss;
		count[ss].frame=Frame;
		count[ss].total=0;
		count[ss].residue.resize(0);
		count[ss].count.resize(0);
	}
	count[found].total+=Count;
	if(Count==0) return;

	int lfound=0;
	for(int i=0;i<count[found].residue.size();i++)
		if(count[found].residue[i]==Residue)
		{
			lfound=1;
			count[found].count[i]+=Count;
			break;
		}
	if(!lfound) 
	{
		count[found].residue.push_back(Residue);
		count[found].count.push_back(Count);
	}
}

void	Counter::print(string filename)
{
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp)
		psf_error("UNABLE TO OPEN "+filename+" FOR WRITTING ");
	Print(fp,printString("#Frame",-8)+" "+printString("TOTAL",7)+" ");
	for(int i=0;i<diff_res.size();i++)
		Print(fp,printString(diff_res[i],8)+" ");
	PrintLine(fp,"");
	for(int i=0;i<count.size();i++)
	{
		Print(fp,printNumber(count[i].frame,8)+" "+printNumber(count[i].total,7)+" ");
		for(int j=0;j<diff_res.size();j++)
		{
			int printCount=0;
			for(int k=0;k<count[i].residue.size();k++)
				if(count[i].residue[k]==diff_res[j])
				{
					printCount=count[i].count[k];
					break;
				}
			Print(fp,printNumber(printCount,8)+" ");
		}
		PrintLine(fp,"");
	}
	fclose(fp);
	
}

void	Counter::printStat(string filename)
{
	if(count.size()==0) return;
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp)
		psf_error("UNABLE TO OPEN "+filename+" FOR WRITTING ");
	vector<DoubleVector> D;
	D.resize(diff_res.size());
	DoubleVector Total;
	Total.resize(count.size());
	for(int i=0;i<D.size();i++) D[i].resize(count.size());
	for(int i=0;i<count.size();i++)
	{
		Total[i]=count[i].total;
		for(int j=0;j<diff_res.size();j++)
		{
			int printCount=0;
			for(int k=0;k<count[i].residue.size();k++)
				if(count[i].residue[k]==diff_res[j])
				{
					printCount=count[i].count[k];
					break;
				}
			D[j][i]=printCount;
		}
	}
	Data dmin,dmax,davg,dstd;
	getVectorStatistics(Total,dmin,dmax,davg,dstd);
	PrintLine(fp,printFrameHeader(" ")+" "+printFrameHeader("MIN")+" "+	
		printFrameHeader("MAX")+" "+printAngleHeader("AVG")+" "+printAngleHeader("STD"));
	PrintLine(fp,printString("TOTAL",7)+" "+printFrame((int)dmin)+" "+
		printFrame((int)dmax)+" "+printAngle(davg)+" "+
		printAngle(dstd));
	for(int i=0;i<D.size();i++)
	{
		getVectorStatistics(D[i],dmin,dmax,davg,dstd);
		PrintLine(fp,printString(diff_res[i],7)+" "+printFrame((int)dmin)+" "+
			printFrame((int)dmax)+" "+printNumber(davg,7,2)+" "+
			printNumber(dstd,7,2));
	}

	fclose(fp);
}

int	Counter::size()
{
	return count.size();
}

void	Counter::printHist(string filename)
{
	if(count.size()==0) return;
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp)
		psf_error("UNABLE TO OPEN "+filename+" FOR WRITTING ");
	DoubleVector Total;
	Total.resize(count.size());
	for(int i=0;i<count.size();i++)
	{
		Total[i]=count[i].total;
	}
	Data dmin,dmax,davg,dstd;
	getVectorStatistics(Total,dmin,dmax,davg,dstd);
	vector<HistStruct> st;
	makeHist(Total,st,dmin-1,dmax,1);

	PrintLine(fp,"#"+printString("Total",DISTANCE_WIDTH-1)+" "+printPercentHeader("frames(%)")+" "+
		printFrameHeader("frames(#)"));
	for(int i=0;i<st.size()-1;i++)
	{
		PrintLine(fp,printNumber((int)(ceil(st[i].value)),DISTANCE_WIDTH)+" "+printPercent(st[i].percent)+" "+
			printFrame(st[i].count));
	}
	fclose(fp);
}

void	Counter::printSmooth(string filename)
{
	IntVector F;
	IntVector sf;
	DoubleVector total;
	DoubleVector Stotal;
	
	vector<DoubleVector> D;
	vector<DoubleVector> SD;
	D.resize(diff_res.size());
	SD.resize(D.size());
	F.resize(count.size());
	total.resize(count.size());
	
	for(int i=0;i<D.size();i++) 
	{
		D[i].resize(F.size());
	}

	for(int i=0;i<count.size();i++)
	{
		F[i]=count[i].frame;
		total[i]=count[i].total;
		for(int j=0;j<diff_res.size();j++)
		{
			int printCount=0;
			for(int k=0;k<count[i].residue.size();k++)
				if(count[i].residue[k]==diff_res[j]) 
				{
					printCount=count[i].count[k];
					break;
				}
			D[j][i]=printCount;
		}
	}
	extern int smooth_start,smooth_step;
	makeSmooth(F,total,sf,Stotal,smooth_start,smooth_step);
	for(int i=0;i<D.size();i++)
	{	
		makeSmooth(F,D[i],sf,SD[i],smooth_start,smooth_step);
	}	
	
	FILE *fp=fopen(filename.c_str(),"w");
	if(!fp)
		psf_error("UNABLE TO OPEN "+filename+" FOR WRITTING ");
	Print(fp,printString("#Frame",-8)+" "+printString("TOTAL",7)+" ");
	for(int i=0;i<diff_res.size();i++)
		Print(fp,printString(diff_res[i],8)+" ");
	PrintLine(fp,"");
	for(int i=0;i<sf.size();i++)
	{
		double ss=0.0;
		for(int j=0;j<SD.size();j++) ss+=SD[j][i];
		if(ss>0) Stotal[i]=ss;
		Print(fp,printNumber(sf[i],8)+" "+printAngle(Stotal[i])+" ");
		for(int j=0;j<SD.size();j++)
			Print(fp,printAngle(SD[j][i])+" ");
		PrintLine(fp,"");
	}
	fclose(fp);
}
