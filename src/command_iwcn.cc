# include <command_iwcn.h>
# include <global.h>

CommandIwcn::CommandIwcn(double p,double d,int m)
	:Command("IWCN")
{
	percent=p;
	distance=d;
	maxn = m;
}


void	CommandIwcn::Run()
{
	Log("#Critical distance="+printDistance(distance));
	Log("#Critcal  percent ="+printPercent(percent));
	IntVector list;
	IntVector f;
	DoubleVector d;
	IntVector OrigWater;
	IntVector TeamAtoms;
	vector<CountFramesStruct> firstCheck;
	if(team) team->enumerateAtoms(table,TeamAtoms);
	for(int i=0;i<table.size();i++)
	{	
		atom a=table[i];
		if(!isWater(a))  continue;
		if(table[i].atom_name!="OH2") continue;
		list.push_back(i+1);
		if(!team) OrigWater.push_back(i+1);
		else
		for(int j=0;j<TeamAtoms.size();j++)
		{
			int ipos=TeamAtoms[j];
			if(table[ipos].atom_name[0]!='H')
			{
				int st=firstCheck.size();
				firstCheck.resize(st+1);
				firstCheck[st].atom1=ipos+1;
				firstCheck[st].atom2=i+1;
			}
		}
	}

	dcd->setFrameStep(smart_skip);
	dcd->countFrames(3 * distance,f,firstCheck);
	for(int i=0;i<firstCheck.size();i++)
	{
		int icount=firstCheck[i].count;
		//if(icount*1.0/f.size()>=percent)
		if(icount>=1)
		{
			if(isin(OrigWater,firstCheck[i].atom2)==-1) 
			{
				OrigWater.push_back(firstCheck[i].atom2);
			}
		}
	}

	
	vector<IwcnStruct> iwcn;
	iwcn.resize(OrigWater.size());
	for(int i=0;i<OrigWater.size();i++)
	{
		iwcn[i].water=OrigWater[i];
//		iwcn[i].list=list;
	}
	for(int i=0;i<OrigWater.size();i++)
	{
		iwcn[i].list=OrigWater;
	}
	dcd->setFrameStep(smart_skip);
	dcd->getIwcn(smart_distance,f,iwcn);

	OrigWater.resize(0);
	for(int i=0;i<iwcn.size();i++)
	{
		int icount=0;
		for(int j=0;j<f.size();j++)
			if(iwcn[i].count[j]<=maxn) icount++;
		if(icount*1.0/f.size()>=percent) OrigWater.push_back(iwcn[i].water);
	}
	dcd->setFrameStep(step);
	iwcn.resize(OrigWater.size());
	for(int i=0;i<OrigWater.size();i++)
	{
		iwcn[i].water=OrigWater[i];
	//	iwcn[i].list=list;
	}
	for(int i=0;i<OrigWater.size();i++)
	{
		iwcn[i].list=OrigWater;
	}
	dcd->getIwcn(distance,f,iwcn);
	int oxygen_pairs=0;
	for(int i=0;i<iwcn.size();i++)
	{
		int icount=0;
		IntVector Acor;
		Acor.resize(f.size());
		for(int j=0;j<f.size();j++)
		{
			if(iwcn[i].count[j]<=maxn) 
			{
				icount++;
				Acor[j]=1;
			}
			else Acor[j]=0;
		}
		if(icount*1.0/f.size()>=percent)
		{
			oxygen_pairs++;
			Log("Found iwcn match for water "+printAtom4(table[iwcn[i].water-1]));
			atom a=table[iwcn[i].water-1];
			string filename="iwcn_"+printAtom1(a)+".dat";
			string statname="iwcn_"+printAtom1(a)+".stat";
			string histname="iwcn_"+printAtom1(a)+".hist";
			string acorname="iwcn_"+printAtom1(a)+".cor";
                        PrintAcor(acorname,f,Acor,kstep,astep);
	
			DoubleVector D;
			D.resize(iwcn[i].count.size());
			FILE *fp=fopen(filename.c_str(),"w");
			if(!fp) Error(WriteError(filename));
			for(int k=0;k<f.size();k++)
			{
				Print(fp,printFrame(f[k])+" "+printNumber(iwcn[i].count[k])+" ");
				for(int i1=0;i1<iwcn[i].atomList[k].size();i1++)
				{
					int id=iwcn[i].atomList[k][i1];
					atom atomid=table[id-1];
					Print(fp,printAtom4(atomid)+" ");
				}
				PrintLine(fp,"");
				D[k]=iwcn[i].count[k];
			}
			fclose(fp);
			Data minDistance,maxDistance,avgDistance,stdDistance;
			string header1="#Iwcn statistics";
                	getVectorStatistics(D,minDistance,maxDistance,avgDistance,stdDistance);
			PrintStatFile(statname,"w",f.size()/statcount,header1,f,D);
			if(histflag)
			{
				int count_freq[16]={0};
				for(int k=0;k<f.size();k++)
				{
					count_freq[(int)D[k]]++;
				}
				fp=fopen(histname.c_str(),"w");
				if(!fp) Error(WriteError(histname));
				PrintLine(fp,"#"+printString("OH2",4)+" "+printFrameHeader("Count")+" "+printPercentHeader("(%)"));
				for(int k=0;k<16;k++)
				{
					PrintLine(fp,printNumber(k,5)+" "+printFrame(count_freq[k])+" "+printPercent(count_freq[k]*1.0/f.size()));
				}
				fclose(fp);
			}
		}
	}
	Log("#Found    "+printNumber(oxygen_pairs)+" oxygen pairs ");
}
