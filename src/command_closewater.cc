# include <command_closewater.h>
# include <global.h>
# include <math.h>
CommandCloseWater::CommandCloseWater():
	Command("HPC")
{
}

typedef struct
{
	int 		atom;
	IntVector 	list;
}wpair;


typedef struct
{
	int		id;
	IntVector	list;
	DoubleVector2	dist;
}cluster;

static	int findCluster(IntVector2 &C,int id)
{
	for(int i=0;i<C.size();i++)
	{
		int ipos=isin(C[i],id);
		if(ipos!=-1) return i;
	}
	return -1;
}

void	CommandCloseWater::Run()
{
	IntVector2 SideList;
	IntVector  WaterList;
	string chain;
	for(int i=0;i<table.size();i++) 
	{
		if(isWater(table[i]) && table[i].atom_name=="OH2") WaterList.push_back(table[i].atom_id);
		if(!team || team->find(table[i],chain))
		{
			if(isSidechain(table[i]) )
			{
				int ipos=-1;
				for(int j=0;j<SideList.size();j++)
				{
					for(int k=0;k<SideList[j].size();k++)
					{
						if(table[SideList[j][k]-1].res_id==table[i].res_id
						 && table[SideList[j][k]-1].chain_id==table[i].chain_id) 
						{
							ipos=j;
							break;
						}
					}
				}
				if(ipos!=-1)
				{
					SideList[ipos].push_back(table[i].atom_id);
				}
				else
				{
					int s=SideList.size();
					SideList.resize(s+1);
					SideList[s].push_back(table[i].atom_id);
				}
			}
		}
	}
	vector<NoeStruct> noe;
	for(int i=0;i<SideList.size();i++)
	{
		for(int j=0;j<SideList[i].size();j++)
		{
			for(int k=0;k<WaterList.size();k++)
			{
				NoeStruct n;
				n.atom1=SideList[i][j];
				n.atom2=WaterList[k];
				noe.push_back(n);
			}
		}
	}
	IntVector F;
	dcd->setFrameStep(smart_skip);
	dcd->getNoe(F,noe);
	DoubleVector D;
	int icount=0;
	IntVector2 WaterList2;
	WaterList2.resize(SideList.size());
	IntVector	SideFlag;
	SideFlag.resize(SideList.size());
	for(int i=0;i<SideList.size();i++) SideFlag[i]=0;
	for(int i=0;i<SideList.size();i++)
	{
		D.resize(F.size());
		for(int k=0;k<F.size();k++) D[k]=1e+100;
		for(int j=0;j<SideList[i].size();j++)
		{
			for(int k=0;k<WaterList.size();k++)
			{
				for(int m=0;m<F.size();m++) 
				{
					if(noe[icount].D[m]<D[m]) D[m]=noe[icount].D[m];
					if(noe[icount].D[m]<smart_distance)
					{
						int ipos=isin(WaterList2[i],WaterList[k]);
						if(ipos==-1) 
							WaterList2[i].push_back(WaterList[k]);
				
						int res_id=table[SideList[i][0]-1].res_id;
						Data p=getBelowPercent(D,smart_distance);
					}
				}
				icount++;
			}
		}
		SideFlag[i]=1;
	}
	noe.resize(0);
	IntVector2 Near;
	Near.resize(SideList.size());
	DoubleVector2 Dist;
	dcd->getNear(F,SideList,Dist);
	icount=0;
	for(int i=0;i<SideList.size();i++) SideFlag[i]=0;
	for(int i=0;i<SideList.size();i++)
	{
		for(int j=0;j<i;j++)
		{
			Data p=getBelowPercent(Dist[i*SideList.size()+j],smart_distance);
			if(p>=critical_percent )
			{
				Near[i].push_back(j);
				SideFlag[i]=1;
				SideFlag[j]=1;
			}
		}
	}
	dcd->setFrameStep(step);
	for(int i=0;i<SideList.size();i++)
	{
		if(SideFlag[i]==0) continue;
		for(int j=0;j<SideList[i].size();j++)
		{
			for(int k=0;k<WaterList2[i].size();k++)
			{
				NoeStruct n;
				n.atom1=SideList[i][j];
				n.atom2=WaterList2[i][k];
				noe.push_back(n);
				icount++;
			}
		}
	}
	icount=0;
	dcd->getNoe(F,noe);
	for(int i=0;i<SideList.size();i++)
	{
		if(SideFlag[i]==0) continue;
		DoubleVector D;
		D.resize(F.size());
		for(int m=0;m<F.size();m++) D[m]=1e+100;
		for(int j=0;j<SideList[i].size();j++)
		{
			for(int k=0;k<WaterList2[i].size();k++)
			{
				for(int m=0;m<F.size();m++)
					if(noe[icount].D[m]<D[m]) D[m]=noe[icount].D[m];
				icount++;
			}
			Data p=getBelowPercent(D,3.5);
			if(p>=critical_percent)
			{
				SideFlag[i]=0;
			}
		}
	}
	
	IntVector2 SideList2;
	for(int i=0;i<SideList.size();i++)
	{
		if(SideFlag[i]) SideList2.push_back(SideList[i]);
	}
	dcd->getNear(F,SideList2,Dist);
	icount=0;
	Near.resize(0);
	Near.resize(SideList2.size());
	SideFlag.resize(SideList2.size());
	for(int i=0;i<Near.size();i++) 
	{	Near[i].resize(0);
		SideFlag[i]=0;
	}
	for(int i=0;i<SideList2.size();i++)
	{		
		Near[i].resize(0);
		for(int j=0;j<i;j++)
		{
			Data p=getBelowPercent(Dist[i*SideList2.size()+j],critical_distance);
			if(p>=critical_percent)
			{
				Near[i].push_back(j);
				SideFlag[j]=1;
				SideFlag[i]=1;
			}
		}	
	}
	icount=0;
	int dcount=0;

	IntVector2 Cluster;
	for(int i=0;i<Near.size();i++)
	{
		if(SideFlag[i]==0 ||  Near[i].size()==0) continue;
		int ipos=findCluster(Cluster,i);
	
		if(ipos==-1) 
		{
			for(int j=0;j<Near[i].size();j++)
			{
				ipos=findCluster(Cluster,Near[i][j]);
				if(ipos!=-1) break;
			}
			if(ipos==-1)
			{
				int s=Cluster.size();
				Cluster.resize(s+1);
				Cluster[s].push_back(i);
				ipos=s;
			}
		}
		for(int j=0;j<Near[i].size();j++) 
		{ 
			if(findCluster(Cluster,Near[i][j])==-1)
			Cluster[ipos].push_back(Near[i][j]);
		}
	}
	for(int i=0;i<Cluster.size();i++)
	{
		int ii=i+1;
		string filename="hpc_"+printNumber(ii)+".dat";
		string statname="hpc_"+printNumber(ii)+".stat";
		FILE *fp=fopen(filename.c_str(),"w");
		Print(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" ");
		Log("Cluster "+printNumber(ii));
		for(int j=0;j<Cluster[i].size();j++)
		{
			string st1=printAtom5(table[SideList2[Cluster[i][j]][0]-1]);
			Log("\tResidue "+st1);
			for(int k=0;k<j;k++)
			{
				string st2=printAtom5(table[SideList2[Cluster[i][k]][0]-1]);
				Print(fp,""+st1+"-"+st2+""+" ");
			}
		}
		PrintLine(fp,"");
		for(int m=0;m<F.size();m++)
		{
			Print(fp,printFrame(F[m])+" ");
			for(int j=0;j<Cluster[i].size();j++)
			{
				for(int k=0;k<j;k++)
				{
					int id1=Cluster[i][j];
					int id2=Cluster[i][k];
					if(id1<id2) 
					{
					int t=id1;
					id1=id2;
					id2=t;
					}
					Data d=Dist[id1*SideList2.size()+id2][m];
					Print(fp,printDistance(d)+" ");
				}
			}
			PrintLine(fp,"");
			
		}
		fclose(fp);
		fp=fopen(statname.c_str(),"w");
		PrintLine(fp,"#Statistics");
		fclose(fp);
		for(int j=0;j<Cluster[i].size();j++)
		{
			for(int k=0;k<j;k++)
			{
				int id1=Cluster[i][j];
				int id2=Cluster[i][k];
				if(id1<id2) 
				{
					int t=id1;
					id1=id2;
					id2=t;
				}
				string header1="#Distance statistics  for "+
					printAtom5(table[SideList2[Cluster[i][j]][0]-1])+" "+
					printAtom5(table[SideList2[Cluster[i][k]][0]-1]);
				PrintStatFile(statname,"a",F.size()/statcount,header1,F,
					Dist[id1*SideList2.size()+id2]);
			
			}
		}
	}
}


CommandCloseWater::~CommandCloseWater()
{

}
