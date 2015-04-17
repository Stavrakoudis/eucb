# include <command_clusterrmsd.h>
# include <global.h>

CommandClusterRmsd::CommandClusterRmsd(double l)
	:Command("CLUSTERRMSD")
{
	limit = l;
}


void	CommandClusterRmsd::Run()
{
	Team *rmsd_team = team;
	if(rmsd_team == NULL)  Error(SeqError("rmsd"));
	IntVector AtomPos;
	string chain;
	vector<string> chainlist;
	IntVector ChainPos;
	
	extern vector<atom> table;
	for(int i=0;i<table.size();i++)
	{
		if(isHydro(table[i]))  continue;
		string c;
		if(!rmsd_team->find(table[i],c)) continue;

		AtomPos.push_back(i);
		int ifound=-1;
		for(int j=0;j<chainlist.size();j++)
		{
			if(chainlist[j]==table[i].chain_id) {ifound=j;break;}
		}
		if(ifound==-1)
		{
			ifound=chainlist.size();
			chainlist.push_back(table[i].chain_id);
		}
		Log(printAtom3(table[i]));
		ChainPos.push_back(ifound);
	}
	if(chainlist.size()==0) return;
	IntVector F;
	vector<DoubleVector> D;
	IntVector AtomPos2;
	int chain_count = chainlist.size();
	D.resize(chain_count);

	IntVector2 FrameFlag;
	IntVector  Pivot;

	for(int i=0;i<chain_count;i++)
	{
		AtomPos2.resize(0);
		for(int j=0;j<AtomPos.size();j++)
			if(table[AtomPos[j]].chain_id==chainlist[i]) AtomPos2.push_back(AtomPos[j]);
		DoubleVector2 DD;
		dcd->getRmsd(AtomPos2,F,DD);
		if(FrameFlag.size()==0)
		{
			Pivot.resize(F.size());
			FrameFlag.resize(F.size());
			for(int j=0;j<FrameFlag.size();j++) 
			{
				FrameFlag[j].resize(0);
				FrameFlag[j].push_back(j);
				Pivot[j]=j;
			}
		}
		for(int j=0;j<DD.size();j++)
		{
			for(int k=0;k<DD[j].size();k++)
			{
				if(k==j) continue;
				
				if(DD[j][k]<=limit && DD[j][k]>=1e-5)
				{
					int ipos=isin(FrameFlag[j],k);
					if(ipos==-1) FrameFlag[j].push_back(k);
				}
			}
		}
	}
	
	if(chainlist.size()>1)	
	{
		DoubleVector2 DD;
		dcd->getRmsd(AtomPos,F,DD);
		if(FrameFlag.size()==0)
		{
			Pivot.resize(F.size());
			FrameFlag.resize(F.size());
			for(int j=0;j<FrameFlag.size();j++) 
			{
				FrameFlag[j].resize(0);
				FrameFlag[j].push_back(j);
				Pivot[j]=j;
			}
		}
		for(int j=0;j<DD.size();j++)
		{
			for(int k=0;k<DD[j].size();k++)
			{
				if(k==j) continue;
				
				if(DD[j][k]<=limit && DD[j][k]>=1e-5)
				{
					int ipos=isin(FrameFlag[j],k);
					if(ipos==-1) FrameFlag[j].push_back(k);
				}
			}
		}
	}

	for(int i=0;i<FrameFlag.size();i++)
	{
		for(int j=0;j<FrameFlag.size()-1;j++)
		{
			if(FrameFlag[j+1].size()>FrameFlag[j].size())
			{
				int t=Pivot[j];
				Pivot[j]=Pivot[j+1];
				Pivot[j+1]=t;
				IntVector temp;
				temp=FrameFlag[j];
				FrameFlag[j]=FrameFlag[j+1];
				FrameFlag[j+1]=temp;
			}
		}
	}
	
	
	for(int i=0;i<FrameFlag.size();i++)
	{
		for(int k=0;k<FrameFlag[i].size();k++)
		{
			for(int j=0;j<FrameFlag.size();j++)
			{
				if(i==j) continue;
				if(Pivot[j]==FrameFlag[i][k]) 
				{
					FrameFlag[j].resize(0);
				}
				else
				if(FrameFlag[j].size())
				{
					int ipos=isin(FrameFlag[j],FrameFlag[i][k]);
					if(ipos!=-1)
					{
						int s=FrameFlag[j].size();
						for(int l=ipos;l<s-1;l++)
							FrameFlag[j][l]=FrameFlag[j][l+1];
						FrameFlag[j].resize(s-1);
					}
				}
			}
		}
	}
	FILE *fp=fopen("clusterrmsd.dat","w");
	if(!fp) Error(WriteError("clusterrmsd.dat"));
	Print(fp,"#");
	Print(fp,printString("frame",FRAME_WIDTH-1)+" "+printFrameHeader("cluster"));
	PrintLine(fp,"");
	for(int i=0;i<FrameFlag.size();i++)
	{
		if(F[i]==0) F[i]=1;
		Print(fp,printFrame(F[i])+" ");
		int flag=0;
		for(int j=0;j<FrameFlag.size();j++)
		{
			for(int k=0;k<FrameFlag[j].size();k++)
			{
				if(FrameFlag[j][k]==i)
				{
					PrintLine(fp,printFrame(Pivot[j]+1));
					flag=1;
					break;
				}
			}
		}
		if(!flag)
			PrintLine(fp,printFrame(Pivot[i]+1));
	}
	fclose(fp);
	for(int i=0;i<FrameFlag.size();i++)
	{
		for(int j=0;j<FrameFlag.size()-1;j++)
		{
			if(FrameFlag[j+1].size()>FrameFlag[j].size())
			{
				int t=Pivot[j];
				Pivot[j]=Pivot[j+1];
				Pivot[j+1]=t;
				IntVector temp;
				temp=FrameFlag[j];
				FrameFlag[j]=FrameFlag[j+1];
				FrameFlag[j+1]=temp;
			}
		}
	}
	
	fp=fopen("clusterrmsd.stat","w");
	PrintLine(fp,"#"+printString("Cluster",FRAME_WIDTH-1)+" "+
		printFrameHeader("Center")+"  "+printFrameHeader("Count"));
	int icount=1;
	for(int i=0;i<FrameFlag.size();i++)
	{
		if(FrameFlag[i].size()==0) continue;	
		PrintLine(fp,printFrame(icount)+" "+printFrame(Pivot[i]+1)+" "+
			printFrame(FrameFlag[i].size()));
		icount++;
	}
	fclose(fp);
}
