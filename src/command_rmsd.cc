# include <command_rmsd.h>
# include <global.h>
# include <pdb.h>

CommandRmsd::CommandRmsd()
	:Command("RMSD")
{
}


void	CommandRmsd::Run()
{
	printf(" %d %d \n",dcd->nframes(),dcd->getnatoms());
	Team *rmsd_team = team;
	if(rmsd_team == NULL)  Error(SeqError("rmsd"));
	if(ref!=0)
	{
		vector<DoubleVector> D;
		if(ref<=0 || ref>dcd->getnatoms())
			psf_error("reference frame out of range");
		dcd->fetchAtoms(ref,D);
		pdbpos.resize(D.size());	
		for(int i=0;i<D.size();i++)
		{
			pdbpos[i].x=D[i][0];
			pdbpos[i].y=D[i][1];
			pdbpos[i].z=D[i][2];	
		}
	}
	else
	if(pdbpos.size()==0)
		psf_error("You must provide either pdb or reference frame.");
	IntVector AtomPos;
	string chain;
	vector<string> chainlist;
	IntVector ChainPos;
	
	extern vector<atom> table;
	for(int i=0;i<table.size();i++)
	{
		if(isHydro(table[i]))  continue;
		if(rmsd_team->find(table[i],chain))
		{
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
	}
	if(chainlist.size()==0) return;
	IntVector F;
	vector<DoubleVector> D;
	IntVector AtomPos2;
	int chain_count = chainlist.size();
	D.resize(chain_count);
	for(int i=0;i<chain_count;i++)
	{
		AtomPos2.resize(0);
		for(int j=0;j<AtomPos.size();j++)
			if(table[AtomPos[j]].chain_id==chainlist[i]) AtomPos2.push_back(AtomPos[j]);
		DoubleVector DD;
		dcd->getRmsd(AtomPos2,pdbpos,F,DD);
		D[i]=DD;
	}
	DoubleVector totalD;
	if(chainlist.size()>1)	
	dcd->getRmsd(AtomPos,pdbpos,F,totalD);
	FILE *fp=fopen("rmsd.dat","w");
	if(!fp) Error(WriteError("rmsd.dat"));
	Print(fp,"#");
	Print(fp,printString("frame",7)+" ");
	for(int i=0;i<chainlist.size();i++)
		Print(fp,printString(chainlist[i],8)+" ");
	if(chainlist.size()>1)
	PrintLine(fp,printDistanceHeader("Total"));
	else
	PrintLine(fp,"");
	for(int i=0;i<F.size();i++)
	{
		Print(fp,printFrame(F[i])+" ");
		for(int j=0;j<chainlist.size();j++)
		{
			Print(fp,printDistance(D[j][i])+" ");
		}
		if(chainlist.size()>1)
		PrintLine(fp,printDistance(totalD[i]));
		else
		PrintLine(fp,"");
	}
	fclose(fp);
	if(smooth_flag)
	{
		IntVector sf;
		vector<DoubleVector> sd;
		DoubleVector smooth_totalD;
		sd.resize(chainlist.size());
		for(int i=0;i<chainlist.size();i++) makeSmooth(F,D[i],sf,sd[i],smooth_start,smooth_step);
		makeSmooth(F,totalD,sf,smooth_totalD,smooth_start,smooth_step);
		fp=fopen("rmsd.sda","w");
		Print(fp,"#");
		Print(fp,printString("frame",7)+" ");
		for(int i=0;i<chainlist.size();i++)
			Print(fp,printString(chainlist[i],8)+" ");
		PrintLine(fp,printDistanceHeader("Total"));
		for(int i=0;i<sf.size();i++)
		{
			Print(fp,printFrame(sf[i])+" ");
			for(int j=0;j<chainlist.size();j++)
			{
				Print(fp,printDistance(sd[j][i])+" ");
			}
			if(chainlist.size()>1)
			PrintLine(fp,printDistance(smooth_totalD[i]));
			else PrintLine(fp,"");
		}
		fclose(fp);
	}
	fp=fopen("rmsd.stat","w");
	fclose(fp);
	for(int i=0;i<chainlist.size();i++)
	{
		PrintStatFile("rmsd.stat","a",F.size()/statcount,"#Rmsd statistics for chain "+chainlist[i],F,D[i]);
	}
	if(chainlist.size()>1)
	PrintStatFile("rmsd.stat","a",F.size()/statcount,"#Total rmsd statistics",F,totalD);
	fp=fopen("rmsd.hist","w");
	typedef vector<HistStruct> histd;
	vector<histd> st;
	st.resize(chainlist.size());
	Data dmin,dmax,davg,dstd;
	for(int i=0;i<chainlist.size();i++)
	{
		getVectorStatistics(D[i],dmin,dmax,davg,dstd);
		makeHist(D[i],st[i],dmin,dmax,bindist);
	}
	
	for(int j=0;j<st[0].size();j++)
	{
		for(int i=0;i<chainlist.size();i++)
		{
		if(i==0) Print(fp,printDistance(st[i][j].value)+" ");
		Print(fp,printPercent(st[i][j].percent)+" "+printFrame(st[i][j].count)+" ");
		}
		PrintLine(fp,"");
	}
	fclose(fp);
}
