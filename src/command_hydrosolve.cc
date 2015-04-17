# include <command_hydrosolve.h>
# include <global.h>
# include <dcd.h>

typedef struct 
{
	int  res_id1;
	string  chain_id1;
	int  res_id2;
	string  chain_id2;
}Cpair;

CommandHydroSolve::CommandHydroSolve()
	:Command("HYDROSOLVE")
{
	percent=0.05;
	distance=3.0;
	angle=145.0;
}

void	CommandHydroSolve::setPda(Data p,Data d,Data a)
{
	percent = p;
	distance = d;
	angle =a;
}

void	CommandHydroSolve::Run()
{
	if(team==NULL)  Error(SeqError("hydrosolve"));
	vector<int> waterlist;
	vector<int> atomlist;
	extern vector<atom> table;
	for(int  i=0;i<table.size();i++)
	{
		string chain;
		if(team->find(table[i],chain) && isHydroPhobic(table[i])) 
			atomlist.push_back(table[i].atom_id);
		if(isWater(table[i]) && table[i].atom_name[0]=='O') 
			waterlist.push_back(table[i].atom_id);
	}
	if(waterlist.size()==0)
	{
		Log("No water atoms found ");
		return;
	}
	vector<HydroSolveStruct> st;
	for(int i=0;i<atomlist.size();i++)
	{
		for(int j=i;j<atomlist.size();j++)
		{
			if(i==j) continue;
			if(table[atomlist[i]-1].res_id==table[atomlist[j]-1].res_id &&
			  table[atomlist[i]-1].chain_id==table[atomlist[j]-1].chain_id) continue;
			int s=st.size();
			st.resize(s+1);
			st[s].atom1=atomlist[i];
			st[s].atom2=atomlist[j];
		}
	}
	IntVector F;
	extern Dcd* dcd;
	IntVector waterlist2;
	dcd->setFrameStep(smart_skip);
	dcd->getHydroSolve(waterlist,F,st,smart_distance,angle);
	vector<HydroSolveStruct> st2;
	for(int i=0;i<st.size();i++)
	{
		int k=st[i].Frame.size();
		if(k) 
		{
			for(int j=0;j<st[i].Water.size();j++)
			{
				int found=0;
				for(int l=0;l<waterlist2.size();l++)
				{
					if(waterlist2[l]==st[i].Water[j]) 
					{
						found=1;
						break;
					}
				}
				if(!found) waterlist2.push_back(st[i].Water[j]);
			}
			st[i].Water.resize(0);
			st[i].Frame.resize(0);
			st[i].Distance1.resize(0);
			st[i].Distance2.resize(0);
			st2.push_back(st[i]);
		}
	}
	dcd->setFrameStep(step);
	dcd->getHydroSolve(waterlist2,F,st2,distance,angle);

	string filename,statname,histname;

	vector<Cpair> pair;
	
	for(int i=0;i<st2.size();i++)
	{
		int icount=0;
		atom a=table[st2[i].atom1-1];
		atom b=table[st2[i].atom2-1];
		int found=0;
		for(int j=0;j<pair.size();j++)
		{
			if(a.chain_id==pair[j].chain_id1 && a.res_id==pair[j].res_id1 &&
			   b.chain_id==pair[j].chain_id2 && b.res_id==pair[j].res_id2)
			{
				found=1;
				break;
			}
		}
		if(found) continue;
		int s=pair.size();
		pair.resize(s+1);
		pair[s].chain_id1=a.chain_id;
		pair[s].res_id1=a.res_id;
		pair[s].chain_id2=b.chain_id;
		pair[s].res_id2=b.res_id;
		
		IntVector   Index;
	
		for(int j=0;j<st2.size();j++)
		{
			atom c=table[st2[j].atom1-1];
			atom d=table[st2[j].atom2-1];
			if(a.res_id==c.res_id && a.chain_id==c.chain_id &&
			   b.res_id==d.res_id && b.chain_id==d.chain_id)
				Index.push_back(j);
		}

		IntVector myFrame;
		for(int j=0;j<Index.size();j++)
		{
			for(int l=0;l<st2[Index[j]].Frame.size();l++)
			{
				int found=0;
				for(int m=0;m<myFrame.size();m++)
				{
					if(myFrame[m]==st2[Index[j]].Frame[l])
					{
						found=1;
						break;
					}
				}
				myFrame.push_back(st2[Index[j]].Frame[l]);
			}
		}
		
		icount=myFrame.size();	
		Data p=icount*1.0/F.size();
		if(p<percent) continue;

		Log("Found water matching for atoms "+
			printAtom4(a)+" "+printAtom4(b)+" "+
			" with percent "+printPercent(p));
		
		filename="hsolve_"+a.chain_id+"_"+printNumber(a.res_id)+"_"+
			b.chain_id+"_"+printNumber(b.res_id)+".dat";
		FILE *fp=fopen(filename.c_str(),"w");
		if(!fp)  Error(WriteError(filename));
		for(int j=0;j<Index.size();j++)
		{
			for(int k=0;k<st2[Index[j]].Frame.size();k++)
			{
				Print(fp,printFrame(st2[Index[j]].Frame[k])+" ");
				atom c=table[st2[Index[j]].atom1-1];
				atom d=table[st2[Index[j]].atom2-1];
				Print(fp,printNumber(st2[Index[j]].atom1,5)+" "+
					 printNumber(st2[Index[j]].atom2,5)+" ");
				Print(fp,printNumber(st2[Index[j]].Water[k],5)+" ");
				PrintLine(fp,printDistance(st2[Index[j]].Distance1[k])+" "+
					 printDistance(st2[Index[j]].Distance2[k]));
			}
		}
		fclose(fp);
		
	}
}


CommandHydroSolve::~CommandHydroSolve()
{
}
