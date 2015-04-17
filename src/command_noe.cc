# include <global.h>
# include <math.h>
# include <command_noe.h>

CommandNoe::CommandNoe()
	:Command("NOE")
{
	noe_atom.resize(0);
	average=6.0;
}

void	CommandNoe::setAverage(Data d)
{	
	if(d>0) average = d;
}

void	CommandNoe::setAtoms(vector<string> &s)
{
	for(int i=0;i<s.size();i++) noe_atom.push_back(s[i]);
}

/*	Implement the noe facility.
 * */
void CommandNoe::Run()
{
	Team *noe_team = team;
	vector<int> atom_pos;
	string chain;
	for(int i=0;i<table.size();i++)
	{
		if(noe_team==NULL || noe_team->find(table[i],chain))
		{
			int found=0;
			for(int j=0;j<noe_atom.size();j++)
			{
				if(table[i].atom_name == noe_atom[j])
					{
						found = 1;
						break;
					}
			}
			if(found) atom_pos.push_back(i);
		}
	}
	vector<NoeStruct> noe;
	for(int i=0;i<atom_pos.size();i++)
	{
		for(int j=0;j<i;j++)
		{
			if(atom_pos[i]==atom_pos[j]) continue;
			if(table[atom_pos[i]].res_id==table[atom_pos[j]].res_id) continue;
			int s=noe.size();
			noe.resize(s+1);
			noe[s].atom1=atom_pos[i]+1;
			noe[s].atom2=atom_pos[j]+1;
		}
	}
	IntVector f;
	DoubleVector d;
	dcd->getNoe(f,noe);
	Log("#Noe critical distance = "+printDistance(noe_critical_distance));
	Log("#Noe critical percent  = "+printPercent(noe_critical_percent));
	for(int i=0;i<noe.size();i++)
	{
		d=noe[i].D;
		int count=0;
		Data minDistance=1e+100;
		Data maxDistance=-1e+100;
		Data avgDistance;
		Data stdDistance;
		Data xx1=0.0;
		Data xx2=0.0;
		int count28=0;
		int count35=0;
		int count5=0;
		Data lmargin = 1.6;
		Data rmargin = 5.5;
		Data delta = bindist;
		int countHisto = 1+(int)((rmargin-lmargin)/delta);
		vector<int>  histo;
		histo.resize(countHisto);
		for(int k=0;k<countHisto;k++) histo[k]=0;
		IntVector2 Acor;
		Acor.resize(4);
		Acor[0].resize(f.size());
		Acor[1].resize(f.size());
		Acor[2].resize(f.size());
		Acor[3].resize(f.size());
		for(int k1=0;k1<4;k1++)	
			for(int k2=0;k2<f.size();k2++)
				Acor[k1][k2]=0;
		for(int k=0;k<f.size();k++)
		{
			if(d[k]<=noe_critical_distance)
			{
				 count++;
				Acor[0][k]=1;
			}
			else Acor[0][k]=0;
			if(d[k]<2.8) 
			{
				count28++;
				Acor[1][k]=1;
			}
			else	Acor[1][k]=0;
		
			if(d[k]>=2.8 && d[k]<=3.5) 
			{
				count35++;
				Acor[2][k]=1;
			}
			else	Acor[2][k]=0;
			if(d[k]>3.5 && d[k]<=5.5) 
			{
				count5++;
				Acor[3][k]=1;
			}
			else	Acor[3][k]=0;
			for(int l=0;l<countHisto;l++)
			 if(d[k]>=lmargin+l*delta && d[k]<=lmargin+(l+1)*delta)
				histo[l]++;
		}
		if(count*1.0/f.size()<noe_critical_percent) continue;
		getVectorStatistics(d,minDistance,maxDistance,avgDistance,stdDistance);
		Data paverage=getPAverage(d,average);
		atom a=table[noe[i].atom1-1];
		atom b=table[noe[i].atom2-1];
		string filename="noe_"+printAtom1(a)+"_"+printAtom1(b)+".dat";
		string acorname="noe_"+printAtom1(a)+"_"+printAtom1(b)+".cor";
                PrintAcor(acorname,f,Acor[0],kstep,astep);

		Log(printAtom3(a)+"  "+printAtom3(b)+"   "+printPercent(count28*1.0/f.size())+"  "+
				printPercent(count35*1.0/f.size())+"  "+printPercent(count5*1.0/f.size())+" "+
				printDistance(paverage)); 	
		PrintDistances(filename,f,d);
		if(smooth_flag)
		{
			string smoothname="noe_"+printAtom1(a)+"_"+printAtom1(b)+".sda";
			IntVector sf;
			DoubleVector sd;
			makeSmooth(f,d,sf,sd,smooth_start,smooth_step);
			PrintDistances(smoothname,sf,sd);
		}
		string statname="noe_"+printAtom1(a)+"_"+printAtom1(b)+".stat";
		string header1="#Distance statistics";
		PrintStatFile(statname,"w",f.size()/statcount,header1,f,d);
		string header2="#Existence statistics";
		string header3="#"+printString("Block",FRAME_WIDTH-1)+" "+
			printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
			printFrameHeader("<5.5")+" "+printFrameHeader("<2.8")+" "+
			printFrameHeader("2.8-3.5")+" "+printFrameHeader("3.5-5.5");
		PrintIntStatFile(statname,"a",f.size()/statcount,header2,header3,f,Acor);

		string histname="noe_"+printAtom1(a)+"_"+printAtom1(b)+".hist";
		FILE *fp=fopen(histname.c_str(),"w");
		if(!fp) Error(WriteError(histname));
		for(int k=0;k<countHisto;k++)
		{
			PrintLine(fp,printDistance(lmargin+k*delta+delta/2.0)+" "+
				printPercent(histo[k]*1.0/f.size()));
		}
		fclose(fp);
	}
}
