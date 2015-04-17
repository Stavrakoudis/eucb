# include <command_salt.h>
# include <pdb.h>
# include <global.h>
# include <tex.h>
# include <counter.h>

CommandSalt::CommandSalt()
	:Command("SALT")
{
	poscsalt.resize(6);
	negcsalt.resize(4);
	poscsalt[0].res_name="LYS";poscsalt[0].atom_name="NZ";
	poscsalt[1].res_name="ARG";poscsalt[1].atom_name="NE";
	poscsalt[2].res_name="ARG";poscsalt[2].atom_name="NH1";
	poscsalt[3].res_name="ARG";poscsalt[3].atom_name="NH2";
	poscsalt[4].res_name="HSP";poscsalt[4].atom_name="ND1";
	poscsalt[5].res_name="HSP";poscsalt[5].atom_name="NE2";

	negcsalt[0].res_name="ASP";negcsalt[0].atom_name="OD1";
	negcsalt[1].res_name="ASP";negcsalt[1].atom_name="OD2";
	negcsalt[2].res_name="GLU";negcsalt[2].atom_name="OE1";
	negcsalt[3].res_name="GLU";negcsalt[3].atom_name="OE2";
}

typedef struct
{
	int id1;
	int id2;
	int count1;
	int count2;
}SaltPair;

/*	Implement the complex salt facility.
 * */
void CommandSalt::Run()
{
	Counter myCounter;
	Tex myTex("salt2.tex");
	vector<string> s;
	s.resize(2);
	s[0]="Positive";
	s[1]="Negative";
	myTex.setHeader(s);
	myTex.setCaption("Salt bridges");


	Team *complexsalt_team = team;
	IntVector posVector,posVector2;
	IntVector negVector,negVector2;
	vector<PosTable2> PosTable;
	vector<PosTable2> PosTable3;

	for(int i=0;i<table.size();i++)
	{
		string chain;
		if(complexsalt_team && !complexsalt_team->find(table[i],chain)) continue;
		for(int j=0;j<poscsalt.size();j++)
			if(table[i].atom_type == "NH3" || 
			  (table[i].res_name == poscsalt[j].res_name && 
                           table[i].atom_name == poscsalt[j].atom_name))
			{
				int index=i;
				while(index>=0 && table[index].res_id == table[i].res_id) 
				{
					index--;
				}
				if(isin(posVector2,index+2)==-1)
				posVector2.push_back(index+2);
			}
		for(int j=0;j<negcsalt.size();j++)
			if(table[i].res_name == negcsalt[j].res_name && 
                           table[i].atom_name == negcsalt[j].atom_name)
			{
				int index=i;
				while(index>=0 && table[index].res_id == table[i].res_id) 
				{
					index--;
				}
				if(isin(negVector2,index+2)==-1)
				negVector2.push_back(index+2);
			}
	}
	IntVector f;
	dcd->setFrameStep(smart_skip);	
	dcd->getPosNeg(f,posVector2,negVector2,poscsalt,negcsalt,PosTable3);
	for(int i=0;i<PosTable3.size();i++)
	{
		int pos_id = i / negVector2.size();
		int neg_id = i % negVector2.size();
		int icount = 0;
		vector<SaltPair> ptable;
		
		for(int j=0;j<f.size();j++) 
		{
			int ipos=-1;
			
			for(int k=0;k<ptable.size();k++)
			{
				if(ptable[k].id1 == PosTable3[i][j].id1 
				&& ptable[k].id2 == PosTable3[i][j].id2) 
				{
					ipos=k;
					break;
				}
			}
			if(ipos==-1)
			{
				ipos=ptable.size();
				ptable.resize(ipos+1);
				ptable[ipos].id1=PosTable3[i][j].id1;
				ptable[ipos].id2=PosTable3[i][j].id2;
				ptable[ipos].count1=0;
				ptable[ipos].count2=0;
			}
			ptable[ipos].count1++;
			if(PosTable3[i][j].dist<=smart_distance) 
			{
				ptable[ipos].count2++;
				icount++;
			}
		}
		if(icount)
		{
			atom a,b;
			int atom_pos1,atom_pos2;
			int ipos=isin(posVector,posVector2[pos_id]);
			if(ipos==-1) posVector.push_back(posVector2[pos_id]);
			ipos=isin(negVector,negVector2[neg_id]);
			if(ipos==-1) negVector.push_back(negVector2[neg_id]);
		}
	}


	dcd->setFrameStep(step);
	dcd->getPosNeg(f,posVector,negVector,poscsalt,negcsalt,PosTable);

	Log("Critical distance = "+printDistance(critical_distance));
	Log("Critical percent  = "+printPercent(critical_percent));
	for(int i=0;i<PosTable.size();i++)
	{
		int pos_id = i / negVector.size();
		int neg_id = i % negVector.size();
		int icount = 0;
		vector<SaltPair> ptable;
		
		IntVector2 Acor;
		Acor.resize(1);
		Acor[0].resize(f.size());
		for(int j=0;j<f.size();j++) 
		{
			int ipos=-1;
			myCounter.add(f[j],"",0);
			
			for(int k=0;k<ptable.size();k++)
			{
				if(ptable[k].id1 == PosTable[i][j].id1 
				&& ptable[k].id2 == PosTable[i][j].id2) 
				{
					ipos=k;
					break;
				}
			}
			if(ipos==-1)
			{
				ipos=ptable.size();
				ptable.resize(ipos+1);
				ptable[ipos].id1=PosTable[i][j].id1;
				ptable[ipos].id2=PosTable[i][j].id2;
				ptable[ipos].count1=0;
				ptable[ipos].count2=0;
			}
			ptable[ipos].count1++;
			if(PosTable[i][j].dist<=critical_distance) 
			{
				ptable[ipos].count2++;
				icount++;
				Acor[0][j]=1;
			}
			else Acor[0][j]=0;
			
		}
		if(icount*1.0/f.size()>=critical_percent)
		{
			atom a,b;
			int atom_pos1,atom_pos2;
			atom_pos1=posVector[pos_id]-1;
			atom_pos2=negVector[neg_id]-1;
			a = table[posVector[pos_id]-1];
			b = table[negVector[neg_id]-1];
			Log("Found salt "+printAtom3(a)+" "+printAtom3(b)+" "+printPercent(icount*1.0/f.size()));

			string filename="salt2_"+printAtom2(a)+"_"+printAtom2(b)+".dat";
			string smoothname="salt2_"+printAtom2(a)+"_"+printAtom2(b)+".sda";
			string statname="salt2_"+printAtom2(a)+"_"+printAtom2(b)+".stat";
			string histname="salt2_"+printAtom2(a)+"_"+printAtom2(b)+".hist";
			string acorname="salt2_"+printAtom2(a)+"_"+printAtom2(b)+".cor";
                        PrintAcor(acorname,f,Acor[0],kstep,astep);

			FILE *fp=fopen(filename.c_str(),"w");
			if(!fp) Error(WriteError(filename));
			DoubleVector Distance;
			Distance.resize(f.size());
			PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
				printDistanceHeader("Dist")+" "+printString("Atom1",5)+" "+
				printString("Atom2",5));
			for(int j=0;j<f.size();j++)
			{
				Distance[j]=PosTable[i][j].dist;
				if(Distance[j]<=critical_distance) myCounter.add(f[j],"",1);
				Print(fp,printFrame(f[j])+" "+printDistance(PosTable[i][j].dist)+" ");
				PrintLine(fp,printString(table[PosTable[i][j].id1-1].atom_name,5)+" "+
				  printString(table[PosTable[i][j].id2-1].atom_name,5));
			}
			fclose(fp);
			if(smooth_flag)
			{
				IntVector sf;
				DoubleVector sd;
				makeSmooth(f,Distance,sf,sd,smooth_start,smooth_step);
				PrintDistances(smoothname,sf,sd);
			}
			string header1="#Distance statistics";
			PrintStatFile(statname,"w",f.size()/statcount,header1,f,Distance);
			Data dmin,dmax,davg,dstd;

			string header3="#Existence statistics";
			string header4="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Count");
			PrintIntStatFile(statname,"a",f.size()/statcount,header3,header4,f,Acor);

			fp=fopen(statname.c_str(),"a");
			string header2="\n#Atom statistics";
			PrintLine(fp,header2);
			PrintLine(fp,"#"+printString("Atom1",4)+" "+printString("Atom2",5)+" "+
				printFrameHeader("Frames")+" "+printPercentHeader("Percent")+" "+
				printFrameHeader("Frames")+" "+printPercentHeader("Percent"));
			for(int j=0;j<ptable.size();j++)
			{
			atom a,b;
			int atom_pos1,atom_pos2;
			atom_pos1=ptable[j].id1-1;
			atom_pos2=ptable[j].id2-1;
			a=table[atom_pos1];
			b=table[atom_pos2];

			Data p=ptable[j].count1*1.0/f.size();

			IntVector atom_list;
			atom_list.resize(2);
			atom_list[0]=atom_pos1+1;
			atom_list[1]=atom_pos2+1;
		
			myTex.addAtomList(atom_list,p);
			PrintLine(fp,printString(table[ptable[j].id1-1].atom_name,5)+" "+
		           printString(table[ptable[j].id2-1].atom_name,5)+" "+
			 printFrame(ptable[j].count1)+" "+
			printPercent(ptable[j].count1*1.0/f.size())+" "+
				printFrame(ptable[j].count2)+" "+
				printPercent(ptable[j].count2*1.0/f.size()));
			}
			fclose(fp);
			if(histflag)
			{
			vector<HistStruct> st_dist;
			Data dmin,dmax,davg,dstd;
			getVectorStatistics(Distance,dmin,dmax,davg,dstd);
			makeHist(Distance,st_dist,dmin,dmax,bindist);
			PrintHist(histname,st_dist);
			}
		}
	}
	myTex.print();
	myCounter.print("salt2_count.dat");
	myCounter.printStat("salt2_count.stat");
	if(smooth_flag) myCounter.printSmooth("salt2_count.sda");
	if(histflag)    myCounter.printHist("salt2_count.hist");
	
}
