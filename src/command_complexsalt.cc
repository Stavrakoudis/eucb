# include <global.h>
# include <pdb.h>
# include <command_complexsalt.h>
# include <tex.h>
# include <counter.h>

CommandComplexSalt::CommandComplexSalt()
	:Command("COMPLEXSALT")
{
	possalt.resize(6);
	negsalt.resize(4);
	possalt[0].res_name="LYS";possalt[0].atom_name="NZ";
	possalt[1].res_name="ARG";possalt[1].atom_name="NE";
	possalt[2].res_name="ARG";possalt[2].atom_name="NH1";
	possalt[3].res_name="ARG";possalt[3].atom_name="NH2";
	possalt[4].res_name="HSP";possalt[4].atom_name="ND1";
	possalt[5].res_name="HSP";possalt[5].atom_name="NE2";

	negsalt[0].res_name="ASP";negsalt[0].atom_name="OD1";
	negsalt[1].res_name="ASP";negsalt[1].atom_name="OD2";
	negsalt[2].res_name="GLU";negsalt[2].atom_name="OE1";
	negsalt[3].res_name="GLU";negsalt[3].atom_name="OE2";
}

typedef struct
{
	int id1;
	int id2;
	int id3;
	int count1;
	int count2;
}Triple;

/*	Implement the complex salt facility.
 * */
void CommandComplexSalt::Run()
{
	Counter myCounter;

	Tex	myTex("salt3.tex");
	vector<string> header;
	header.resize(3);
	header[0]="Positive1";
	header[1]="Negative";
	header[2]="Positive2";
	myTex.setHeader(header);
	myTex.setCaption("Complex salt bridges");

	Team *salt_team = team;
	vector<int> posVector;
	vector<int> negVector;
	vector<int> posVector2;
	vector<int> negVector2;
	for(int i=0;i<table.size();i++)
	{
		string chain;
		if(salt_team && !salt_team->find(table[i],chain)) continue;
		for(int j=0;j<possalt.size();j++)
			if(table[i].atom_type == "NH3" || 
			  (table[i].res_name == possalt[j].res_name && 
                           table[i].atom_name == possalt[j].atom_name))
			{
				int index=i;
				while(index>=0 && table[index].res_id == table[i].res_id) 
				{
					index--;
				}
				if(isin(posVector2,index+2)==-1)
				posVector2.push_back(index+2);
			}
		for(int j=0;j<negsalt.size();j++)
			if(table[i].res_name == negsalt[j].res_name && 
                           table[i].atom_name == negsalt[j].atom_name)
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
	vector<NoeStruct> noe;
	for(int i=0;i<posVector2.size();i++)
	{
		for(int j=0;j<negVector2.size();j++)
		{
			NoeStruct n;
			n.atom1=posVector2[i];
			n.atom2=negVector2[j];
			noe.push_back(n);
		}
	}
	dcd->setFrameStep(smart_skip);
	dcd->getNoe(f,noe);
	int icount=0;
	for(int i=0;i<posVector2.size();i++)
	{
		for(int j=0;j<negVector2.size();j++)
		{
			for(int k=0;k<f.size();k++)
			{
				double d=noe[icount].D[k];
				if(d<=smart_distance)
				{
					int ipos=isin(posVector,posVector2[i]);
					if(ipos==-1) posVector.push_back(posVector2[i]);
					ipos=isin(negVector,negVector2[j]);
					if(ipos==-1) negVector.push_back(negVector2[j]);
					break;
				}
			}
			icount++;
		}
	}

	vector<PosNegPosTable> PosTable;
	dcd->setFrameStep(step);
	dcd->getPosNegPos(f,posVector,negVector,possalt,negsalt,PosTable);
	Log("#Critical distance = "+printDistance(critical_distance));
	Log("#Critical percent  = "+printPercent(critical_percent));
	IntVector2 Acor;
	Acor.resize(3);
	Acor[0].resize(f.size());
	Acor[1].resize(f.size());
	Acor[2].resize(f.size());
	for(int i=0;i<PosTable.size();i++)
	{
		int pos_id1 = PosTable[i][0].pos_id1;
		int neg_id =  PosTable[i][0].neg_id;
		int pos_id2 = PosTable[i][0].pos_id2;
		if(pos_id1 == pos_id2) continue;
		int icount = 0;
		int icount1 = 0;
		int icount2 = 0;
		vector<Triple> ptable;
		
		
		for(int j=0;j<f.size();j++) 
		{
			
			int ipos=-1;
			for(int k=0;k<ptable.size();k++)
			{
				if(ptable[k].id1 == PosTable[i][j].id1 
				&& ptable[k].id2 == PosTable[i][j].id2
				&& ptable[k].id3 == PosTable[i][j].id3)
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
				ptable[ipos].id3=PosTable[i][j].id3;
				ptable[ipos].count1=0;
				ptable[ipos].count2=0;
			}
			ptable[ipos].count1++;
			

			if(PosTable[i][j].dist1<=critical_distance) icount1++;
			if(PosTable[i][j].dist2<=critical_distance) icount2++;
			if(PosTable[i][j].dist1<=critical_distance && PosTable[i][j].dist2<=critical_distance) 
			{
				ptable[ipos].count2++;
				icount++;
			}
			
		}
	
		if(icount*1.0/f.size()>=critical_percent)
		{
			int atom_pos1,atom_pos2,atom_pos3;
			string filename,statname,histname,smoothname,acorname;	
			atom_pos1=posVector[pos_id1]-1;
			atom_pos2=negVector[neg_id]-1;
			atom_pos3=posVector[pos_id2]-1;
			atom a = table[posVector[pos_id1]-1];
			atom b = table[negVector[neg_id]-1];
			atom c = table[posVector[pos_id2]-1];
			if(c.res_id == a.res_id) continue;
			Log(printAtom3(a)+" "+printAtom3(b)+" "+printAtom3(c)+" "+printPercent(icount*1.0/f.size()));
			filename="salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+printAtom2(c)+".dat";
			statname="salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+printAtom2(c)+".stat";
			histname="salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+printAtom2(c)+".hist";
			smoothname="salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+printAtom2(c)+".sda";
			acorname="salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+printAtom2(c)+".cor";

			FILE *fp=fopen(filename.c_str(),"w");
			if(!fp) Error(WriteError(filename));
			DoubleVector Distance1;
			Distance1.resize(f.size());
			DoubleVector Distance2;
			Distance2.resize(f.size());
			for(int k1=0;k1<3;k1++) 
				for(int k2=0;k2<f.size();k2++) Acor[k1][k2]=0;
			for(int j=0;j<f.size();j++)
			{
				Distance1[j]=PosTable[i][j].dist1;
				Distance2[j]=PosTable[i][j].dist2;

				Print(fp,printFrame(f[j])+" "+printDistance(PosTable[i][j].dist1)+" "+
					printDistance(PosTable[i][j].dist2)+" ");
				PrintLine(fp,printString(table[PosTable[i][j].id1-1].atom_name,5)+" "+
				  printString(table[PosTable[i][j].id2-1].atom_name,5)+" "+
				  printString(table[PosTable[i][j].id3-1].atom_name,5));
				if(PosTable[i][j].dist1<=critical_distance && PosTable[i][j].dist2<=
					critical_distance)
				{
					myCounter.add(f[j],printAtom2(table[PosTable[i][j].id2-1]),1);
					Acor[0][j]=1;
				}
				else
				{
					myCounter.add(f[j],"",0);
					Acor[0][j]=0;
				}
				if(PosTable[i][j].dist1<=critical_distance)
					Acor[1][j]=1; 
				else
					Acor[1][j]=0;
				if(PosTable[i][j].dist2<=critical_distance)
					Acor[2][j]=1;
				else
					Acor[2][j]=0;
			}
			fclose(fp);
			 PrintAcor(acorname,f,Acor[0],kstep,astep);

			if(smooth_flag)
			{
				IntVector sf;
				DoubleVector sd1,sd2;
				makeSmooth(f,Distance1,sf,sd1,smooth_start,smooth_step);
				makeSmooth(f,Distance2,sf,sd2,smooth_start,smooth_step);
				fp=fopen(smoothname.c_str(),"w");
				if(!fp) Error(WriteError(smoothname));
				for(int j=0;j<sf.size();j++)
					PrintLine(fp,printFrame(sf[j])+" "+
					printDistance(sd1[j])+" "+printDistance(sd2[j]));
				fclose(fp);
			}
			
			Data dmin1,dmax1,davg1,dstd1;
			Data dmin2,dmax2,davg2,dstd2;
			getVectorStatistics(Distance1,dmin1,dmax1,davg1,dstd1);
			getVectorStatistics(Distance2,dmin2,dmax2,davg2,dstd2);
			string header1="#Distance statistics "+printAtom5(a)+"-"+printAtom5(b);
			string header2="#Distance statistics "+printAtom5(b)+"-"+printAtom5(c);
			string header3="#Existence statistics";
			string header4="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Total")+" "+printFrameHeader("Count1")+" "+
				printFrameHeader("Count2");
			PrintStatFile(statname,"w",f.size()/statcount,header1,f,Distance1);	
			PrintStatFile(statname,"a",f.size()/statcount,header2,f,Distance2);	
			PrintIntStatFile(statname,"a",f.size()/statcount,header3,header4,f,Acor);

			fp=fopen(statname.c_str(),"a");
			PrintLine(fp,"#Atom statistics");
			for(int j=0;j<ptable.size();j++)
			{
			int atom_pos1,atom_pos2,atom_pos3;
			atom_pos1=ptable[j].id1-1;
			atom_pos2=ptable[j].id2-1;
			atom_pos3=ptable[j].id3-1;
			atom a,b,c;
			a=table[atom_pos1];
			b=table[atom_pos2];
			c=table[atom_pos3];
			IntVector atom_pos;
			atom_pos.resize(3);
			atom_pos[0]=atom_pos1+1;
			atom_pos[1]=atom_pos2+1;
			atom_pos[2]=atom_pos3+1;
			Data p=ptable[j].count1*1.0/f.size();
			myTex.addAtomList(atom_pos,p);

			PrintLine(fp,printString(table[ptable[j].id1-1].atom_name,5)+" "+
				             printString(table[ptable[j].id2-1].atom_name,5)+" "+
					     printString(table[ptable[j].id3-1].atom_name,5)+" "+
					 printFrame(ptable[j].count1)+" "+
					printPercent(ptable[j].count1*1.0/f.size())+" "+
					printFrame(ptable[j].count2)+" "+
					printPercent(ptable[j].count2*1.0/f.size()));
			}
			fclose(fp);
			
			if(histflag)
			{
			vector<HistStruct> st_dist1;
			vector<HistStruct> st_dist2;
			Data dmin1,dmax1,davg1,dstd1;
			Data dmin2,dmax2,davg2,dstd2;
			getVectorStatistics(Distance1,dmin1,dmax1,davg1,dstd1);
			getVectorStatistics(Distance2,dmin2,dmax2,davg2,dstd2);
			if(dmin1<=dmin2) dmin2=dmin1; else dmin1=dmin2;
			if(dmax1>=dmax2) dmax2=dmax1; else dmax1=dmax2;
			makeHist(Distance1,st_dist1,dmin1,dmax1,bindist);
			makeHist(Distance2,st_dist2,dmin2,dmax2,bindist);
			FILE *fout = fopen(histname.c_str(),"w");
			if(!fout) Error(WriteError(histname));
			for(int l=0;l<st_dist1.size();l++)
				PrintLine(fout,printNumber(
					st_dist1[l].value,8,2)+" "+
					printPercent(st_dist1[l].percent)+" "+
					printNumber(st_dist1[l].count,5)+" "+
					printNumber(st_dist2[l].value,8,2)+" "+
					printPercent(st_dist2[l].percent)+" "+
					printNumber(st_dist2[l].count,5));
			fclose(fout);
			
			}
		}
		
	}
	vector<NegPosNegTable> PosTable2;
	dcd->getNegPosNeg(f,posVector,negVector,possalt,negsalt,PosTable2);
	for(int i=0;i<PosTable2.size();i++)
	{
		int pos_id1 = PosTable2[i][0].neg_id1;
		int neg_id =  PosTable2[i][0].pos_id;
		int pos_id2 = PosTable2[i][0].neg_id2;
		if(pos_id1 == pos_id2) continue;
		int icount = 0;
		int icount1 = 0;
		int icount2 = 0;
		vector<Triple> ptable;
		
		for(int k1=0;k1<3;k1++) 
			for(int k2=0;k2<f.size();k2++) Acor[k1][k2]=0;
		
		for(int j=0;j<f.size();j++) 
		{
			
			int ipos=-1;
			for(int k=0;k<ptable.size();k++)
			{
				if(ptable[k].id1 == PosTable2[i][j].id1 
				&& ptable[k].id2 == PosTable2[i][j].id2
				&& ptable[k].id3 == PosTable2[i][j].id3)
				{
					ipos=k;
					break;
				}
			}
			if(ipos==-1)
			{
				ipos=ptable.size();
				ptable.resize(ipos+1);
				ptable[ipos].id1=PosTable2[i][j].id1;
				ptable[ipos].id2=PosTable2[i][j].id2;
				ptable[ipos].id3=PosTable2[i][j].id3;
				ptable[ipos].count1=0;
				ptable[ipos].count2=0;
			}
			ptable[ipos].count1++;
			
			if(PosTable2[i][j].dist1<=critical_distance) 
			{
				icount1++;
				Acor[1][j]=1;
			}
			else	Acor[1][j]=0;
			if(PosTable2[i][j].dist2<=critical_distance) 
			{
				icount2++;
				Acor[2][j]=1;
			}
			else
				Acor[2][j]=0;
			if(PosTable2[i][j].dist1<=critical_distance && PosTable2[i][j].dist2<=critical_distance) 
			{
				ptable[ipos].count2++;
				icount++;
				Acor[0][j]=1;
			}
			else	
				Acor[0][j]=0;
			
		}
	
		if(icount*1.0/f.size()>=critical_percent)
		{
			atom a = table[negVector[pos_id1]-1];
			atom b = table[posVector[neg_id]-1];
			atom c = table[negVector[pos_id2]-1];
			if(a.res_id == b.res_id || b.res_id == c.res_id 
			|| c.res_id == a.res_id) continue;
			string filename = "salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+
			  		printAtom2(c)+".dat";
			string statname = "salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+
			  		printAtom2(c)+".stat";
			string smoothname = "salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+
			  		printAtom2(c)+".sda";
			string histname = "salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+
			  		printAtom2(c)+".hist";
			string acorname = "salt3_"+printAtom2(a)+"_"+printAtom2(b)+"_"+
			  		printAtom2(c)+".cor";
			PrintAcor(acorname,f,Acor[0],kstep,astep);
			Log(printAtom3(a)+" "+printAtom3(b)+" "+printAtom3(c)+" "+printPercent(icount*1.0/f.size()));

			FILE *fp=fopen(filename.c_str(),"w");
			if(!fp) Error(WriteError(filename));
			DoubleVector Distance1;
			Distance1.resize(f.size());
			DoubleVector Distance2;
			Distance2.resize(f.size());
			for(int j=0;j<f.size();j++)
			{
				Distance1[j]=PosTable2[i][j].dist1;
				Distance2[j]=PosTable2[i][j].dist2;

				Print(fp,printFrame(f[j])+" "+printDistance(Distance1[j])+" "+
					printDistance(Distance2[j])+" ");
				PrintLine(fp,printString(table[PosTable2[i][j].id1-1].atom_name,5)+" "+
				  printString(table[PosTable2[i][j].id2-1].atom_name,5)+" "+
				  printString(table[PosTable2[i][j].id3-1].atom_name,5));
			}
			fclose(fp);
			if(smooth_flag)
			{
				IntVector sf;
				DoubleVector sd1,sd2;
				makeSmooth(f,Distance1,sf,sd1,smooth_start,smooth_step);
				makeSmooth(f,Distance2,sf,sd2,smooth_start,smooth_step);
				fp=fopen(smoothname.c_str(),"w");
				if(!fp) Error(WriteError(smoothname));
				for(int j=0;j<sf.size();j++)
					PrintLine(fp,printFrame(sf[j])+" "+printDistance(sd1[j])+
					" "+printDistance(sd2[j]));
				fclose(fp);
			}
			
			Data dmin1,dmax1,davg1,dstd1;
			Data dmin2,dmax2,davg2,dstd2;
			getVectorStatistics(Distance1,dmin1,dmax1,davg1,dstd1);
			getVectorStatistics(Distance2,dmin2,dmax2,davg2,dstd2);
			string header1="#Distance statistics "+printAtom5(a)+"-"+printAtom5(b);
			string header2="#Distance statistics "+printAtom5(b)+"-"+printAtom5(c);
			string header3="#Existence statistics";
			string header4="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Total")+" "+printFrameHeader("Count1")+" "+
				printFrameHeader("Count2");
			PrintStatFile(statname,"w",f.size()/statcount,header1,f,Distance1);	
			PrintStatFile(statname,"a",f.size()/statcount,header2,f,Distance2);	
			PrintIntStatFile(statname,"a",f.size()/statcount,header3,header4,f,Acor);

			fp=fopen(statname.c_str(),"a");
			PrintLine(fp,"#Atom statistics");
			PrintLine(fp,printFrame(icount1)+" "+printPercent(icount1*1.0/f.size()));
			PrintLine(fp,printFrame(icount2)+" "+printPercent(icount2*1.0/f.size()));
			PrintLine(fp,printFrame(icount)+" "+printPercent(icount*1.0/f.size()));
			for(int j=0;j<ptable.size();j++)
			{
				PrintLine(fp,printString(table[ptable[j].id1-1].atom_name,5)+" "+
				             printString(table[ptable[j].id2-1].atom_name,5)+" "+
					     printString(table[ptable[j].id3-1].atom_name,5)+" "+
					 printFrame(ptable[j].count1)+" "+
					printPercent(ptable[j].count1*1.0/f.size())+" "+
					printFrame(ptable[j].count2)+" "+
					printPercent(ptable[j].count2*1.0/f.size()));
			}
			fclose(fp);
			
			if(histflag)
			{
			vector<HistStruct> st_dist1;
			vector<HistStruct> st_dist2;
			Data dmin1,dmax1,davg1,dstd1;
			Data dmin2,dmax2,davg2,dstd2;
			getVectorStatistics(Distance1,dmin1,dmax1,davg1,dstd1);
			getVectorStatistics(Distance2,dmin2,dmax2,davg2,dstd2);
			if(dmin1<=dmin2) dmin2=dmin1; else dmin1=dmin2;
			if(dmax1>=dmax2) dmax2=dmax1; else dmax1=dmax2;
			makeHist(Distance1,st_dist1,dmin1,dmax1,bindist);
			makeHist(Distance2,st_dist2,dmin2,dmax2,bindist);
			FILE *fout = fopen(histname.c_str(),"w");
			if(!fout) Error(WriteError(histname));
			for(int l=0;l<st_dist1.size();l++)
				PrintLine(fout,printDistance(st_dist1[l].value)+" "+
					printPercent(st_dist1[l].percent)+" "+
					printFrame(st_dist1[l].count)+" "+
					printPercent(st_dist2[l].value)+" "+
					printPercent(st_dist2[l].percent)+" "+
					printFrame(st_dist2[l].count));
			fclose(fout);
			
			}
		}
		
	}
	myTex.print();
	myCounter.print("salt3_count.dat");
	myCounter.printStat("salt3_count.stat");
	if(smooth_flag) myCounter.printSmooth("salt3_count.sda");
	if(histflag)    myCounter.printHist("salt3_count.hist");
}
