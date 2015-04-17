# include <global.h>
# include <command_sidedist.h>
# include <tex.h>
# include <counter.h>

typedef struct
{
	int atom1;
	int atom2;
	int count1;
	int count2;
}PairAtom;

CommandSideDist::CommandSideDist()
	:Command("SIDEDIST")
{
}

void	CommandSideDist::setName1(string s)
{
	if(s=="aromatic") team1=makeAromatic();
	else
	if(s=="aliphatic") team1=makeAliphatic();
	else
	if(s=="positive") team1=makePositive();
	else
	if(s=="negative") team1=makeNegative();
	else
	if(s=="hydrophobic") team1=makeHydrophobic();
	else
	{
		Log("Unknown keyword for side1 ");
		psf_error("Unknown keyword for side1 ");
	}
}

void	CommandSideDist::setName2(string s)
{
	if(s=="aromatic") team2=makeAromatic();
	else
	if(s=="aliphatic") team2=makeAliphatic();
	else
	if(s=="positive") team2=makePositive();
	else
	if(s=="negative") team2=makeNegative();
	else
	if(s=="hydrophobic") team2=makeHydrophobic();
	else
	{
		Log("Unknown keyword for side2 ");
		psf_error("Unknown keyword for side2 ");
	}
}

void CommandSideDist::Run()
{
	vector<int> list1;
	vector<int> list2;
	Counter myCounter;
	Tex	myTex("side.tex");
	vector<string>  s; s.resize(2);
	s[0]="Side1";
	s[1]="Side2";
	myTex.setHeader(s);
	myTex.setCaption("Sidechain");
	myTex.enableAtoms(0);
	
	Team *sidedist_team1=team1;
	Team *sidedist_team2=team2;
	if(sidedist_team1==NULL || sidedist_team2==NULL)  Error(SeqError("sidedist"));
	sidedist_team1->enumerateAtoms(table,list1);
	sidedist_team2->enumerateAtoms(table,list2);
	Log("Critical distance = "+printDistance(critical_distance));
	Log("Critical percent  = "+printPercent(critical_percent));

	IntVector  frame;
	DoubleVector2 distance;
	IntVector2 group1;
	IntVector2 group2;
	IntVector2 atom1;
	IntVector2 atom2;
	
	for(int i=0;i<list1.size();i++)
	{
		IntVector group1_1;
		IntVector group2_1;
		group1_1.resize(0);
		int start_id = table[list1[i]].res_id;
		int start_I = i;
		while(table[list1[i]].res_id==start_id) 
		{
			group1_1.push_back(list1[i]+1);
			i++;
		}
		i--;
		for(int j=0;j<list2.size();j++)
		{
			group2_1.resize(0);
			int start2_id = table[list2[j]].res_id;
			int start_J = j;
			if(start_I == start_J) continue;
			if(table[list1[start_I]].res_id == start2_id 
			&& table[list1[start_I]].chain_id == 
		           table[list2[start_J]].chain_id) continue;
			while(table[list2[j]].res_id==start2_id)
			{
				group2_1.push_back(list2[j]+1);
				j++;
			}
			group1.push_back(group1_1);
			group2.push_back(group2_1);
			j--;
		}
	}
	
	dcd->setFrameStep(smart_skip);
	dcd->getGroupDistance(group1,group2,frame,distance,atom1,atom2);
	IntVector2 Group1;
	IntVector2 Group2;
	for(int i=0;i<group1.size();i++)
	{
		int frames_below=getBelowCount(distance[i],smart_distance);
		if(frames_below)
		{
			Group1.push_back(group1[i]);
			Group2.push_back(group2[i]);
		}
	}
	dcd->setFrameStep(step);
	dcd->getGroupDistance(Group1,Group2,frame,distance,atom1,atom2);
	for(int i=0;i<Group1.size();i++)
	{
		Data pt=getBelowPercent(distance[i],critical_distance);
		int frames_below=getBelowCount(distance[i],critical_distance);
		if(pt<critical_percent) continue;
		
		atom Atom1=table[Group1[i][0]-1];
		atom Atom2=table[Group2[i][0]-1];
	

		IntVector atom_pos;
		atom_pos.resize(2);
		atom_pos[0]=Group1[i][0];
		atom_pos[1]=Group2[i][0];
		myTex.addAtomList(atom_pos,pt);

		Log("Found side "+printAtom3(Atom1)+" "+printAtom3(Atom2)+" "+printPercent(pt));
		string filename,statname,histname,smoothname,acorname;
		if(Atom1.res_id<=Atom2.res_id)
		{
			filename="side_"+printAtom2(Atom1)+"_"+printAtom2(Atom2)+".dat";
			statname="side_"+printAtom2(Atom1)+"_"+printAtom2(Atom2)+".stat";
			histname="side_"+printAtom2(Atom1)+"_"+printAtom2(Atom2)+".hist";
			smoothname="side_"+printAtom2(Atom1)+"_"+printAtom2(Atom2)+".sda";
			acorname="side_"+printAtom2(Atom1)+"_"+printAtom2(Atom2)+".cor";
		}
		else
		{
			filename="side_"+printAtom2(Atom2)+"_"+printAtom2(Atom1)+".dat";
			statname="side_"+printAtom2(Atom2)+"_"+printAtom2(Atom1)+".stat";
			histname="side_"+printAtom2(Atom2)+"_"+printAtom2(Atom1)+".hist";
			smoothname="side_"+printAtom2(Atom2)+"_"+printAtom2(Atom1)+".sda";
			acorname="side_"+printAtom2(Atom2)+"_"+printAtom2(Atom1)+".cor";
		}
		FILE *fp=fopen(filename.c_str(),"w");
		if(!fp) Error(WriteError(filename));
		vector<PairAtom> pair;
		PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
			printDistanceHeader("Dist")+" "+
			printString("Atom1",5)+" "+printString("Atom2",5));
		IntVector2 Acor;
		Acor.resize(1);
		Acor[0].resize(frame.size());
		for(int k=0;k<frame.size();k++)
		{
			PrintLine(fp,printFrame(frame[k])+" "+
			  printDistance(distance[i][k])+" "+
			  printString(table[atom1[i][k]-1].atom_name,5)+" "+
			  printString(table[atom2[i][k]-1].atom_name,5));
			if(distance[i][k]<=critical_distance)
			{
				myCounter.add(frame[k],printAtom2(table[atom1[i][k]-1]),1);	
				Acor[0][k]=1;
			}
			else
			{
				myCounter.add(frame[k],"",0);
				Acor[0][k]=0;
			}
			int ifound=0;
			for(int l=0;l<pair.size();l++)
			{
				if(pair[l].atom1==atom1[i][k] && pair[l].atom2==atom2[i][k])
				{
					ifound=l+1;
					break;
					}
				}
				if(!ifound)
				{
					PairAtom t;
					t.atom1=atom1[i][k];
					t.atom2=atom2[i][k];
					t.count1=1;
					if(distance[i][k]<=critical_distance) t.count2=1;
					else t.count2=0;
					pair.push_back(t);
				}
				else
				{
					pair[ifound-1].count1++;
					if(distance[i][k]<=critical_distance)
						pair[ifound-1].count2++;
				}
			}
			fclose(fp);
			if(smooth_flag)
			{
				IntVector sf;
				DoubleVector sd;
				makeSmooth(frame,distance[i],sf,sd,smooth_start,smooth_step);
				PrintDistances(smoothname,sf,sd);
			}
                        PrintAcor(acorname,frame,Acor[0],kstep,astep);

			string header1="#Distance statistics";
			PrintStatFile(statname,"w",frame.size()/statcount,header1,frame,distance[i]);
			string header2="#Existence statistics";
			string header3="#"+printString("Block",FRAME_WIDTH-1)+" "+
				printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" "+
				printFrameHeader("Count");
			
			PrintIntStatFile(statname,"a",frame.size()/statcount,header2,header3,frame,Acor);

			fp=fopen(statname.c_str(),"a");
			PrintLine(fp,"#Atom statistics");
			PrintLine(fp,printString("Atom1",5)+" "+printString("Atom2",5)+" "+
				printFrameHeader("Frames")+" "+printPercentHeader("Percent")+" "+
				printFrameHeader("Frames")+" "+printPercentHeader("Percent"));
			for(int k=0;k<pair.size();k++)
			{
				Data p1=pair[k].count1*1.0/frame.size();
				Data p2=pair[k].count2*1.0/frame.size();
				PrintLine(fp,printString(table[pair[k].atom1-1].atom_name,5)+" "+
					printString(table[pair[k].atom2-1].atom_name,5)+" "+
					printFrame(pair[k].count1)+" "+printPercent(p1)+" "+
					printFrame(pair[k].count2)+" "+printPercent(p2));
			}
			fclose(fp);

			if(histflag)
			{
				Data dmin,dmax,davg,dstd;
				getVectorStatistics(distance[i],dmin,dmax,davg,dstd);
				vector<HistStruct> st;
				makeHist(distance[i],st,dmin,dmax,bindist);
				PrintHist(histname,st);
			}
	}
	myTex.print();
	if(myCounter.size()==0) return;
	myCounter.print("side_count.dat");
	myCounter.printStat("side_count.stat");
	if(smooth_flag)
		myCounter.printSmooth("side_count.sda");
	if(histflag)
		myCounter.printHist("side_count.hist");
}


CommandSideDist::~CommandSideDist()
{
}
