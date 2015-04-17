# include <global.h>
# include <command_stack.h>

CommandStack::CommandStack(Data sp,Data sd,Data sa )
	:Command("STACK")
{
	stack_percent = sp;
	stack_distance = sd;
	stack_angle = sa;
	SideChain.resize(24);
	SideChain[0].res_name="ARG";SideChain[0].atom_name="NE";
	SideChain[1].res_name="ARG";SideChain[1].atom_name="NH1";
	SideChain[2].res_name="ARG";SideChain[2].atom_name="NH2";

	SideChain[3].res_name="TYR";SideChain[3].atom_name="CD1";
	SideChain[4].res_name="TYR";SideChain[4].atom_name="CD2";
	SideChain[5].res_name="TYR";SideChain[5].atom_name="CZ";

	SideChain[6].res_name="PHE";SideChain[6].atom_name="CD1";
	SideChain[7].res_name="PHE";SideChain[7].atom_name="CD2";
	SideChain[8].res_name="PHE";SideChain[8].atom_name="CZ";

	SideChain[9].res_name="TRP";SideChain[9].atom_name="CD1";
	SideChain[10].res_name="TRP";SideChain[10].atom_name="CE3";
	SideChain[11].res_name="TRP";SideChain[11].atom_name="CZ2";

	SideChain[12].res_name="HSD";SideChain[12].atom_name="CG";
	SideChain[13].res_name="HSD";SideChain[13].atom_name="CE1";
	SideChain[14].res_name="HSD";SideChain[14].atom_name="NE2";

	SideChain[15].res_name="HSE";SideChain[15].atom_name="CG";
	SideChain[16].res_name="HSE";SideChain[16].atom_name="CE1";
	SideChain[17].res_name="HSE";SideChain[17].atom_name="NE2";

	SideChain[18].res_name="HSP";SideChain[18].atom_name="CG";
	SideChain[19].res_name="HSP";SideChain[19].atom_name="CE1";
	SideChain[20].res_name="HSP";SideChain[20].atom_name="NE2";

	SideChain[21].res_name="PRO";SideChain[21].atom_name="CB";
	SideChain[22].res_name="PRO";SideChain[22].atom_name="CG";
	SideChain[23].res_name="PRO";SideChain[23].atom_name="CD";

}

void CommandStack::Run()
{
	Team *stack_team = team;
	vector<stack_table> stable;
	string chain;
	int icount=0;

	for(int i=0;i<table.size();i++)
	{
		for(int j=0;j<SideChain.size();j++)
		{
			if(table[i].res_name==SideChain[j].res_name 
			  && table[i].atom_name==SideChain[j].atom_name &&
			  (stack_team==NULL || stack_team->find(table[i],chain)))
			{
				icount++;
				if(icount==1) {int s=stable.size();stable.resize(s+1);}
				if(icount==1) 
				{
					stable[stable.size()-1].id1=i+1;
				}
				if(icount==2) 
				{
					stable[stable.size()-1].id2=i+1;
				}
				if(icount==3) 
				{
					stable[stable.size()-1].id3=i+1;
				}
				if(icount==3) icount=0;
			}
		}
	}
	
	vector<StackAtom> stack_atom;
	IntVector F;
	dcd->getStackAverages(F,stable,stack_atom);
	Log("Critical percent  = "+printPercent(stack_percent));
	Log("Critical distance = "+printDistance(stack_distance));
	Log("Critical angle    = "+printAngle(stack_angle));
	IntVector2 Acor;
	Acor.resize(3);
	Acor[0].resize(F.size());
	Acor[1].resize(F.size());
	Acor[2].resize(F.size());
	for(int i=0;i<stack_atom.size();i++)
	{
		int icount=0;
		int icount1=0;
		int icount2=0;
		for(int j=0;j<F.size();j++)
		{
			if(stack_atom[i].D[j]<=stack_distance) 
			{
				icount1++;
				Acor[1][j]=1;
			}
			else	Acor[1][j]=0;
			if(stack_atom[i].A[j]<=stack_angle) 
			{
				icount2++;
				Acor[2][j]=1;
			}
			else	Acor[2][j]=0;
			if(stack_atom[i].D[j]<=stack_distance && stack_atom[i].A[j]<=stack_angle)
			{
				icount++;
				Acor[0][j]=1;
			}
			else    Acor[0][j]=0;
		}
		if(icount*1.0/F.size()<stack_percent) continue;
		atom atom1,atom2;
		atom1=table[stack_atom[i].atom1_id1-1];
		atom2=table[stack_atom[i].atom2_id1-1];
		Log("Atom1 "+printAtom3(atom1)+" "+"Atom2" +printAtom3(atom2));
		string filename="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".dat";
		string statname="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".stat";
		string histnamed="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".histd";
		string histnamea="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".hista";
		string histname2="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".hist2";
		string smoothname="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".sda";
		string acorname="stack_"+printAtom2(atom1)+"_"+printAtom2(atom2)+".cor";
		PrintDistAndAngle(filename,F,stack_atom[i].D,stack_atom[i].A);
		if(smooth_flag)
		{
			IntVector sf;
			DoubleVector sd;
			DoubleVector sa;
			makeSmooth(F,stack_atom[i].D,sf,sd,smooth_start,smooth_step);
			makeAngleSmooth(F,stack_atom[i].A,sf,sa,smooth_start,smooth_step);
			PrintDistAndAngle(smoothname,sf,sd,sa);
		}
                PrintAcor(acorname,F,Acor[0],kstep,astep);

		string header1="#Statistics for the distance";
		PrintStatFile(statname,"w",F.size()/statcount,header1,F,stack_atom[i].D);
		string header2="#Statistics for the angle";
		PrintAngleStatFile(statname,"a",F.size()/statcount,header2,F,stack_atom[i].A);
		string header3="#Existence statistics";
		string header4="#"+printString("Block",FRAME_WIDTH-1)+" "+printFrameHeader("Frame1")+" "+
			printFrameHeader("Frame2")+" "+printFrameHeader("Total")+" "+
			printFrameHeader("Distance")+" "+printFrameHeader("Angle");
		PrintIntStatFile(statname,"a",F.size()/statcount,header3,header4,F,Acor);

		if(histflag)
		{
			vector<HistStruct> st_dist;
			double dmin,dmax;
			dmin=getMinimum(stack_atom[i].D,F.size());
			dmax=getMaximum(stack_atom[i].D,F.size());
			makeHist(stack_atom[i].D,st_dist,dmin,dmax,bindist);
			PrintHist(histnamed,st_dist);
			vector<HistStruct> st_angle;
			makeHist(stack_atom[i].A,st_angle,0.0,90.0,binangle);
			PrintAngleHist(histnamea,st_angle);
			FILE *fp=fopen(histname2.c_str(),"w");
			if(!fp) Error(WriteError(histname2));
			PrintLine(fp,"#"+
				printString("Dist",DISTANCE_WIDTH-1)+" "+printTorsionHeader("Tors")+" "+
				printFrameHeader("Common"));
			for(int k=0;k<st_dist.size();k++)
			{
				for(int m=0;m<st_angle.size();m++)
				{
					int ct=0;
					for(int l=0;l<stack_atom[i].D.size();l++)
					{
						Data vd=stack_atom[i].D[l];
						Data va=stack_atom[i].A[l];
						if(vd>=st_dist[k].value-bindist/2.0 &&
						vd<=st_dist[k].value+bindist/2.0 &&
						va>=st_angle[m].value-bindihe/2.0 && 
						va<=st_angle[m].value+bindihe/2.0) ct++;
					}
					PrintLine(fp,printDistance(st_dist[k].value)+" "+
				  	 printAngle(st_angle[m].value)+" "+
					 printFrame(ct));
				}
			}
			fclose(fp);
		}
	}
}
