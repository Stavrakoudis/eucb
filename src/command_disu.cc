# include <command_disu.h>
# include <global.h>
# include <pdb.h>
CommandDisu::CommandDisu()
	:Command("DISU")
{
}

void	CommandDisu::Run()
{
	IntVector Index1,Index2,Index3,Index4;	
	string chain;
	for(int i=0;i<table.size();i++)
	{
		if(!team || team->find(table[i],chain))
		{
			if(table[i].res_name=="CYS" && table[i].atom_name=="SG")
			{
				for(int j=0;j<table[i].bond.size();j++)
				{
					int id=table[i].bond[j];
					if(table[id-1].res_name=="CYS" && table[id-1].atom_name=="SG")
					{
						if(isin(Index2,id)!=-1) continue;
						for(int k=i-10;k<=i+10;k++)
						{
							if(k<0) k==0;
							if(k>=table.size()) k=table.size()-1;
							if(table[k].res_name=="CYS" && table[k].atom_name=="CB" && table[k].res_id==table[i].res_id)
						{Index1.push_back(k+1);break;}
						}
						for(int k=id-10;k<=id+10;k++)
						{
							if(k<0) k==0;
							if(k>=table.size()) k=table.size()-1;
							if(table[k].res_name=="CYS" && table[k].atom_name=="CB" && table[k].res_id==table[id-1].res_id)
						{Index4.push_back(k+1);break;}
						}
						Index2.push_back(i+1);
						Index3.push_back(id);
					}
				}
			}
		}	
	}	
	DoubleVector2 D,A;
	IntVector F;
	dcd->getDisu(F,Index1,Index2,Index3,Index4,D,A);
	string filename,statname,histnamea,histnamed,histname2;
	FILE *fptex=fopen("tableDISU.tex","w");

	PrintLine(fptex,"\\documentclass{article}");
        PrintLine(fptex,"\\begin{document}");
	PrintLine(fptex,"\\noindent");
	PrintLine(fptex,"\\textbf{Table X.}");
	PrintLine(fptex,"");
        PrintLine(fptex,"\\begin{table}");
        PrintLine(fptex,"\\caption{Disulfide results}");
	Print(fptex,"\\begin{tabular}{");
	if(pdbfile=="")
	for(int i=0;i<5;i++)
		Print(fptex," c ");
	else
	for(int i=0;i<7;i++)
		Print(fptex," c ");
	PrintLine(fptex,"}\n\\\\\\hline");
	vector<string> header;
	header.resize(7);
	header[0]="Cys";header[1]="Cys";
	header[2]="PDB-dist";header[3]="PDB-dihe";
	header[4]="MD-dist";header[5]="MD-dihe";header[6]="+/-";
	for(int i=0;i<7;i++)
	{
		if((i==2 || i==3 ) && pdbfile=="") continue;
		Print(fptex,"\\textbf{"+header[i]+"}");
		if(i!=6) Print(fptex,"&");
	}

	PrintLine(fptex,"\\\\\\hline");	
	for(int i=0;i<Index1.size();i++)
	{
		atom a = table[Index1[i]-1];
		atom b = table[Index4[i]-1];
		filename="disu_"+printAtom2(a)+"_"+printAtom2(b)+".dat";
		FILE *fp=fopen(filename.c_str(),"w");
		PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
			printDistanceHeader("Dist")+" "+
			printAngleHeader("Dihe"));
		int icount=0;
		for(int j=0;j<F.size();j++)
		{
			PrintLine(fp,printFrame(F[j])+" "+
				printDistance(D[i][j])+" "+
				printAngle(A[i][j]));
			if(j && A[i][j]*A[i][j-1]<0) icount++;
		}
		fclose(fp);
		statname="disu_"+printAtom2(a)+"_"+printAtom2(b)+".stat";
		string header1="#Distance statistics";
		PrintStatFile(statname,"w",F.size()/statcount,header1,F,D[i]);
		string header2="#Dihedral statistics";
		PrintAngleStatFile(statname,"a",F.size()/statcount,header2,F,A[i]);
		histnamed="disu_"+printAtom2(a)+"_"+printAtom2(b)+".histd";
		histnamea="disu_"+printAtom2(a)+"_"+printAtom2(b)+".hista";
		histname2="disu_"+printAtom2(a)+"_"+printAtom2(b)+".hist2";
		Print(fptex,a.chain_id+":"+printNumber(a.res_id)+"&"+
			b.chain_id+":"+printNumber(b.res_id)+"&");
		if(pdbfile!="")
		{
			PdbAtom pa,pb,pc,pd;
			pa=pdb_table[Index1[i]-1];
			pb=pdb_table[Index2[i]-1];
			pc=pdb_table[Index3[i]-1];
			pd=pdb_table[Index4[i]-1];
			
			Data dist=getPdbDistance(pa,pd);
			Data dihe=torsion(pa.x,pa.y,pa.z,
					  pb.x,pb.y,pb.z,
					  pc.x,pc.y,pc.z,
					  pd.x,pd.y,pd.z);	
			Print(fptex,printNumber(dist,7,1)+" & ");
			Print(fptex,printNumber(dihe,7,1)+" & ");
		}
		Data dmin,dmax,davg,dstd;
		getVectorStatistics(D[i],dmin,dmax,davg,dstd);
		Data amin,amax,aavg,astd;
		getAngleVectorStatistics(A[i],amin,amax,aavg,astd);
		Print(fptex,printNumber(davg,7,1)+"("+printNumber(dstd,7,1)+")&");
		Print(fptex,printNumber(aavg,7,1)+"("+printNumber(astd,7,1)+")&");
		Print(fptex,printNumber(icount));
			
		Print(fptex,"\\\\\n");
		if(histflag)
		{
			vector<HistStruct> st;
			makeHist(D[i],st,dmin,dmax,bindist);
			PrintHist(histnamed,st);

			vector<HistStruct> st_a;
			makeHist(A[i],st_a,amin,amax,bindihe);
			PrintAngleHist(histnamea,st_a);
			fp=fopen(histname2.c_str(),"w");
			for(int k=0;k<st.size();k++)
			{
				for(int m=0;m<st_a.size();m++)
				{
				int icount=0;
				for(int l=0;l<F.size();l++)
				{
					if(D[i][l]>=st[k].value-bindist/2.0 
					 && D[i][l]<=st[k].value+bindist/2.0 
					 && A[i][l]>=st_a[m].value-binangle/2.0
					 && A[i][l]<=st_a[m].value+binangle/2.0)
							icount++;
						}
					PrintLine(fp,printDistance(st[k].value)+" "+
					 printAngle(st_a[m].value)+" "+printFrame(icount));
					}
				   PrintLine(fp,"");
				}
				fclose(fp);

		}
	}
	PrintLine(fptex,"\\hline");
	PrintLine(fptex,"\\end{tabular}");
        PrintLine(fptex,"\\end{table}");
        PrintLine(fptex,"\\end{document}");
	fclose(fptex);
}

CommandDisu::~CommandDisu()
{
}
