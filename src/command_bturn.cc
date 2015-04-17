# include <command_bturn.h>
# include <global.h>
# include <math.h>

CommandBturn::CommandBturn(Data P,Data D,Data A)
	:Command("BTURN")
{
	Beta.resize(8);
	Beta[0].type="I";Beta[0].phi2=-60.0;Beta[0].psi2=-30.0;
	Beta[0].phi3=-90.0;Beta[0].psi3=0.0;Beta[0].distance=4.6;
	Beta[1].type="I'";Beta[1].phi2=60.0;Beta[1].psi2=30.0;
	Beta[1].phi3=90.0;Beta[1].psi3=0.0;Beta[1].distance=4.6;
	Beta[2].type="II";Beta[2].phi2=-60.0;Beta[2].psi2=120.0;
	Beta[2].phi3=80.0;Beta[2].psi3=0.0;Beta[2].distance=4.6;
	Beta[3].type="II'";Beta[3].phi2=60.0;Beta[3].psi2=-120.0;
	Beta[3].phi3=-80.0;Beta[3].psi3=0.0;Beta[3].distance=4.6;
	Beta[4].type="VIa1";Beta[4].phi2=-60.0;Beta[4].psi2=120.0;
	Beta[4].phi3=-90.0;Beta[4].psi3=0.0;Beta[4].distance=3.4;

	Beta[5].type="VIa2";Beta[5].phi2=-120.0;Beta[5].psi2=120.0;
	Beta[5].phi3=-60.0;Beta[5].psi3=0.0;Beta[5].distance=3.7;
	Beta[6].type="VIb";Beta[6].phi2=-135.0;Beta[6].psi2=135.0;
	Beta[6].phi3=-75.0;Beta[6].psi3=160.0;Beta[6].distance=6.0;
	Beta[7].type="VIII";Beta[7].phi2=-60.0;Beta[7].psi2=-30.0;
	Beta[7].phi3=-120.0;Beta[7].psi3=120.0;Beta[7].distance=6.3;

	bturn_percent = P;
	bturn_distance = D;
	bturn_angle = A;
}

/*	Implement the beta turn search facility.
 * */
void CommandBturn::Run()
{
	Team *bturn_team = team;
	if(bturn_team == NULL) return ;
	vector<string>  chain;
	bturn_team->enumerateChains(chain);

	vector<int> frame;
	vector<Data> distance;
	vector<Data> angle;
	vector<Data> hbd_distance;
	vector<Data> hba_angle;
	vector<Data> phi2_angle;
	vector<Data> psi2_angle;
	vector<Data> phi3_angle;
	vector<Data> psi3_angle;

	vector<NoeStruct>    Noe;
	vector<NoeAngleStruct> NoeAngle;
	vector<NoeDihedralStruct> NoeDihedral;

	//pass 1
	for(int i=0;i<chain.size();i++)
	{
		vector<int> res_id;
		bturn_team->enumerateResId(chain[i],res_id);
		if(res_id.size()<4)  continue;
		for(int j=0;j<=res_id.size()-4;j++)
		{
			int id1=res_id[j];
			int id2=res_id[j+1];
			int id3=res_id[j+2];
			int id4=res_id[j+3];
			if(id2-id1!=1 && id3-id2!=1 && id4-id3!=1) continue;
			int atom1_CA=-1,atom2_CA=-1,atom3_CA=-1,atom4_CA=-1;
			int atom1_O=-1,atom4_HN=-1,atom4_N=-1;
			int atom2_C=-1,atom2_N=-1,atom1_C=-1;
			int atom3_N=-1,atom3_C=-1;
			for(int k=0;k<table.size();k++)
			{
				string mychain;
				if(bturn_team->find(table[k],mychain)==0) continue;
				if(table[k].res_id==id1 && table[k].atom_name=="CA") 
					atom1_CA=table[k].atom_id;
				if(table[k].res_id==id2 && table[k].atom_name=="CA") 
					atom2_CA=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="CA") 
					atom3_CA=table[k].atom_id;
				if(table[k].res_id==id4 && table[k].atom_name=="CA") atom4_CA=table[k].atom_id;
				if(table[k].res_id==id4 && table[k].atom_name=="N")  atom4_N=table[k].atom_id;
				if(table[k].res_id==id4 && table[k].atom_name=="HN") atom4_HN=table[k].atom_id;
				if(table[k].res_id==id1 && table[k].atom_name=="O")  atom1_O=table[k].atom_id;	
				if(table[k].res_id==id2 && table[k].atom_name=="N")  atom2_N=table[k].atom_id;
				if(table[k].res_id==id2 && table[k].atom_name=="CA") atom2_CA=table[k].atom_id;
				if(table[k].res_id==id2 && table[k].atom_name=="C")  atom2_C=table[k].atom_id;
				if(table[k].res_id==id1 && table[k].atom_name=="C")  atom1_C=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="N")  atom3_N=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="CA") atom3_CA=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="C")  atom3_C=table[k].atom_id;
			}
			if(atom1_CA==-1 || atom2_CA==-1 || atom3_CA==-1 || atom4_CA==-1) continue;
			//perform queries
			NoeStruct noe;
			noe.atom1=atom1_CA;
			noe.atom2=atom4_CA;
			Noe.push_back(noe);

			NoeAngleStruct noea;
			NoeDihedralStruct noed;
			noed.atom1=atom1_CA;
			noed.atom2=atom2_CA;
			noed.atom3=atom3_CA;
			noed.atom4=atom4_CA;
			NoeDihedral.push_back(noed);

			//printf header
			if(atom4_HN!=-1 && atom1_O!=-1)
			{
				noe.atom1=atom4_HN;
				noe.atom2=atom1_O;
				Noe.push_back(noe);
			}
			if(atom4_N!=-1 && atom1_O!=-1 && atom4_HN!=-1) 	
			{
				noea.atom1=atom4_N;
				noea.atom2=atom4_HN;
				noea.atom3=atom1_O;
				NoeAngle.push_back(noea);
			}
			if(atom1_C!=-1 && atom2_N!=-1 && atom2_CA!=-1 && atom2_C!=-1)
			{
				noed.atom1=atom1_C;
				noed.atom2=atom2_N;
				noed.atom3=atom2_CA;
				noed.atom4=atom2_C;
				NoeDihedral.push_back(noed);
			}
			if(atom2_N!=-1 && atom2_CA!=-1 && atom2_C!=-1 && atom3_N!=-1)
			{
				noed.atom1=atom2_N;
				noed.atom2=atom2_CA;
				noed.atom3=atom2_C;
				noed.atom4=atom3_N;
				NoeDihedral.push_back(noed);
			}
			if(atom2_C!=-1 && atom3_N!=-1 && atom3_CA!=-1 && atom3_C!=-1)
			{
				noed.atom1=atom2_C;
				noed.atom2=atom3_N;
				noed.atom3=atom3_CA;
				noed.atom4=atom3_C;
				NoeDihedral.push_back(noed);
			}
			if(atom3_N!=-1 && atom3_CA!=-1 && atom3_C!=-1 && atom4_N!=-1)
			{
				noed.atom1=atom3_N;
				noed.atom2=atom3_CA;
				noed.atom3=atom3_C;
				noed.atom4=atom4_N;
				NoeDihedral.push_back(noed);
			}
		}
	}
	if(Noe.size()) 
	dcd->getNoe(frame,Noe);
	if(NoeAngle.size())
	dcd->getNoeAngle(frame,NoeAngle);
	if(NoeDihedral.size())
	dcd->getNoeDihedral(frame,NoeDihedral);
	if(!frame.size()) return;
	int  icount_distance=0;
	int  icount_angle=0;
	int  icount_dihedral=0;

	
	//pass2
	for(int i=0;i<chain.size();i++)
	{
		vector<int> res_id;
		bturn_team->enumerateResId(chain[i],res_id);
		if(res_id.size()<4)  continue;
		for(int j=0;j<=res_id.size()-4;j++)
		{
			int id1=res_id[j];
			int id2=res_id[j+1];
			int id3=res_id[j+2];
			int id4=res_id[j+3];
			if(id2-id1!=1 && id3-id2!=1 && id4-id3!=1) continue;
		string filename="bturn_"+chain[i]+"_"+printNumber(id1)+"_"+printNumber(id4)+".dat";
		string statname="bturn_"+chain[i]+"_"+printNumber(id1)+"_"+printNumber(id4)+".stat";
		string histnamea="bturn_"+chain[i]+"_"+printNumber(id1)+"_"+printNumber(id4)+".hista";
		string histnamed="bturn_"+chain[i]+"_"+printNumber(id1)+"_"+printNumber(id4)+".histd";
		string histname2="bturn_"+chain[i]+"_"+printNumber(id1)+"_"+printNumber(id4)+".hist2";
			
			
			int atom1_CA=-1,atom2_CA=-1,atom3_CA=-1,atom4_CA=-1;
			int atom1_O=-1,atom4_HN=-1,atom4_N=-1;
			int atom2_C=-1,atom2_N=-1,atom1_C=-1;
			int atom3_N=-1,atom3_C=-1;
			for(int k=0;k<table.size();k++)
			{
				string mychain;
				if(bturn_team->find(table[k],mychain)==0) continue;
				if(table[k].res_id==id1 && table[k].atom_name=="CA") atom1_CA=table[k].atom_id;
				if(table[k].res_id==id2 && table[k].atom_name=="CA") atom2_CA=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="CA") atom3_CA=table[k].atom_id;
				if(table[k].res_id==id4 && table[k].atom_name=="CA") atom4_CA=table[k].atom_id;
				if(table[k].res_id==id4 && table[k].atom_name=="N")  atom4_N=table[k].atom_id;
				if(table[k].res_id==id4 && table[k].atom_name=="HN") atom4_HN=table[k].atom_id;
				if(table[k].res_id==id1 && table[k].atom_name=="O")  atom1_O=table[k].atom_id;	
				if(table[k].res_id==id2 && table[k].atom_name=="N")  atom2_N=table[k].atom_id;
				if(table[k].res_id==id2 && table[k].atom_name=="CA") atom2_CA=table[k].atom_id;
				if(table[k].res_id==id2 && table[k].atom_name=="C")  atom2_C=table[k].atom_id;
				if(table[k].res_id==id1 && table[k].atom_name=="C")  atom1_C=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="N")  atom3_N=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="CA") atom3_CA=table[k].atom_id;
				if(table[k].res_id==id3 && table[k].atom_name=="C")  atom3_C=table[k].atom_id;
			}
			if(atom1_CA==-1 || atom2_CA==-1 || atom3_CA==-1 || atom4_CA==-1) continue;
			//perform queries
			distance = Noe[icount_distance].D;icount_distance++;
			//dcd->printDistance(frame,distance,atom1_CA,atom4_CA);
			angle = NoeDihedral[icount_dihedral].D;icount_dihedral++;
			//dcd->printDAngle(frame,angle,atom1_CA,atom2_CA,atom3_CA,atom4_CA);
			int icount=0;
			for(int k=0;k<frame.size();k++)
			{
				icount+=(distance[k]<bturn_distance && fabs(angle[k])<bturn_angle)?1:0;
			}
			int toprint=1;
			if(icount*1.0/frame.size()<bturn_percent) toprint=0;
			//printf header
			FILE *fp;
			FILE *stat;
			if(toprint)
			{
			fp=fopen(filename.c_str(),"w");
			if(!fp) continue;
			string header1="#Distance statistics";
			string header2="#Angle Statistics";
			PrintStatFile(statname,"w",frame.size()/statcount,header1,frame,distance);
			PrintAngleStatFile(statname,"a",frame.size()/statcount,header2,frame,angle);
			string header="#"+printString("Frame",(FRAME_WIDTH-1))+" "+
				printFrameHeader("exist")+" "+
				printFrameHeader("type")+" "+printDistanceHeader("dist")+" "+
				printAngleHeader("tors")+" "+printDistanceHeader("hbd")+" "+
				printAngleHeader("hba")+" "+printAngleHeader("phi2")+" "+
				printAngleHeader("psi2")+" "+printAngleHeader("phi3")+" "+
				printAngleHeader("psi3");
			PrintLine(fp,header);
				
			}
			if(atom4_HN!=-1 && atom1_O!=-1)
			{
				hbd_distance = Noe[icount_distance].D;icount_distance++;
			}
			if(atom4_N!=-1 && atom1_O!=-1 && atom4_HN!=-1) 
			{
				hba_angle = NoeAngle[icount_angle].D;icount_angle++;
			}
			if(atom1_C!=-1 && atom2_N!=-1 && atom2_CA!=-1 && atom2_C!=-1)
			{
				phi2_angle=NoeDihedral[icount_dihedral].D;icount_dihedral++;
			}
			if(atom2_N!=-1 && atom2_CA!=-1 && atom2_C!=-1 && atom3_N!=-1)
			{
				psi2_angle = NoeDihedral[icount_dihedral].D;icount_dihedral++;
			}
			if(atom2_C!=-1 && atom3_N!=-1 && atom3_CA!=-1 && atom3_C!=-1)
			{
				phi3_angle=NoeDihedral[icount_dihedral].D;icount_dihedral++;
			}
			if(atom3_N!=-1 && atom3_CA!=-1 && atom3_C!=-1 && atom4_N!=-1)
			{
				psi3_angle = NoeDihedral[icount_dihedral].D;icount_dihedral++;
			}
			if(!toprint) continue;
			int existcount=0;
				int	ivcount=0;
				for(int l=0;l<Beta.size();l++) Beta[l].counter=0;
			IntVector2 CountTurn;
			CountTurn.resize(9);
			for(int k=0;k<9;k++) CountTurn[k].resize(frame.size());
			for(int k1=0;k1<9;k1++) 
				for(int k2=0;k2<frame.size();k2++) 
					CountTurn[k1][k2]=0;
			for(int k=0;k<frame.size();k++)
			{
				Data phi2=(atom1_C!=-1 && atom2_N!=-1 && atom2_CA!=-1 && atom2_C!=-1)?
					phi2_angle[k]:999.9;
				Data  psi2=(atom2_N!=-1 && atom2_CA!=-1 && atom2_C!=-1 && atom3_N!=-1)?
					psi2_angle[k]:999.9;
				Data phi3=(atom2_C!=-1 && atom3_N!=-1 && atom3_CA!=-1 && atom3_C!=-1)?
					phi3_angle[k]:999.9;
				Data psi3=(atom3_N!=-1 && atom3_CA!=-1 && atom3_C!=-1 && atom4_N!=-1)?
					psi3_angle[k]:999.9;
				Data print_distance=(atom4_HN!=-1 && atom1_O!=-1)?hbd_distance[k]:-1;

				int exist=(distance[k]<bturn_distance && fabs(angle[k])<bturn_angle)?1:0;
				existcount+=exist;
				//if(!exist) continue;
				Print(fp,printFrame(frame[k])+" ");
				/**/
				string betaType="null";
				
				for(int l=0;l<Beta.size();l++)
				{
					if(distance[k]<=bturn_distance)
					{
					if(phi2>=Beta[l].phi2-THETA1 && phi2<=Beta[l].phi2+THETA1
			                && psi2>=Beta[l].psi2-THETA1 && psi2<=Beta[l].psi2+THETA1
					&& phi3>=Beta[l].phi3-THETA1 && phi3<=Beta[l].phi3+THETA1
					&& psi3>=Beta[l].psi3-THETA1 && psi3<=Beta[l].psi3+THETA1) 
					{
						betaType=Beta[l].type;
						Beta[l].counter++;
						CountTurn[l][k]=1;
						break;
					}

					if(phi2>=Beta[l].phi2-THETA2 && phi2<=Beta[l].phi2+THETA2
			                && psi2>=Beta[l].psi2-THETA1 && psi2<=Beta[l].psi2+THETA1
					&& phi3>=Beta[l].phi3-THETA1 && phi3<=Beta[l].phi3+THETA1
					&& psi3>=Beta[l].psi3-THETA1 && psi3<=Beta[l].psi3+THETA1) 
					{
						betaType=Beta[l].type;
						Beta[l].counter++;
						CountTurn[l][k]=1;
						break;
					}

					if(phi2>=Beta[l].phi2-THETA1 && phi2<=Beta[l].phi2+THETA1
			                && psi2>=Beta[l].psi2-THETA2 && psi2<=Beta[l].psi2+THETA2
					&& phi3>=Beta[l].phi3-THETA1 && phi3<=Beta[l].phi3+THETA1
					&& psi3>=Beta[l].psi3-THETA1 && psi3<=Beta[l].psi3+THETA1) 
					{
						betaType=Beta[l].type;
						Beta[l].counter++;
						CountTurn[l][k]=1;
						break;
					}

					if(phi2>=Beta[l].phi2-THETA1 && phi2<=Beta[l].phi2+THETA1
			                && psi2>=Beta[l].psi2-THETA1 && psi2<=Beta[l].psi2+THETA1
					&& phi3>=Beta[l].phi3-THETA2 && phi3<=Beta[l].phi3+THETA2
					&& psi3>=Beta[l].psi3-THETA1 && psi3<=Beta[l].psi3+THETA1)
					{
						betaType=Beta[l].type;
						Beta[l].counter++;
						CountTurn[l][k]=1;
						break;
					}

					if(phi2>=Beta[l].phi2-THETA1 && phi2<=Beta[l].phi2+THETA1
			                && psi2>=Beta[l].psi2-THETA1 && psi2<=Beta[l].psi2+THETA1
					&& phi3>=Beta[l].phi3-THETA1 && phi3<=Beta[l].phi3+THETA1
					&& psi3>=Beta[l].psi3-THETA2 && psi3<=Beta[l].psi3+THETA2) 
					{
						betaType=Beta[l].type;
						Beta[l].counter++;
						CountTurn[l][k]=1;
						break;
					}
					}
				}
				if(distance[k]<=bturn_distance && betaType=="null")
				{
					betaType="IV";
					CountTurn[8][k]=1;
					ivcount++;
				}
				if(!exist) betaType="null";
				Print(fp,printFrame(exist)+" "+printFrameHeader(betaType)+" ");
				/**/
				Print(fp,printDistance(distance[k])+" "+printAngle(fabs(angle[k]))+" ");
				Data print_angle=(atom4_HN!=-1 && atom1_O!=-1 && atom4_N!=-1)?
					hba_angle[k]:999.9;
				Print(fp,printDistance(print_distance)+" "+printAngle(print_angle)+" ");
				PrintLine(fp,printAngle(phi2)+" "+printAngle(psi2)+" "+
				             printAngle(phi3)+" "+printAngle(psi3));
			}
			fclose(fp);
			stat=fopen(statname.c_str(),"a");
			PrintLine(stat,"#Bturn statistics");
			int dd=frame.size()/statcount;
			if(frame.size() % statcount!=0) dd++;
			PrintLine(stat,printString("#Begin block",FRAME_WIDTH)+" "+printFrame(dd));
			Print(stat,"#"+printString("Block",FRAME_WIDTH-1)+" "+
					printFrameHeader("Frame1")+" "+printFrameHeader("Frame2")+" ");
			IntVector2 res;
			res.resize(9);
			IntVector start;
			IntVector last;
			Print(stat,printFrameHeader("Bcount")+" ");
			for(int k1=0;k1<9;k1++)
			{
				SplitTable(frame.size()/statcount,start,last,frame,CountTurn[k1],res[k1]);
				if(k1!=8) Print(stat,printString(Beta[k1].type,FRAME_WIDTH)+" ");
				else	  Print(stat,printString("IV",FRAME_WIDTH)+" ");
			}
			PrintLine(stat,"");
			for(int k1=0;k1<start.size();k1++)
			{
				int ii=k1+1;
				int low=start[k1];
				int upper=last[k1];
				Print(stat,printFrame(ii)+" "+printFrame(low)+" "+printFrame(upper)+" ");
				int bcount=0;
				for(int k2=0;k2<9;k2++) bcount+=res[k2][k1];
				Print(stat,printFrame(bcount)+" ");
				for(int k2=0;k2<9;k2++)
					Print(fp,printFrame(res[k2][k1])+" ");
				PrintLine(stat,"");
			}
			PrintLine(stat,printString("#End block",FRAME_WIDTH)+" "+printFrame(dd));
			PrintLine(stat,"");
			PrintLine(stat,"#Complete trajectory");
			int low=frame[0];
			int upper=frame[frame.size()-1];
			Print(stat,printFrameHeader(" ")+" "+printFrame(low)+" "+printFrame(upper)+" ");
			int bcount=0;
			for(int k2=0;k2<8;k2++) bcount+=Beta[k2].counter;
			bcount+=ivcount;
				Print(stat,printFrame(bcount)+" ");
			for(int k1=0;k1<8;k1++)
				Print(stat,printFrame(Beta[k1].counter)+" ");
			PrintLine(stat,printFrame(ivcount));
			fclose(stat);
			Data dmin,dmax,davg,dstd;
			getVectorStatistics(distance,dmin,dmax,davg,dstd);
			Data amin,amax,aavg,astd;
			getAngleVectorStatistics(angle,amin,amax,aavg,astd);
			vector<HistStruct> st_dist;
			vector<HistStruct> st_angle;
			makeHist(distance,st_dist,dmin,dmax,bindist);
			makeHist(angle,st_angle,amin,amax,bindihe);

			fp=fopen(histnamed.c_str(),"w");
			for(int l=0;l<st_dist.size();l++)
			{
				PrintLine(fp,printDistance(st_dist[l].value)+" "+
					printPercent(st_dist[l].percent)+" "+
					printFrame(st_dist[l].count));
			}
			fclose(fp);

			fp=fopen(histnamea.c_str(),"w");
			for(int l=0;l<st_angle.size();l++)
			{
				PrintLine(fp,printAngle(st_angle[l].value)+" "+
					printPercent(st_angle[l].percent)+" "+
					printFrame(st_angle[l].count));
			}
			fclose(fp);

			fp=fopen(histname2.c_str(),"w");
			for(int k=0;k<st_dist.size();k++)
			{
				for(int m=0;m<st_angle.size();m++)
				{
					int icount =0;
					for(int l=0;l<distance.size();l++)
					{
						if(distance[l]>=st_dist[k].value-bindist/2.0 && 
						 distance[l]<=st_dist[k].value+bindist/2.0 && 
						 angle[l]>=st_angle[m].value-bindihe/2.0 &&
						 angle[l]<=st_angle[m].value+bindihe/2.0)
							icount++;
					}
					PrintLine(fp,printDistance(st_dist[k].value)+" "+
					 printAngle(st_angle[m].value)+" "+printFrame(icount));	
				}
				PrintLine(fp,"");
			}
			fclose(fp);
		}
	}
}

