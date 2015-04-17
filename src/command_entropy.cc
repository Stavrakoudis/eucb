# include <command_entropy.h>		
# include <pdb.h>
# include <global.h>

#define CONST1   (double)(48.50867485) /* (HBAR^2 / Kb) / (AMU*1e-20) */
#define AVOGADRO (double)(6.02214150e23)
#define KBOL     (double)(1.3806505e-23)
#define CONST2   (double)(0.020614869)/* (AMU*1e-20) / (HBAR^2 / Kb) */


static void sortArray(DoubleVector &x)
{
	for(int i=0;i<x.size();i++)
		for(int j=0;j<x.size()-1;j++)
		{
			Data t;
			if(x[j+1]<x[j]) 
			{
				t=x[j];
				x[j]=x[j+1];
				x[j+1]=t;
			}
		}
}


CommandEntropy::CommandEntropy(int n,Data t,int f,int fit,int pair)
	:Command("ENTROPY")
{
	nframes = n;
	temp =t;
	perresidue=f;
	fitflag=fit;
	pairflag=pair;
}

void	CommandEntropy::getEntropies(int eigendim,DoubleVector &evalue,Data &entropy1,Data &entropy2)
{
	sortArray(evalue);
	for(int j=0;j<eigendim/2;j++)
	{
		Data temp=evalue[j+1-1];
		evalue[j+1-1]=evalue[eigendim-j-1];
		evalue[eigendim-j-1]=temp;
	}
	Data entropy = 0.0;
	int ipos=1;
	for(ipos=1;ipos<=eigendim-6;ipos++)
	{
		if(fabs(evalue[ipos-1])<1e-10) break;// || evalue[ipos-1]<=0) break;
	}
	if(ipos==eigendim-5)
	{	
		  entropy = 0.0;
       		  for ( int j=1 ; j <= eigendim -6 ; j++ )
               	  	entropy += (sqrt(CONST1/(temp*fabs(evalue[j-1])))/expm1(sqrt(CONST1/(temp* fabs(evalue[j-1]) ) ))
                         - log( -expm1( -sqrt( CONST1 / ( temp * fabs(evalue[j-1] )))) ));
       		entropy *= (AVOGADRO * KBOL);
		entropy1=entropy;
	}
	entropy = 0.0;
       	for ( int j=1 ; j <= eigendim -6 ; j++ )
               entropy += log1p( CONST2 * temp * M_E * M_E * evalue[j-1] );
       	entropy *= (AVOGADRO * 0.50 * KBOL );
	entropy2=entropy;
}

void	CommandEntropy::getSubEntropy(IntVector &posAtom,IntVector &startFrame,
				DoubleVector &entropy1,DoubleVector &entropy2)
{
	DoubleVector3 data;
	DoubleVector2 evector;
	DoubleVector  evalue;
	IntVector F;
	startFrame.resize(0);
	dcd->getEntropyData(F,data,posAtom,nframes,fitflag,pdbpos);
	int count =data.size();
	int dd=nframes-1;
	for(int i=0;i<count;i++)
	{
		int eigendim = 3*posAtom.size();
		Data e1,e2;
		if(!eigen(3*posAtom.size(),3*posAtom.size(),data[i],evector,evalue))
		{
			e1=0.0;
			e2=0.0;
		}
		else
		getEntropies(eigendim,evalue,e1,e2);
		if(isnan(e1) || isinf(e1)) e1=0.0;
		if(isnan(e2) || isinf(e2)) e2=0.0;
		entropy1.push_back(e1);
		entropy2.push_back(e2);
		int a=dd<F.size()?F[dd]:F[F.size()-1];
		startFrame.push_back(a);
		dd+=nframes;
	}
}


void CommandEntropy::Run()
{
	DoubleVector3 data;
	DoubleVector2 evector;
	DoubleVector  evalue;
	IntVector F;
	IntVector posAtom;
	string chain;
	for(int i=0;i<table.size();i++) 
	{
		if(!team || team->find(table[i],chain)) 
			posAtom.push_back(i);
	}

	if(pairflag && team)
	{
		FILE *fp=fopen("entropy.dat","w");
		PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
			printString("Res1",5)+" "+
			printString("Res2",5)+" "+
			printDistanceHeader("Schlitter")+" "+
			printDistanceHeader("Andricioaei"));
		IntVector posAtom1,posAtom2,posAtom3;
		vector<string> chainlist;
		team->enumerateChains(chainlist);
		for(int i=0;i<chainlist.size();i++)
		{
			vector<int> reslist;
			team->enumerateResId(chainlist[i],reslist);
			for(int j=0;j<reslist.size();j++)
			{
				posAtom1.resize(0);
				posAtom2.resize(0);
				posAtom3.resize(0);
				for(int l=0;l<posAtom.size();l++)
					if(table[posAtom[l]].chain_id==chainlist[i]
					 && table[posAtom[l]].res_id==reslist[j])
					{
						posAtom1.push_back(l);
					}
				for(int k=0;k<j;k++)
				{
					posAtom2.resize(0);
					posAtom3.resize(0);
					for(int l=0;l<posAtom.size();l++)
					{
						if(table[posAtom[l]].chain_id==chainlist[i]
						&& table[posAtom[l]].res_id==reslist[k])
						{
							posAtom2.push_back(l);
						}
					}
					for(int l=0;l<posAtom1.size();l++)
						posAtom3.push_back(posAtom1[l]);
					for(int l=0;l<posAtom2.size();l++)
						posAtom3.push_back(posAtom2[l]);
				if(posAtom1.size()==0 || posAtom2.size()==0 || posAtom3.size()==0) 
					continue;
				IntVector startFrame;
				DoubleVector entropy1_1,entropy1_2;
				DoubleVector entropy2_1,entropy2_2;
				DoubleVector entropy3_1,entropy3_2;
				getSubEntropy(posAtom1,startFrame,entropy1_1,entropy1_2);
				getSubEntropy(posAtom2,startFrame,entropy2_1,entropy2_2);
				getSubEntropy(posAtom3,startFrame,entropy3_1,entropy3_2);
				for(int l=0;l<startFrame.size();l++)
				{
					Data aa=entropy3_1[l]-entropy1_1[l]-entropy2_1[l];
					Data bb=entropy3_2[l]-entropy1_2[l]-entropy2_2[l];
					PrintLine(fp,printFrame(startFrame[l])+" "+
						chainlist[i]+printNumber(reslist[j],5)+" "+
						chainlist[i]+printNumber(reslist[k],5)+" "+
						printDistance(aa)+" "+printDistance(bb));	
				}
				}
			}
		}
		fclose(fp);
		return;
	}

	FILE *fp=fopen("entropy.dat","w");
	if(perresidue==0)
	{
		dcd->getEntropyData(F,data,posAtom,nframes,fitflag,pdbpos);
		int count =data.size();
		PrintLine(fp,"#"+printString("Frame",FRAME_WIDTH-1)+" "+
			printDistanceHeader("Schlitter")+" "+
			printDistanceHeader("Andricioaei"));
			int dd=nframes-1;
		for(int i=0;i<count;i++)
		{
			int eigendim = 3*posAtom.size();
			eigen(3*posAtom.size(),3*posAtom.size(),data[i],evector,evalue);
			Data entropy1,entropy2;
			getEntropies(eigendim,evalue,entropy1,entropy2);
			if(isnan(entropy1) || isinf(entropy1)) entropy1=0.0;
			if(isnan(entropy2) || isinf(entropy2)) entropy2=0.0;
			PrintLine(fp,printFrame(dd<F.size()?F[dd]:F[F.size()-1])+" "+
			printDistance(entropy1)+" "+printDistance(entropy2));	
			dd+=nframes;
		}
	}
	else
	if(team)
	{


		int oldLast=dcd->getFrameEnd();
		dcd->setFrameEnd(nframes);
		dcd->getEntropyData(F,data,posAtom,nframes,fitflag,pdbpos);
		FILE *fnorm=fopen("entropy_norm.dat","w");

		PrintLine(fp,"#"+printString("ChainId",8)+" "+printString("ResId",8)+" "+
			printString("Sch(all)",9)+" "+
			printString("And(all)",9)+" "+
			printString("Sch(noh)",9)+" "+
			printString("And(noh)",9)+" "+
			printString("Sch(b)",9)+" "+
			printString("And(b)",9)+" "+
			printString("Sch(s)",9)+" "+
			printString("And(s)",9)
			);

		PrintLine(fnorm,"#"+printString("ChainId",8)+" "+printString("ResId",8)+" "+
			printString("Sch(all)",9)+" "+
			printString("And(all)",9)+" "+
			printString("Sch(noh)",9)+" "+
			printString("And(noh)",9)+" "+
			printString("Sch(b)",9)+" "+
			printString("And(b)",9)+" "+
			printString("Sch(s)",9)+" "+
			printString("And(s)",9)
			);

		for(int i=0;i<posAtom.size();i++)
		{
			DoubleVector2 Pie;
			int resid=table[posAtom[i]].res_id;
			string chain_id=table[posAtom[i]].chain_id;
			int startI=i;
			int nbackbone=0;
			int nsidechain=0;
			int noh=0;
			while(i<posAtom.size() && table[posAtom[i]].res_id==resid
				&& table[posAtom[i]].chain_id==chain_id) 
			{
				if(isBackbone(table[posAtom[i]])) nbackbone++;
				if(isSidechain(table[posAtom[i]])) nsidechain++;
				if(table[posAtom[i]].atom_name[0]!='H') noh++;
				i++;
			}
			int endI=i-1;
			int diff=(endI-startI+1);
			int eigendim=3*diff;
			Pie.resize(eigendim);
			for(int j=0;j<eigendim;j++) Pie[j].resize(eigendim);
			int ipos=0;
			for(int j=3*startI;j<=3*endI+2;j++)
			{
				int jpos=0;
				for(int k=3*startI;k<=3*endI+2;k++)
				{		
					Pie[ipos][jpos]=data[0][j][k];
					jpos++;
				}
				ipos++;
			}
			DoubleVector2 evector;
			DoubleVector  evalue;
			Data entropy1,entropy2;

			eigen(eigendim,eigendim,Pie,evector,evalue);
			getEntropies(eigendim,evalue,entropy1,entropy2);
			if(isnan(entropy1) || isinf(entropy1)) entropy1=0.0;
			if(isnan(entropy2) || isinf(entropy2)) entropy2=0.0;

			Print(fp,printString(chain_id,8)+" "+
				printNumber(resid,8)+" "+
				printNumber(entropy1,9,2)+" "+
				printNumber(entropy2,9,2)+" ");

			Print(fnorm,printString(chain_id,8)+" "+
				printNumber(resid,8)+" "+
				printNumber(entropy1/diff,9,2)+" "+
				printNumber(entropy2/diff,9,2)+" ");

			Pie.resize(0);
			eigendim=3*noh;
			Pie.resize(eigendim);
			for(int j=0;j<eigendim;j++) Pie[j].resize(eigendim);
			ipos=0;
			for(int j=3*startI;j<=3*endI+2;j++)
			{
				int jpos=0;
				if(table[posAtom[j/3]].atom_name[0]!='H')
				{
					for(int k=3*startI;k<=3*endI+2;k++)
					{
						if(table[posAtom[k/3]].atom_name[0]!='H')
						{
							Pie[ipos][jpos]=data[0][j][k];
							jpos++;
						}
					}
					ipos++;
				}
			}
			evector.resize(0);evalue.resize(0);
			eigen(eigendim,eigendim,Pie,evector,evalue);
			getEntropies(eigendim,evalue,entropy1,entropy2);
			if(isnan(entropy1) || isinf(entropy1)) entropy1=0.0;
			if(isnan(entropy2) || isinf(entropy2)) entropy2=0.0;
			Print(fp,printNumber(entropy1,9,2)+" "+printNumber(entropy2,9,2)+" ");
			Print(fnorm,
				printNumber(entropy1/noh,9,2)+" "+printNumber(entropy2/noh,9,2)+" ");

			Pie.resize(0);
			eigendim=3*nbackbone;
			Pie.resize(eigendim);
			for(int j=0;j<eigendim;j++) Pie[j].resize(eigendim);
			ipos=0;
			for(int j=3*startI;j<=3*endI+2;j++)
			{
				int jpos=0;
				if(isBackbone(table[posAtom[j/3]]))
				{
					for(int k=3*startI;k<=3*endI+2;k++)
					{
						if(isBackbone(table[posAtom[k/3]]))
						{
							Pie[ipos][jpos]=data[0][j][k];
							jpos++;
						}
					}
					ipos++;
				}
			}
			evector.resize(0);evalue.resize(0);
			eigen(eigendim,eigendim,Pie,evector,evalue);
			getEntropies(eigendim,evalue,entropy1,entropy2);
			if(isnan(entropy1) || isinf(entropy1)) entropy1=0.0;
			if(isnan(entropy2) || isinf(entropy2)) entropy2=0.0;
			Print(fp,printNumber(entropy1,9,2)+" "+printNumber(entropy2,9,2)+" ");
			Print(fnorm,
				printNumber(entropy1/nbackbone,9,2)+" "+printNumber(entropy2/nbackbone,9,2)+" ");


			if(nsidechain==0) 	
			{
				PrintLine(fp,printNumber(0.0,9,2)+" "+printNumber(0.0,9,2));
				PrintLine(fnorm,printNumber(0.0,9,2)+" "+printNumber(0.0,9,2));
				i--;
				continue;
			}
			Pie.resize(0);
			eigendim=3*nsidechain;
			Pie.resize(eigendim);
			for(int j=0;j<eigendim;j++) Pie[j].resize(eigendim);
			ipos=0;
			for(int j=3*startI;j<=3*endI+2;j++)
			{
				int jpos=0;
				if(isSidechain(table[posAtom[j/3]]))
				{
					for(int k=3*startI;k<=3*endI+2;k++)
					{
						if(isSidechain(table[posAtom[k/3]]))
						{
							Pie[ipos][jpos]=data[0][j][k];
							jpos++;
						}
					}
					ipos++;
				}
			}
			evector.resize(0);evalue.resize(0);
			eigen(eigendim,eigendim,Pie,evector,evalue);
			getEntropies(eigendim,evalue,entropy1,entropy2);
			if(isnan(entropy1) || isinf(entropy1)) entropy1=0.0;
			if(isnan(entropy2) || isinf(entropy2)) entropy2=0.0;
			PrintLine(fp,printNumber(entropy1,9,2)+" "+printNumber(entropy2,9,2));
			PrintLine(fnorm,
				printNumber(entropy1/nsidechain,9,2)+" "+printNumber(entropy2/nsidechain,9,2)+" ");
			i--;
		}
	
		dcd->setFrameEnd(oldLast);
		fclose(fnorm);
	}
	fclose(fp);
}


CommandEntropy::~CommandEntropy()
{
}
