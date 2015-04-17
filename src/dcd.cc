# include <dcd.h>
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include <math.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <fcntl.h>
# include <limits.h>
# include <psf.h>
# include <parfile.h>

#define NR_END	1

static int imin(int a,int b)
{
	return a<b?a:b;
}

static int imax(int a,int b)
{
	return a>b?a:b;
}

static float *make_vector(long nl,long nh)
{
	float *v;
	v=(float*) malloc((size_t)((nh-nl+1+NR_END)*sizeof(float)));
	return v-nl+NR_END;
}

/*	Skip a frame from the dcd file.
 * */
int	Dcd::skipFrame(int &frame)
{
	 lseek(dcd,(off_t)(bytes_per_set*(frameStep-1)),SEEK_CUR);
	 frame+=frameStep-1;
	 if(frame>frameEnd || frame>headerframes)  return 0;
	 return 1;
}

Point	Dcd::makePoint(int index)
{
	Point pt;
	pt.x=CAs[index][0];
	pt.y=CAs[index][1];
	pt.z=CAs[index][2];
	return pt;
}

Point	Dcd::makeCenter(IntVector &list,int wflag)
{
	Point center;
	center.x=0.0;
	center.y=0.0;
	center.z=0.0;
	Data sum_weight=0.0;	
	extern vector<atom> table;
	for(int i=0;i<list.size();i++)
	{
		Data w=1.0;
		atom a=table[list[i]-1];
		if(list.size()>1 && wflag) w=a.dnum2;
		sum_weight+=w;
		center.x+=w*CAs[list[i]-1][0];
		center.y+=w*CAs[list[i]-1][1];
		center.z+=w*CAs[list[i]-1][2];
	}
	center.x/=sum_weight;
	center.y/=sum_weight;
	center.z/=sum_weight;
	return center;
}

/*	Read the atoms' position in the current frame.
 * */
int	Dcd::readCAs()
{
 	unsigned int read_bytes=(unsigned int)(read( dcd, (void *)(&dcd_frame[1]), (size_t)(bytes_per_set) ));
 	if(read_bytes!= bytes_per_set ) return 0;
    	build_CAs( (int)(natoms) + 2);
	return 1;
}

Data	Dcd::getDistance(Point a,Point b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

Data	Dcd::getDistance(int DISTFIRST,int DISTSECOND)
{
	return getDistance(makePoint(DISTFIRST),makePoint(DISTSECOND));
}

Data	Dcd::getAngle(int DISTFIRST,int DISTSECOND,int DISTTHIRD)
{
	return getAngle(makePoint(DISTFIRST),makePoint(DISTSECOND),makePoint(DISTTHIRD));
}

Data	Dcd::getAngle(Point p1,Point p2,Point p3)
{
	Data xab = p1.x-p2.x;
       	Data yab = p1.y-p2.y;
       	Data zab = p1.z-p2.z;;
       	Data xcb = p3.x-p2.x;
       	Data ycb = p3.y-p2.y;
       	Data zcb = p3.z-p2.z;
       	Data rab2 = xab*xab + yab*yab + zab*zab;
       	Data rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
       	Data  rabc = sqrt(rab2 * rcb2);
	if(rabc!=0) 
	{
           		Data cosine = (xab*xcb + yab*ycb + zab*zcb) / rabc;
          	 	cosine = dmin(1.0,dmax(-1.0,cosine));
           		return  180.0/M_PI * acos(cosine);
	}
	return 0.0;
}

static int filesize(char *filename)
{
	struct stat filestatus;
	stat( filename, &filestatus );
	return  filestatus.st_size;
}

Dcd::Dcd(char *s)
{

	strcpy(filename, s);
	have_cell = 0;
	cell_offset = 0;
	headerframes = 0;
	natoms = 0;
	dcd=open(filename,O_RDONLY);
	int bytes_so_far=0;
	if(dcd == -1) 
	{
		psf_error("Can not open DCD file for read. Abort.");
	}
	int k=read(dcd,dcd_head,(size_t)(92));
	bytes_so_far+=92;
	if(k!=92)
	{
		close(dcd);
		psf_error("Premature EOF for DCD file ? Abort.");
	}
       if ( strncasecmp( (char *)(&dcd_head[4]), "CORD", 4) != 0 )
       {
	       close(dcd);
	       psf_error("Premature EOF for DCD file ? Abort.");
       }
       if ( (unsigned int)(*((unsigned int *)(&dcd_head[48]))) == 1 )
       {
	       have_cell=1;
	       cell_offset=14;
       }
       if ( (unsigned int)(*((unsigned int *)(&dcd_head[52]))) == 1 )
       {
	       close(dcd);
	       psf_error("It appears that this is a CHARMM-specific DCD file. Abort.");
       }
       headerframes = (unsigned int)(*((unsigned int *)(&dcd_head[8])));
       if((int)(*((int *)(&dcd_head[8]))) < 0 ||
	 (int)(*((int *)(&dcd_head[12]))) < 0 ||
	 (int)(*((int *)(&dcd_head[16]))) < 0 )
       {
	       close(dcd);
	       psf_error("Negative number of frames ? Wrong endian for DCD file ?");
       }
       if(read( dcd, &auint, (size_t)(4) ) != 4 )
       {
		  close(dcd);
                  psf_error("Premature EOF for DCD file ? Abort.");
       }
       if(read( dcd, &auint, (size_t)(4) ) != 4 )
       {
		 close(dcd);
                 psf_error("Premature EOF for DCD file ? Abort.");
       }
	bytes_so_far+=8;
       for(int i=0;i<auint;i++)
       {
	      if ( read( dcd, &title, (size_t)(80) ) != 80 )
	      {
		    close(dcd);
	            psf_error("Premature EOF for DCD file ? Abort.");
	      }
	bytes_so_far+=80;
	      title[80] = 0;
       }
       if(read( dcd, dcd_head2, (size_t)(16) ) != 16 )
       {
	   close(dcd);
           psf_error("Premature EOF for DCD file ? Abort.");
       }
	bytes_so_far+=16;
       natoms = (unsigned int)(*((unsigned int *)(&dcd_head2[8])));
       if(have_cell)
       {
	       bytes_per_set = (unsigned int)(3 * ( 4 * natoms + 8 ) + 56);
	       wbytes_per_set = (unsigned int)(3 * ( 4 *natoms + 8 )  + 56);
       }
       else
       {
	        bytes_per_set = (unsigned int)(3 * ( 4 * natoms + 8 ));
		wbytes_per_set = (unsigned int)(3 * ( 4 * natoms + 8 ));

       }
       dcd_frame = make_vector( 1, bytes_per_set / 4 );
	//GIANNIS
	//headerframes=(filesize(filename)-bytes_so_far)/bytes_per_set; // here cancel this line, what was for?
	//END OF GIANNIS
       wdcd_frame = make_vector( 1, wbytes_per_set / 4 );
       if(have_cell)
       {
	     *((unsigned int *)(&wdcd_frame[1])) = (unsigned int)(48);
	     *((unsigned int *)(&wdcd_frame[14])) = (unsigned int)(48);
             *((unsigned int *)(&wdcd_frame[14 + 1])) = (unsigned int)( natoms*4);
	     *((unsigned int *)(&wdcd_frame[14 + 2+natoms])) = (unsigned int)( natoms*4);
       	     *((unsigned int *)(&wdcd_frame[14 + 3+natoms])) = (unsigned int)( natoms*4);
	     *((unsigned int *)(&wdcd_frame[14 + 4+2*natoms])) = (unsigned int)( natoms*4);
	     *((unsigned int *)(&wdcd_frame[14 + 5+2*natoms])) = (unsigned int)( natoms*4);
	     *((unsigned int *)(&wdcd_frame[14 + 6+3*natoms])) = (unsigned int)( natoms*4);
       }
       else
       {
        *((unsigned int *)(&wdcd_frame[1])) = (unsigned int)( natoms*4);
        *((unsigned int *)(&wdcd_frame[2+natoms])) = (unsigned int)( natoms*4);
        *((unsigned int *)(&wdcd_frame[3+natoms])) = (unsigned int)( natoms*4);
        *((unsigned int *)(&wdcd_frame[4+2*natoms])) = (unsigned int)( natoms*4);
        *((unsigned int *)(&wdcd_frame[5+2*natoms])) = (unsigned int)( natoms*4);
        *((unsigned int *)(&wdcd_frame[6+3*natoms])) = (unsigned int)( natoms*4);
       }
       posCAs.resize(natoms);
       CAs.resize(natoms);
       for(int i=0;i<natoms;i++)
       {
	       posCAs[i]=i;
	       CAs[i].resize(3);
       }
       frameStart = 1;
       frameEnd=headerframes;//INT_MAX;
       frameStep=1;
}

Data	Dcd::dihedral(int x1,int  x2,int  x3,int x4)
{
	return torsion(CAs[x1][0],CAs[x1][1],CAs[x1][2],
			CAs[x2][0],CAs[x2][1],CAs[x2][2],
			CAs[x3][0],CAs[x3][1],CAs[x3][2],
			CAs[x4][0],CAs[x4][1],CAs[x4][2]);
}

Data	Dcd::dihedral(DoubleVector x1,DoubleVector x2,DoubleVector x3,DoubleVector x4)
{
	return torsion(x1[0],x1[1],x1[2],
			x2[0],x2[1],x2[2],
			x3[0],x3[1],x3[2],
			x4[0],x4[1],x4[2]);
}

Data	Dcd::dihedral(Point x1,Point x2,Point x3,Point x4)
{
	return torsion(x1.x,x1.y,x1.z,
			x2.x,x2.y,x2.z,
			x3.x,x3.y,x3.z,
			x4.x,x4.y,x4.z);
}


void	Dcd::setFrameStart(int s)
{
	if(s>0 && s<=headerframes ) frameStart = s;
	else
		psf_error("FRAME START OUT OF LIMIT");
}

void	Dcd::setFrameEnd(int s)
{
	if(s>0) frameEnd=imin(s,headerframes);
	else
		psf_error("FRAME END OUT OF LIMIT");
}

void	Dcd::setFrameStep(int s)
{
	if(s>0 && s<=headerframes) frameStep=s;
	else
		psf_error("FRAME STEP OUT OF LIMIT");
}

const 	int Dcd::getFrameStart() const
{
	return frameStart;
}

const	int Dcd::getFrameEnd() const
{
	return frameEnd;
}


const	int Dcd::getFrameStep() const
{
	return frameStep;
}


void	Dcd::build_CAs(int offset)
{
	int i,atomno;
	if(have_cell)
	{
		cell_a     = (Data)(*((Data *)(&dcd_frame[2] +  0)));
                cell_gamma = (Data)(*((Data *)(&dcd_frame[2] +  2)));
                cell_b     = (Data)(*((Data *)(&dcd_frame[2] +  4)));
                cell_beta  = (Data)(*((Data *)(&dcd_frame[2] +  6)));
                cell_alpha = (Data)(*((Data *)(&dcd_frame[2] +  8)));
                cell_c     = (Data)(*((Data *)(&dcd_frame[2] + 10)));
		if ( cell_alpha <= 1.0 && cell_beta <= 1.0 && cell_gamma <= 1.0)
		 {
	              if ( cell_alpha == 0.0 )
	                    cell_alpha = 90.00;
	              else
                      cell_alpha = 180.0 * acos( cell_alpha ) / M_PI;
		      if ( cell_beta == 0.0 )
		                 cell_beta = 90.00;
		       else
		       cell_beta = 180.0 * acos( cell_beta ) / M_PI;
		     if ( cell_gamma == 0.0 )
		      cell_gamma = 90.00;
		    else
		    cell_gamma = 180.0 * acos( cell_gamma ) / M_PI;
		   }
	}
	
	for ( i=0 ; i < natoms ; i++ )
	{
	      const int CELL_OFFSET = cell_offset;
	      atomno = posCAs[i];
	      CAs[i][0] = dcd_frame[ CELL_OFFSET + 2+atomno];
	      CAs[i][1] = dcd_frame[ CELL_OFFSET + 2+atomno+offset];
	      CAs[i][2] = dcd_frame[ CELL_OFFSET + 2+atomno+2*offset];
	}
}

const	unsigned int Dcd::getnatoms() const
{
	return natoms;
}

const	int Dcd::getframes() const
{
	return headerframes;
}

void	Dcd::printRMSF(vector<int> &res_number,vector<Data> &rmsf)
{
	vector<Data> xmean;
	vector<Data> ymean;
	vector<Data> zmean;
	vector<Data> variance;

	variance.resize(natoms);
	xmean.resize(natoms);
	ymean.resize(natoms);
	zmean.resize(natoms);
	for(int i=0;i<natoms;i++)
	{
		variance[i]=xmean[i]=ymean[i]=zmean[i]=0;
	}

	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;

	/**/

	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<res_number.size();i++)
		{
			int index=res_number[i]-1;
			xmean[index]+=CAs[index][0];
			ymean[index]+=CAs[index][1];
			zmean[index]+=CAs[index][2];
		}
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	for(int i=0;i<xmean.size();i++)
	{
		xmean[i]/=icount;
		ymean[i]/=icount;
		zmean[i]/=icount;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
	FIRST = frameStart;
	LAST  = frameEnd;
	istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	frame= FIRST;
	firstframe = FIRST;
	STEP=frameStep;
	itotalframes=0;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<res_number.size();i++)
		{
			int index=res_number[i]-1;
			variance[index]+=sqrt(
			pow(CAs[index][0]-xmean[index],2.0)+
			pow(CAs[index][1]-ymean[index],2.0)+
			pow(CAs[index][2]-zmean[index],2.0));
		}
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
	for(int i=0;i<res_number.size();i++) rmsf[i]=variance[res_number[i]-1]/icount;
}

int	Dcd::nframes()
{
	int frame= frameStart;
	int icount=1;
	printf ("End=%d, Header=%d\n", frameEnd, headerframes); //here

	while(frame<=frameEnd)
	{
	 if((frame-frameStart)%frameStep==0)
	 {
	 	frame+=frameStep-1;
		frame++;
	 	if(frame>frameEnd || frame>headerframes)  break;
		icount++;
	 }
	}
	return icount;
}

void	Dcd::getStackAverages(IntVector &F,vector<stack_table> &stable,vector<StackAtom> 
		&stack_atom)
{
	int icount =0;
	F.resize(nframes());
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;
	vector<IntVector> stack_flag;
	stack_flag.resize(stable.size());
	for(int i=0;i<stack_flag.size();i++) 
	{
		stack_flag[i].resize(stable.size());	
		for(int j=0;j<stack_flag[i].size();j++) stack_flag[i][j]=0;
	}
	/**/

	Point x1,y1,z1;
	Point x2,y2,z2;
	Point avg1;
	Point avg2;
	extern vector<atom> table;

	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
		IntVector list;
		list.resize(3);
		for(int i=0;i<stable.size();i++)
		{
			list[0]=stable[i].id1-1;list[1]=stable[i].id2-1;list[2]=stable[i].id3-1;
			avg1=makeCenter(list,0);
			for(int j=0;j<i;j++)
			{
				list[0]=stable[j].id1-1;list[1]=stable[j].id2-1;list[2]=stable[j].id3-1;
				avg2=makeCenter(list,0);
				Data dist=getPointDistance(avg1,avg2);
				x1.x=CAs[stable[i].id1-1][0];x1.y=CAs[stable[i].id2-1][0];x1.z=CAs[stable[i].id3-1][0];
				y1.x=CAs[stable[i].id1-1][1];y1.y=CAs[stable[i].id2-1][1];y1.z=CAs[stable[i].id3-1][1];
				z1.x=CAs[stable[i].id1-1][2];z1.y=CAs[stable[i].id2-1][2];z1.z=CAs[stable[i].id3-1][2];
				x2.x=CAs[stable[j].id1-1][0];x2.y=CAs[stable[j].id2-1][0];x2.z=CAs[stable[j].id3-1][0];
				y2.x=CAs[stable[j].id1-1][1];y2.y=CAs[stable[j].id2-1][1];y2.z=CAs[stable[j].id3-1][1];
				z2.x=CAs[stable[j].id1-1][2];z2.y=CAs[stable[j].id2-1][2];z2.z=CAs[stable[j].id3-1][2];
				Data angle=getPlaneAngle(x1,y1,z1,x2,y2,z2);
				int ipos=-1;
				for(int k=0;k<stack_atom.size();k++)
				{
					if(
					 stack_atom[k].atom1_id1==stable[i].id1  && 
					 stack_atom[k].atom1_id2==stable[i].id2 && 
					 stack_atom[k].atom1_id3==stable[i].id3 &&
					 stack_atom[k].atom2_id1==stable[j].id1 &&
					 stack_atom[k].atom2_id2==stable[j].id2 && 
					 stack_atom[k].atom2_id3==stable[j].id3)
					 {ipos=k;break;}
				}
				if(ipos==-1)
				{
					int s=stack_atom.size();
					stack_atom.resize(s+1);	
					stack_atom[s].D.resize(F.size());
					stack_atom[s].A.resize(F.size());
					stack_atom[s].atom1_id1=stable[i].id1;
					stack_atom[s].atom1_id2=stable[i].id2; 
					stack_atom[s].atom1_id3=stable[i].id3;
					stack_atom[s].atom2_id1=stable[j].id1;
					stack_atom[s].atom2_id2=stable[j].id2; 
					stack_atom[s].atom2_id3=stable[j].id3;
					ipos=s;
				}
				stack_atom[ipos].D[icount]=dist;
				stack_atom[ipos].A[icount]=angle;
			}
		}
		F[icount]=frame;
 		icount++;
 		frame++;
 	}
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::acceptors(IntVector &f,vector<AcceptorStruct> &st,IntVector &Donor)
{
	f.resize(nframes());
	for(int i=0;i<st.size();i++)
	{
		st[i].dist.resize(st[i].Donor.size());
		st[i].angle.resize(st[i].Donor.size());
		st[i].HydroFlag.resize(f.size());
		st[i].HydroDonor.resize(st[i].Donor.size());
		for(int j=0;j<st[i].dist.size();j++) st[i].dist[j].resize(f.size());
		for(int j=0;j<st[i].angle.size();j++) st[i].angle[j].resize(f.size());
		for(int j=0;j<st[i].HydroDonor.size();j++) 
		{
			st[i].HydroDonor[j].resize(st[i].Hydro[j].size());
			for(int k=0;k<st[i].HydroDonor[j].size();k++)
				st[i].HydroDonor[j][k].resize(f.size());
		}
	}
	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
	 	Data d,d1;
	 	int DISTFIRST,DISTSECOND,DISTTHIRD,DISTFOURTH;
		for(int i=0;i<st.size();i++)
		{
			DISTTHIRD=st[i].acceptor_id-1;
			for(int j=0;j<st[i].Donor.size();j++)
			{
				DISTFIRST=Donor[st[i].Donor[j]]-1;
				Data d=getDistance(DISTFIRST,DISTTHIRD);
				st[i].dist[j][icount]=d;
				st[i].angle[j][icount]=-1e+100;
				st[i].HydroFlag[icount]=0;
				for(int k=0;k<st[i].Hydro[j].size();k++)
				{
					DISTSECOND=st[i].Hydro[j][k]-1;
					Data dd=getAngle(DISTFIRST,DISTSECOND,DISTTHIRD);
					if(dd>=st[i].angle[j][icount])
					{
						st[i].angle[j][icount]=dd;
						st[i].HydroFlag[icount]=st[i].Hydro[j][k];
						st[i].HydroDonor[j][k][icount]=getDistance(DISTSECOND,DISTTHIRD);

					}
				}
			}
		}
 	f[icount]=frame;
 	icount++;
 	frame++;
	}
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::hbondsCheck(Data cdist,Data cangle,Data cpercent,IntVector &X1,
		         vector<IntVector> &X2,IntVector &X3,vector<IntVector> &Flag)
{
	Flag.resize(X1.size());
	for(int i=0;i<Flag.size();i++)
	{
		 Flag[i].resize(X3.size());
		for(int j=0;j<Flag[i].size();j++) Flag[i][j]=0;
	}
	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
	 	Data d,d1;
	 	int DISTFIRST,DISTSECOND,DISTTHIRD,DISTFOURTH;
		for(int i=0;i<X1.size();i++)
		{
			for(int j=0;j<X3.size();j++)
			{
				DISTFIRST=X1[i]-1;
				DISTTHIRD=X3[j]-1;
				Data d=getDistance(DISTFIRST,DISTTHIRD);
				if(d<=cdist)
				{
					for(int k=0;k<X2[i].size();k++)
					{
						DISTSECOND=X2[i][k]-1;
						Data ang=getAngle(DISTFIRST,DISTSECOND,DISTTHIRD);
						if(ang>=cangle)
						{
							Flag[i][j]++;
							break;
						}
					}
				}
			}
		}		
 	icount++;
 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	for(int i=0;i<Flag.size();i++)
		for(int j=0;j<Flag[i].size();j++)
		{
			if(Flag[i][j]>=1) Flag[i][j]=1;
			//if(Flag[i][j]*1.0/icount>=cpercent) Flag[i][j]=1;
			else			          Flag[i][j]=0;
		}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::hbonds(IntVector &f,vector<HbondStruct> &st)
{
	f.resize(nframes());
	printf ( "f-size : %d\n", f.size() );  // here
	printf ( "st-size: %d\n", st.size() ); // here
	for(int i=0;i<st.size();i++)
	{
		st[i].dist.resize(st[i].Acceptor.size());
		st[i].angle.resize(st[i].Acceptor.size());
		st[i].HydroFlag.resize(f.size());
		st[i].HydroAcc.resize(st[i].Acceptor.size());
		for(int j=0;j<st[i].dist.size();j++) {
			st[i].dist[j].resize(f.size()); 
			//if (i%100==0 && j%100==0) printf ( "i=,%d j= %d\n", i,j );  // here
		}
		for(int j=0;j<st[i].angle.size();j++) st[i].angle[j].resize(f.size());
		for(int j=0;j<st[i].HydroAcc.size();j++) st[i].HydroAcc[j].resize(f.size());
	}
	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
	 	Data d,d1;
	 	int DISTFIRST,DISTSECOND,DISTTHIRD,DISTFOURTH;
		for(int i=0;i<st.size();i++)
		{
			DISTFIRST=st[i].donor_id-1;
			for(int j=0;j<st[i].Acceptor.size();j++)
			{
				DISTTHIRD=st[i].Acceptor[j]-1;
				Data d=getDistance(DISTFIRST,DISTTHIRD);
				st[i].dist[j][icount]=d;
				st[i].angle[j][icount]=-1e+100;
				st[i].HydroFlag[icount]=0;
				for(int k=0;k<st[i].Hydro.size();k++)
				{
					DISTSECOND=st[i].Hydro[k]-1;
					Data dd=getAngle(DISTFIRST,DISTSECOND,DISTTHIRD);
					if(dd>=st[i].angle[j][icount])
					{
						st[i].angle[j][icount]=dd;
						st[i].HydroFlag[icount]=st[i].Hydro[k];
						st[i].HydroAcc[j][icount]=getDistance(DISTSECOND,DISTTHIRD);

					}
				}
			}
		}
 	f[icount]=frame;
 	icount++;
 	frame++;
	}
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::distAndAngle(int DISTFIRST,int DISTSECOND,int DISTTHIRD,Data &dist,Data &angle)
{
	angle = getAngle(DISTFIRST,DISTSECOND,DISTTHIRD);
	DISTFIRST = DISTSECOND;
	DISTSECOND = DISTTHIRD;
	dist=getDistance(DISTFIRST,DISTSECOND);
}

void	Dcd::getNegPosNeg(IntVector &f,IntVector &pos,IntVector &neg,
					vector<ComplexSaltStruct> &poscsalt,
					vector<ComplexSaltStruct> &negcsalt,
					vector<NegPosNegTable> &PosTable)
{
	PosTable.resize(neg.size()*pos.size()*neg.size());
	f.resize(nframes());
	for(int i=0;i<PosTable.size();i++)
	{
		PosTable[i].resize(f.size());
		for(int j=0;j<PosTable[i].size();j++) PosTable[i][j].id1=PosTable[i][j].id2=PosTable[i][j].id3=0;
		for(int j=0;j<PosTable[i].size();j++) PosTable[i][j].dist1=PosTable[i][j].dist2=1e+100;
	}

	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;

	extern vector<atom> table;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
		int ipos = -1;
		for(int i=0;i<neg.size();i++)
		{
			int pos1=neg[i]-1;
			for(int j=pos1;table[j].res_id == table[pos1].res_id;j++)
			{
				 int ifound=0;
				 for(int i1=0;i1<negcsalt.size();i1++)
				 {
					if(negcsalt[i1].res_name==table[j].res_name
					&& negcsalt[i1].atom_name == table[j].atom_name) 
					{
						ifound=1;
						break;
					}
				 }
				 if(!ifound) continue;

				for(int k=0;k<pos.size();k++)
				{
					int pos2=pos[k]-1;
					if(table[pos1].res_id==table[pos2].res_id) continue;
					for(int l=pos2;table[l].res_id==table[pos2].res_id;l++)
					{
						if(table[l].atom_type != "NH3")
						{
				 			int ifound=0;
				 			for(int i1=0;i1<poscsalt.size();i1++)
				 			{
							if(poscsalt[i1].res_name==table[l].res_name
							&& poscsalt[i1].atom_name == table[l].atom_name) 
							{
								ifound=1;
								break;
							}
							}
				 		}
						if(!ifound) continue;

						for(int m=0;m<neg.size();m++)
						{
							int pos3=neg[m]-1;
							if(table[pos1].res_id==table[pos2].res_id) 
								continue;
							if(table[pos1].res_id==table[pos3].res_id) 
								continue;
							if(table[pos2].res_id==table[pos3].res_id) 
								continue;
							if(pos1==pos3) continue;
							for(int m1=pos3;table[m1].res_id == table[pos3].res_id;m1++)
							{
								ipos = (i+1)*(k+1)*(m+1)-1;
								ipos = i * (neg.size()*pos.size())+
								k * neg.size()+m;
								if(table[m1].chain_id == table[j].chain_id &&
								table[m1].res_id==table[j].res_id) continue;
				 				int ifound=0;
				 				for(int i1=0;i1<negcsalt.size();i1++)
				 				{
									if(negcsalt[i1].res_name==table[m1].res_name
									&& negcsalt[i1].atom_name == table[m1].atom_name) 
									{
										ifound=1;
										break;
									}
				 				}				
				 				if(!ifound) continue;


								int DISTFIRST,DISTSECOND,DISTTHIRD;
								DISTFIRST=j;
								DISTSECOND=l;
								DISTTHIRD=m1;
								if(j==m1) continue;
								Data d1=getDistance(DISTFIRST,DISTSECOND);
								Data d2=getDistance(DISTSECOND,DISTTHIRD);
								if(d1<PosTable[ipos][icount].dist1)
								{
									PosTable[ipos][icount].dist1=d1;
									PosTable[ipos][icount].id1=j+1;
									PosTable[ipos][icount].id2=l+1;
								}
								if(d2<PosTable[ipos][icount].dist2)
								{
									PosTable[ipos][icount].dist2=d2;
									PosTable[ipos][icount].id2=l+1;
									PosTable[ipos][icount].id3=m1+1;
								}
								PosTable[ipos][icount].neg_id1=i;
								PosTable[ipos][icount].pos_id=k;
								PosTable[ipos][icount].neg_id2=m;
							}
						}
					}
				}
			}
		}
		f[icount]=frame;
 		frame++;
 		icount++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getPosNegPos(IntVector &f,IntVector &pos,IntVector &neg,
					vector<ComplexSaltStruct> &poscsalt,
					vector<ComplexSaltStruct> &negcsalt,
					vector<PosNegPosTable> &PosTable)
{
	PosTable.resize(pos.size()*neg.size()*pos.size());
	f.resize(nframes());
	for(int i=0;i<PosTable.size();i++)
	{
		PosTable[i].resize(f.size());
		for(int j=0;j<PosTable[i].size();j++) PosTable[i][j].pos_id1=PosTable[i][j].id2=PosTable[i][j].id3=0;
		for(int j=0;j<PosTable[i].size();j++) PosTable[i][j].dist1=PosTable[i][j].dist2=1e+100;
	}

	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;

	extern vector<atom> table;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
		int ipos = -1;
		for(int i=0;i<pos.size();i++)
		{
			int pos1=pos[i]-1;
			for(int j=pos1;table[j].res_id == table[pos1].res_id;j++)
			{
				if(table[j].atom_type != "NH3")
				{
				 int ifound=0;
				 for(int i1=0;i1<poscsalt.size();i1++)
				 {
					if(poscsalt[i1].res_name==table[j].res_name
					&& poscsalt[i1].atom_name == table[j].atom_name) 
					{
						ifound=1;
						break;
					}
				 }
				 if(!ifound) continue;
				}
				for(int k=0;k<neg.size();k++)
				{
					int pos2=neg[k]-1;
					if(table[pos1].res_id == table[pos2].res_id) 
					{
						//ADDITION
						continue;
					}
					for(int l=pos2;table[l].res_id==table[pos2].res_id;l++)
					{
						int ifound=0;
						for(int i1=0;i1<negcsalt.size();i1++)
						{
							if(negcsalt[i1].res_name==table[l].res_name
							&& negcsalt[i1].atom_name == table[l].atom_name) 
							{
								ifound=1;
								break;
							}
						}
						if(!ifound) continue;

						for(int m=0;m<pos.size();m++)
						{
							int pos3=pos[m]-1;
							if(pos1==pos3) continue;
							if(table[pos1].res_id==table[pos2].res_id
							 || table[pos2].res_id==table[pos3].res_id
							 || table[pos3].res_id==table[pos1].res_id)
							continue;
							for(int m1=pos3;table[m1].res_id == table[pos3].res_id;m1++)
							{
								ipos = (i+1)*(k+1)*(m+1)-1;
								ipos = i * (neg.size()*pos.size())+
								k * pos.size()+m;
								if(table[m1].chain_id == table[j].chain_id &&
								table[m1].res_id==table[j].res_id) continue;
								if(table[m1].atom_type != "NH3")
								{
				 				int ifound=0;
				 				for(int i1=0;i1<poscsalt.size();i1++)
				 				{
								if(poscsalt[i1].res_name==table[m1].res_name
								&& poscsalt[i1].atom_name == table[m1].atom_name) 
								{
									ifound=1;
									break;
								}
				 				}
				 				if(!ifound) continue;
								}
								int DISTFIRST,DISTSECOND,DISTTHIRD;
								DISTFIRST=j;
								DISTSECOND=l;
								DISTTHIRD=m1;
								if(j==m1) continue;
								Data d1=getDistance(DISTFIRST,DISTSECOND);
								Data d2=getDistance(DISTSECOND,DISTTHIRD);
								if(d1<PosTable[ipos][icount].dist1)
								{
									PosTable[ipos][icount].dist1=d1;
									PosTable[ipos][icount].id1=j+1;
									PosTable[ipos][icount].id2=l+1;
								}
								if(d2<PosTable[ipos][icount].dist2)
								{
									PosTable[ipos][icount].dist2=d2;
									PosTable[ipos][icount].id2=l+1;
									PosTable[ipos][icount].id3=m1+1;
								}
								PosTable[ipos][icount].pos_id1=i;
								PosTable[ipos][icount].neg_id=k;
								PosTable[ipos][icount].pos_id2=m;
							}
						}
					}
				}
			}
		}
		f[icount]=frame;
 		frame++;
 		icount++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getPosNeg(IntVector &f,IntVector &pos,IntVector &neg,
	vector<ComplexSaltStruct> &poscsalt,
	vector<ComplexSaltStruct> &negcsalt,
	vector<PosTable2> &PosTable)
{
	PosTable.resize(pos.size()*neg.size());
	f.resize(nframes());
	for(int i=0;i<PosTable.size();i++)
	{
		PosTable[i].resize(f.size());
		for(int j=0;j<PosTable[i].size();j++) 
		{
			PosTable[i][j].dist=100.0;
			PosTable[i][j].id1=0;
			PosTable[i][j].id2=0;
		}
	}

	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;

	extern vector<atom> table;
	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<pos.size();i++)
		{
			int pos1=pos[i]-1;
			for(int j=pos1;table[j].res_id == table[pos1].res_id;j++)
			{
				if(table[j].atom_type != "NH3")
				{
				 int ifound=0;
				 for(int i1=0;i1<poscsalt.size();i1++)
				 {
					if(poscsalt[i1].res_name==table[j].res_name
					&& poscsalt[i1].atom_name == table[j].atom_name) 
					{
						ifound=1;
						break;
					}
				 }
				 if(!ifound) continue;
				}
				for(int k=0;k<neg.size();k++)
				{
					int pos2=neg[k]-1;
					int negfound=0;
					for(int l=pos2;table[l].res_id==table[pos2].res_id;l++)
					{
						int ifound=0;
						for(int i1=0;i1<negcsalt.size();i1++)
						{
							if(negcsalt[i1].res_name==table[l].res_name
							&& negcsalt[i1].atom_name == table[l].atom_name) 
							{
								ifound=1;
								break;
							}
						}
						if(!ifound) continue;
						negfound=1;
						int DISTFIRST,DISTSECOND;
						DISTFIRST=j;
						DISTSECOND=l;
						Data d=getDistance(DISTFIRST,DISTSECOND);
						if(PosTable[i*neg.size()+k][icount].dist<0 || d<PosTable[i*neg.size()+k][icount].dist)
						{
							PosTable[i*neg.size()+k][icount].dist=d;
							PosTable[i*neg.size()+k][icount].id1=j+1;
							PosTable[i*neg.size()+k][icount].id2=l+1;
						}
					}
				}
			}
		}
		f[icount]=frame;
 		frame++;
 		icount++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

/*	Return the minimum distance between the atoms belonging in the first group (group1)
 *	and the atoms belonging in the second group (group2)
 * */
void	Dcd::getGroupDistance(IntVector2 &group1,IntVector2 &group2,
				IntVector &f,DoubleVector2 &d,
				IntVector2 &atom1,IntVector2 &atom2)
{
	f.resize(nframes());
	d.resize(group1.size());
	atom1.resize(group1.size());
	atom2.resize(group1.size());
	for(int i=0;i<group1.size();i++)
	{
		d[i].resize(f.size());
		atom1[i].resize(f.size());
		atom2[i].resize(f.size());
	}

	int icount =0;
	int FIRST = frameStart;
	int LAST  = frameEnd;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int STEP=frameStep;
	int itotalframes=0;

	/**/

	while(frame<=LAST)
	{
	 if((frame-firstframe)%STEP==0)
	 {
		if(frame%1000==0) printf("frame = %d\n",frame);
		if(!readCAs()) break;
		for(int k=0;k<group1.size();k++)
		{
	 		int DISTFIRST,DISTSECOND;
			Data ddmin=1e+100;
			int    minatom1=0,minatom2=0;
			for(int i=0;i<group1[k].size();i++)
			{
				for(int j=0;j<group2[k].size();j++)
				{
					DISTFIRST=group1[k][i]-1;
					DISTSECOND=group2[k][j]-1;
					Data value=getDistance(DISTFIRST,DISTSECOND);
					if(value <= ddmin) 
					{
						ddmin = value;
						minatom1 = group1[k][i];
						minatom2 = group2[k][j];
					}
				}
			}
	 		d[k][icount]=ddmin;
			atom1[k][icount] = minatom1;
			atom2[k][icount] = minatom2;
		}
	 	f[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else
	 if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}


Data	Dcd::countAtoms(int start,Data distance,vector<int> &list,vector<Data> &vframe)
{
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	Data total_icount=0;
	/**/

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
	 	Data d,d1;
	 	int DISTFIRST,DISTSECOND,DISTTHIRD,DISTFOURTH;
	 	DISTFIRST=start-1;
		int dicount=0;
		for(int i=0;i<list.size();i++)
		{
			DISTSECOND=list[i]-1;	
	         	d=getDistance(DISTFIRST,DISTSECOND);
			if(d<distance) dicount++;
		}
		Data rho = list.size()*1.0/vframe[frame-1];
		total_icount+=dicount/rho;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
	return total_icount/icount;	
}

Dcd::~Dcd()
{
	if(dcd!=-1) 
	{
		close(dcd);
		free(dcd_frame);
		free(wdcd_frame);
	}
}

void	Dcd::getAronh(IntVector &f,vector<DoubleVector> &Dist,vector<DoubleVector> &Dist2,
			vector<DoubleVector> &Angle,
			vector<stack_table> &stable,vector<NHStruct> &nhpos,
			Data stack_distance,Data stack_angle)
{
	f.resize(nframes());
	Dist.resize(nhpos.size() * stable.size());
	Dist2.resize(nhpos.size() * stable.size());
	Angle.resize(nhpos.size() * stable.size());
	for(int i=0;i<Dist.size();i++)
	{
		Dist[i].resize(f.size());
		Dist2[i].resize(f.size());
		Angle[i].resize(f.size());
	}
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	Data total_icount=0;
	Point avg;
	Point hpoint;
	Point npoint;
	Point pt1;
	Point pt2;
	Point pt3;
	Point inter;
	Point normal;
	/**/

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int j=0;j<nhpos.size();j++)
		{
			if(nhpos[j].hpos==-1 || nhpos[j].npos==-1) continue;
			hpoint=makePoint(nhpos[j].hpos);
			npoint=makePoint(nhpos[j].npos);	
			for(int i=0;i<stable.size();i++)
			{
			Data A,B,C,D;
			pt1=makePoint(stable[i].id1-1);
			pt2=makePoint(stable[i].id2-1);
			pt3=makePoint(stable[i].id3-1);

			getPlaneElements(pt1,pt2,pt3,A,B,C,D);
			getPlaneNormal(pt1,pt2,pt3,normal);
			getInterSectionPoint(npoint,hpoint,A,B,C,D,inter);
			IntVector list;
			list.resize(3);list[0]=stable[i].id1-1;list[1]=stable[i].id2-1;list[2]=stable[i].id3-1;	
			avg=makeCenter(list,0);
			Data d=getPointDistance(npoint,avg);
			Data d1=getPointDistance(hpoint,avg);
			Data a=getAngle(npoint,hpoint,avg);		
			//Data a=fabs(getAngleBetweenLines(npoint,hpoint,inter,normal));
			
			Dist[j*stable.size()+i][icount]=d;
			Dist2[j*stable.size()+i][icount]=d1;
			Angle[j*stable.size()+i][icount]=a;
			}
		}
	
		f[icount]=frame;	
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getAropos(IntVector &f,vector<DoubleVector> &Dist,vector<DoubleVector> &Angle,
			vector<stack_table> &stable,vector<AroposStruct> &nhpos,
			Data stack_distance,Data stack_angle)
{
	f.resize(nframes());
	Dist.resize(nhpos.size() * stable.size());
	Angle.resize(nhpos.size() * stable.size());
	for(int i=0;i<Dist.size();i++)
	{
		Dist[i].resize(f.size());
		Angle[i].resize(f.size());
	}
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	Data total_icount=0;
	Point avg;
	Point hpoint;
	Point npoint;
	Point pt1;
	Point pt2;
	Point pt3;
	Point inter;
	Point normal;
	/**/

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int j=0;j<nhpos.size();j++)
		{
			npoint = makePoint(nhpos[j].npos);
			for(int i=0;i<stable.size();i++)
			{
			IntVector list;
			list.resize(3);list[0]=stable[i].id1-1;list[1]=stable[i].id2-1;list[2]=stable[i].id3-1;	
			avg=makeCenter(list,0);
			Data d=getPointDistance(npoint,avg);
			Dist[j*stable.size()+i][icount]=d;

			Data A,B,C,D;
			pt1=makePoint(stable[i].id1-1);
			pt2=makePoint(stable[i].id2-1);
			pt3=makePoint(stable[i].id3-1);

			getPlaneElements(pt1,pt2,pt3,A,B,C,D);
			getPlaneNormal(pt1,pt2,pt3,normal);

			Angle[j*stable.size()+i][icount] = 1000;
			for(int k=0;k<nhpos[j].hpos.size();k++)
			{
			hpoint = makePoint(nhpos[j].hpos[k]);
			Data a=getAngle(npoint,hpoint,avg);		
			/*
			getInterSectionPoint(npoint,hpoint,A,B,C,D,inter);
			Data a=fabs(getAngleBetweenLines(npoint,hpoint,inter,normal));
			*/
			if(a<Angle[j*stable.size()+i][icount]) Angle[j*stable.size()+i][icount]=a;
			}
			}
		}
	
		f[icount]=frame;	
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNoe(IntVector &f,vector<NoeStruct> &noe)
{
	f.resize(nframes());
	for(int i=0;i<noe.size();i++)
	{
		noe[i].D.resize(f.size());
	}
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<noe.size();i++)
		{
			int DISTFIRST=noe[i].atom1-1;
			int DISTSECOND=noe[i].atom2-1;
			Data d=getDistance(DISTFIRST,DISTSECOND);
			noe[i].D[icount]=d;
		}
		f[icount]=frame;	
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

extern Data rmsd(Data *v1, Data *v2, int N, Data *mtx,Data *v3);
typedef struct
{
  float m[4][4];
} MATRIX;

void	Dcd::getRmsd(IntVector &AtomPos,PointVector &pdbpos,IntVector &f,DoubleVector &d)
{
	f.resize(nframes());
	printf ("F=%d\n", f.size() ); //here
	d.resize(f.size());
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	int N=AtomPos.size();
	Data *v1=new Data[N * 3];
	Data *v2=new Data[N * 3];
	Data *v3=new Data[N * 3];
	Data centerB[3]={0,0,0};
	
	for(int i=0;i<N;i++)
	{
		v2[i*3+0]=pdbpos[AtomPos[i]].x;
		v2[i*3+1]=pdbpos[AtomPos[i]].y;
		v2[i*3+2]=pdbpos[AtomPos[i]].z;
		centerB[0]+=v2[i*3+0];
		centerB[1]+=v2[i*3+1];
		centerB[2]+=v2[i*3+2];
	}

	centerB[0]/=N;centerB[1]/=N;centerB[2]/=N;
	MATRIX mtx;

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		Data distance=0;
		Data centerA[3]={0,0,0};
		for(int i=0;i<N;i++)
		{
			v1[i*3+0]=CAs[AtomPos[i]][0];
			v1[i*3+1]=CAs[AtomPos[i]][1];
			v1[i*3+2]=CAs[AtomPos[i]][2];
			centerA[0]+=v1[i*3+0];
			centerA[1]+=v1[i*3+1];
			centerA[2]+=v1[i*3+2];
		}
		centerA[0]/=N;centerA[1]/=N;centerA[2]/=N;
			
		for(int i=0;i<N;i++)
		{
			v1[i*3+0]=v1[i*3+0]+(centerB[0]-centerA[0]);
			v1[i*3+1]=v1[i*3+1]+(centerB[1]-centerA[1]);
			v1[i*3+2]=v1[i*3+2]+(centerB[2]-centerA[2]);
		}
		
		double old_distance=1e+10;
		double diff=0.0;
		for(int iters=1;iters<=10;iters++)
		{
		distance=rmsd(v1, v2, N, NULL,v3);
		diff=fabs(distance-old_distance);
		if(diff<0.1) break;
		old_distance=distance;
		}
			
		f[icount]=frame;	
		d[icount]=distance;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
	delete[] v1;
	delete[] v2;
	delete[] v3;
}

void	Dcd::getRmsd(IntVector &AtomPos,IntVector &f,DoubleVector2 &d)
{
	f.resize(nframes());
	d.resize(f.size());
	for(int i=0;i<f.size();i++) d[i].resize(f.size());

	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	int N=AtomPos.size();

	Data *v1=new Data[N * 3];

	Data **v2=new Data*[f.size()];
	Data	*v3=new Data[N*3];
	for(int i=0;i<f.size();i++)
		v2[i]=new Data[N * 3];
	Data **centerB=new Data*[f.size()];
	for(int i=0;i<f.size();i++)
		centerB[i]=new Data[3];
	
	
	int iicount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
	
		centerB[iicount][0]=0;centerB[iicount][1]=0;centerB[iicount][2]=0;
		for(int i=0;i<N;i++)
		{
			v2[iicount][i*3+0]=CAs[AtomPos[i]][0];
			v2[iicount][i*3+1]=CAs[AtomPos[i]][1];
			v2[iicount][i*3+2]=CAs[AtomPos[i]][2];
			centerB[iicount][0]+=v2[iicount][i*3+0];
			centerB[iicount][1]+=v2[iicount][i*3+1];
			centerB[iicount][2]+=v2[iicount][i*3+2];
		}
		centerB[iicount][0]/=N;centerB[iicount][1]/=N;centerB[iicount][2]/=N;
	 }
	 else if(!skipFrame(frame)) break;
	 iicount++;
	 frame++;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);

	MATRIX mtx;

	frame= frameStart;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		Data distance=0;
		
		for(int ik=0;ik<icount;ik++)
		{
			Data centerA[3]={0,0,0};
			for(int i=0;i<N;i++)
			{
			v1[i*3+0]=CAs[AtomPos[i]][0];
			v1[i*3+1]=CAs[AtomPos[i]][1];
			v1[i*3+2]=CAs[AtomPos[i]][2];
			centerA[0]+=v1[i*3+0];
			centerA[1]+=v1[i*3+1];
			centerA[2]+=v1[i*3+2];
			}
			centerA[0]/=N;centerA[1]/=N;centerA[2]/=N;
			
			for(int i=0;i<N;i++)
			{
			v1[i*3+0]=v1[i*3+0]+(centerB[ik][0]-centerA[0]);
			v1[i*3+1]=v1[i*3+1]+(centerB[ik][1]-centerA[1]);
			v1[i*3+2]=v1[i*3+2]+(centerB[ik][2]-centerA[2]);
			}
		
			double old_distance=1e+10;
			double diff=0.0;
			for(int iters=1;iters<=10;iters++)
			{
			distance=rmsd(v1, v2[ik], N, NULL,v3);
			diff=fabs(distance-old_distance);
			if(diff<0.1) break;
			old_distance=distance;
			}
			
			f[icount]=frame;	
			d[icount][ik]=distance;
		}
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
	delete[] v1;
	for(int i=0;i<f.size();i++)
	{
		delete[] v2[i];
		delete[] centerB[i];
	}
	delete[] v2;
	delete[] v3;
	delete[] centerB;
}
void	Dcd::getNoeAngle(IntVector &f,vector<NoeAngleStruct> &noe)
{
	f.resize(nframes());
	for(int i=0;i<noe.size();i++)
	{
		noe[i].D.resize(f.size());
	}
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<noe.size();i++)
		{
			int DISTFIRST=noe[i].atom1-1;
			int DISTSECOND=noe[i].atom2-1;
			int DISTTHIRD=noe[i].atom3-1;
			Data d=getAngle(DISTFIRST,DISTSECOND,DISTTHIRD);
			noe[i].D[icount]=d;
		}
		f[icount]=frame;	
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNoeDihedral(IntVector &f,vector<NoeDihedralStruct> &noe)
{
	f.resize(nframes());
	for(int i=0;i<noe.size();i++)
	{
		noe[i].D.resize(f.size());
	}
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<noe.size();i++)
		{
			int DISTFIRST=noe[i].atom1-1;
			int DISTSECOND=noe[i].atom2-1;
			int DISTTHIRD=noe[i].atom3-1;
			int DISTFOURTH=noe[i].atom4-1;
			Data d=dihedral(CAs[DISTFIRST],CAs[DISTSECOND],CAs[DISTTHIRD],CAs[DISTFOURTH]);
			noe[i].D[icount]=d;
		}
		f[icount]=frame;	
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNoeCenter(IntVector &atom1,IntVector &atom2,IntVector &F,DoubleVector &D)
{
	F.resize(nframes());
	D.resize(F.size());
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	extern int weight_flag;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		Data distance;
		Point center1,center2;
		center1=makeCenter(atom1,weight_flag);
		center2=makeCenter(atom2,weight_flag);
		distance=getPointDistance(center1,center2);	
		F[icount]=frame;	
		D[icount]=distance;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNoeCenterAngle(IntVector &atom1,IntVector &atom2,IntVector &atom3,
		IntVector &F,DoubleVector &D)
{
	F.resize(nframes());
	D.resize(F.size());
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;

	extern int weight_flag;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		Data distance;
		Point center1,center2,center3;
		center1=makeCenter(atom1,0);
		center2=makeCenter(atom2,0);
		center3=makeCenter(atom3,0);
		distance=getAngle(center1,center2,center3);	
		F[icount]=frame;	
		D[icount]=distance;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNoeCenterDihedral(IntVector &atom1,IntVector &atom2,IntVector &atom3,
		IntVector &atom4,IntVector &F,DoubleVector &D)
{
	F.resize(nframes());
	D.resize(F.size());
	int icount =0;
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	extern int weight_flag;

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		Data distance;
		Point center1,center2,center3,center4;
		center1=makeCenter(atom1,0);
		center2=makeCenter(atom2,0);
		center3=makeCenter(atom3,0);
		center4=makeCenter(atom4,0);
		distance=dihedral(center1,center2,center3,center4);	
		F[icount]=frame;	
		D[icount]=distance;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::fetchAtoms(int frame,vector<DoubleVector> &D)
{
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	lseek(dcd,(off_t)(bytes_per_set*(frame-1)),SEEK_CUR);
	readCAs();
	D.resize(natoms);
	
	for(int i=0;i<natoms;i++)
	{
		D[i].resize(3);
		D[i][0]=CAs[i][0];
		D[i][1]=CAs[i][1];
		D[i][2]=CAs[i][2];
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getPdo(IntVector &F,vector<PdoStruct> &Pdo)
{
	F.resize(nframes());
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	for(int i=0;i<Pdo.size();i++)
	{
		Pdo[i].D1.resize(F.size());
		Pdo[i].D2.resize(F.size());
		Pdo[i].A1.resize(F.size());
		Pdo[i].A2.resize(F.size());
	}
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<Pdo.size();i++)
		{
			Pdo[i].D1[icount]=getDistance(Pdo[i].dist_atom1-1,Pdo[i].dist_atom2-1);
			Pdo[i].D2[icount]=getDistance(Pdo[i].dist_atom3-1,Pdo[i].dist_atom4-1);
			Pdo[i].A1[icount]=dihedral(Pdo[i].dihe_atom1-1,Pdo[i].dihe_atom2-1,
				Pdo[i].dihe_atom3-1,Pdo[i].dihe_atom4-1);
			Pdo[i].A2[icount]=dihedral(Pdo[i].dihe_atom5-1,Pdo[i].dihe_atom6-1,
				Pdo[i].dihe_atom7-1,Pdo[i].dihe_atom8-1);
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getHydroSolve(IntVector &waterlist,IntVector &F,vector<HydroSolveStruct> &st,
		Data distance,Data angle)
{
	F.resize(nframes());
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<st.size();i++)
		{
			int DISTFIRST=st[i].atom1-1;
			int DISTTHIRD=st[i].atom2-1;
			double d=getDistance(DISTFIRST,DISTTHIRD);
			if(d>=2.5 * distance) continue;
			for(int j=0;j<waterlist.size();j++)
			{
				int DISTSECOND=waterlist[j]-1;
				double d1=getDistance(DISTFIRST,DISTSECOND);
				double d2=getDistance(DISTSECOND,DISTTHIRD);
				double a=getAngle(DISTFIRST,DISTSECOND,DISTTHIRD);
				if(d1<=distance && d2<=distance && fabs(a)>=angle)
				{
					int k=st[i].Frame.size();
					st[i].Frame.resize(k+1);
					st[i].Water.resize(k+1);
					st[i].Distance1.resize(k+1);
					st[i].Distance2.resize(k+1);
					st[i].Frame[k]=icount;
					st[i].Water[k]=waterlist[j];
					st[i].Distance1[k]=d1;
					st[i].Distance2[k]=d2;
				}
			}
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}
void	Dcd::getCenterNear(IntVector &F,IntVector &atom1,IntVector &atom2,
				IntVector &Index,DoubleVector &distance)
{
	F.resize(nframes());
	distance.resize(F.size());
	Index.resize(F.size());
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	Point center;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		center=makeCenter(atom1,0);
		Point pt;
		distance[icount]=1e+100;
		for(int i=0;i<atom2.size();i++)
		{
			pt=makePoint(atom2[i]-1);
			double d=getPointDistance(center,pt);
			if(d<distance[icount])
			{
				distance[icount]=d;
				Index[icount]=atom2[i];
			}
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getIwcn(double d,IntVector &F,vector<IwcnStruct> &iwcn)
{
	int n=nframes();
	for(int i=0;i<iwcn.size();i++)
	{
		iwcn[i].count.resize(n);
		iwcn[i].atomList.resize(n);
		for(int j=0;j<n;j++) 
		{
			iwcn[i].count[j]=0;
			iwcn[i].atomList[j].resize(0);
		}
	}
	F.resize(n);
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<iwcn.size();i++)
		{
			int DISTFIRST=iwcn[i].water-1;
			for(int j=0;j<iwcn[i].list.size();j++)
			{
				int DISTSECOND=iwcn[i].list[j]-1;
				if(DISTFIRST == DISTSECOND) continue;
				double dist=getDistance(DISTFIRST,DISTSECOND);
				if(dist<=d) 
				{
					iwcn[i].count[icount]++;
					iwcn[i].atomList[icount].push_back(DISTSECOND+1);
				}
			}
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::countFrames(double d,IntVector &F,vector<CountFramesStruct> &count)
{
	for(int i=0;i<count.size();i++) count[i].count=0;
	F.resize(nframes());
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<count.size();i++)
		{
			int DISTFIRST=count[i].atom1-1;
			int DISTSECOND=count[i].atom2-1;
			double dist=getDistance(DISTFIRST,DISTSECOND);
			if(dist<=d) count[i].count++;
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNear(IntVector &F,IntVector2 &List,DoubleVector2 &Dist)
{
	F.resize(nframes());
	Dist.resize(List.size()  * List.size());
	for(int i=0;i<Dist.size();i++) Dist[i].resize(F.size());
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		int dcount=0;
		for(int i=0;i<List.size();i++)
		{
			for(int j=0;j<i;j++)
			{
				dcount=i*List.size()+j;
				Dist[dcount][icount]=1e+100;
				for(int k1=0;k1<List[i].size();k1++)
				{
					for(int k2=0;k2<List[j].size();k2++)
					{
						int DISTFIRST=List[i][k1]-1;
						int DISTSECOND=List[j][k2]-1;
						double d=getDistance(DISTFIRST,DISTSECOND);
						if(d<Dist[dcount][icount]) 
							Dist[dcount][icount]=d;
					}
				}
			}
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getNear(IntVector &F,IntVector2 &List1,IntVector &List2,DoubleVector2 &Dist,IntVector &Vote)
{
	F.resize(nframes());
	Dist.resize(List1.size());
	for(int i=0;i<Dist.size();i++) Dist[i].resize(F.size());

	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<List1.size();i++)
		{
			Dist[i][icount]=1e+100;
			int imax=0;
			for(int k1=0;k1<List1.size();k1++)
			{
				for(int j=0;j<List2.size();j++)
				{
					int DISTFIRST=List1[i][k1]-1;
					int DISTSECOND=List2[j]-1;
					double d=getDistance(DISTFIRST,DISTSECOND);
					if(d<Dist[i][icount]) 
					{
						imax=j;
						Dist[i][icount]=d;
					}
				}
			}
			Vote[imax]++;
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getAvgPos(IntVector &F,IntVector &posAtom,DoubleVector2 &D,int N,
		int fitflag,PointVector &pdbpos)
{
	F.resize(nframes());
	int count;
	if(F.size() % N !=0) 
	count=1+F.size() / N;
	else	count=F.size() /N;
	D.resize(count);
	for(int k=0;k<D.size();k++)
	{
		D[k].resize(3*posAtom.size());
		for(int i=0;i<posAtom.size();i++)	
		{
			D[k][i*3+0]=0.0;
			D[k][i*3+1]=0.0;
			D[k][i*3+2]=0.0;
		}
	}
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	int ik=0;
	int slice=0;

	DoubleVector2 Dest;
	Dest.resize(CAs.size());
	for(int i=0;i<CAs.size();i++) Dest[i].resize(3);

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		if(pdbpos.size()==0 &&fitflag)
		{
			pdbpos.resize(CAs.size());
			for(int i=0;i<CAs.size();i++) 
			{
				pdbpos[i].x=CAs[i][0];
				pdbpos[i].y=CAs[i][1];
				pdbpos[i].z=CAs[i][2];
			}
		}
		/**FITTED FRAME**/
		if(fitflag) getFitted(fitflag,posAtom,pdbpos,CAs,Dest);
		for(int k=count-1;k>=slice;k--)
		for(int i=0;i<posAtom.size();i++)	
		{
		D[k][i*3+0]+=fitflag?Dest[posAtom[i]][0]:CAs[posAtom[i]][0];
		D[k][i*3+1]+=fitflag?Dest[posAtom[i]][1]:CAs[posAtom[i]][1];
		D[k][i*3+2]+=fitflag?Dest[posAtom[i]][2]:CAs[posAtom[i]][2];
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
		ik++;
		if(ik==N) {slice++;ik=0;}
	 }
	 else if(!skipFrame(frame)) break;
	}
	for(int k=0;k<count;k++)
	for(int i=0;i<posAtom.size();i++)	
	{
		D[k][i*3+0]/=((k+1)*N);
		D[k][i*3+1]/=((k+1)*N);
		D[k][i*3+2]/=((k+1)*N);
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getEntropyData(IntVector &F,DoubleVector3 &D,
			    	IntVector &posAtom,int N,
				int fitflag,PointVector &pdbpos)
{
	F.resize(nframes());
	int count;
	if(F.size() % N !=0) 
	count=1+F.size() / N;
	else	count=F.size() /N;
	D.resize(count);
	DoubleVector2 avg;
	avg.resize(count);
	getAvgPos(F,posAtom,avg,N,fitflag,pdbpos);
	for(int i=0;i<count;i++)
	{
		D[i].resize(3*posAtom.size());
		for(int j=0;j<D[i].size();j++)
		{
			D[i][j].resize(3*posAtom.size());
			for(int k=0;k<D[i][j].size();k++) 
			{
				D[i][j][k]=0.0;
			}
		}
	}

	int NN=posAtom.size();

	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	DoubleVector CAxyz;
	CAxyz.resize(3 * posAtom.size());
	int ik=0;
	int slice=0;
	DoubleVector2 Dest;
	Dest.resize(CAs.size());
	for(int i=0;i<CAs.size();i++) Dest[i].resize(3);
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		if(pdbpos.size()==0 && fitflag)
		{
			pdbpos.resize(CAs.size());
			for(int i=0;i<CAs.size();i++) 
			{
				pdbpos[i].x=CAs[i][0];
				pdbpos[i].y=CAs[i][1];
				pdbpos[i].z=CAs[i][2];
			}
		}
		/**FITTED FRAME**/
		if(fitflag) getFitted(fitflag,posAtom,pdbpos,CAs,Dest);
		for(int i=0;i<posAtom.size();i++)
		{
		CAxyz[i*3+0]=fitflag?Dest[posAtom[i]][0]:CAs[posAtom[i]][0];
		CAxyz[i*3+1]=fitflag?Dest[posAtom[i]][1]:CAs[posAtom[i]][1];
		CAxyz[i*3+2]=fitflag?Dest[posAtom[i]][2]:CAs[posAtom[i]][2];
		}

		for(int k=count-1;k>=slice;k--)
		for(int i=0;i<3*posAtom.size();i++)
		{
			for(int j=i;j<3*posAtom.size();j++)
			{
				float v=CAxyz[i]*CAxyz[j]+avg[k][i]*avg[k][j]
						-CAxyz[i]*avg[k][j]
						-avg[k][i]*CAxyz[j];
				D[k][i][j]+=v;	
			}
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
		ik++;
		if(ik==N) {slice++;ik=0;}
	 }
	 else if(!skipFrame(frame)) break;
	}

	extern vector<atom> table;
	Data dmin=1e+10;
	Data dmax=-1e+10;
	for(int k=0;k<D.size();k++)
	{
		for(int i=0;i<D[k].size();i++)
		{
			for(int j=i;j<D[k][i].size();j++)
			{
				Data wi=table[posAtom[i/3]].dnum2;
				Data wj=table[posAtom[j/3]].dnum2;
				D[k][i][j]*=sqrt(wi)*sqrt(wj);
				D[k][i][j]/=((k+1)*N);
				D[k][j][i]=D[k][i][j];
				if(k==count-1 && D[k][i][j]<dmin) dmin=D[k][i][j];
				if(k==count-1 && D[k][i][j]>dmax) dmax=D[k][i][j];
			}
		}
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::writeDcd(string fname,IntVector &posAtom,int fit)
{
	int wdcd;
	wdcd=open(fname.c_str(),O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP);
	if(wdcd == -1) 
	{
		psf_error("Can not open DCD file for write. Abort.");
	}
	unsigned int wauint,w_bytes_per_set;
	int whave_cell=0;
	int wcell_offset=0;
	char wdcd_head[92];
	char wdcd_head2[16];
	float *w_dcd_frame;
	char wtitle[81];

	for(int i=0;i<92;i++) wdcd_head[i]=dcd_head[i];
	for(int i=0;i<16;i++) wdcd_head2[i]=dcd_head2[i];
	if ( (unsigned int)(*((unsigned int *)(&wdcd_head[48]))) == 1 )
       {
               whave_cell=1;
               wcell_offset=14;
       }

	unsigned int ff=nframes();
	
	 *((unsigned int *)(&wdcd_head[8])) = (unsigned int)(ff);
	//write(wdcd,wdcd_head,(size_t)92);
	wauint=auint;
	//write(wdcd,&wauint,(size_t)4);
	//write(wdcd,&wauint,(size_t)4);
	strcpy(wtitle,title);
       	for(int i=0;i<wauint;i++)
       	{
	///	write(wdcd,&wtitle,(size_t)80);
       	}
	
	unsigned int pp=posAtom.size();
	 *((unsigned int *)(&wdcd_head2[8])) = (unsigned int)(pp);
	//write(wdcd,wdcd_head2,(size_t)16);
       if(whave_cell)
       {
	       w_bytes_per_set = (unsigned int)(3 * ( 4 * posAtom.size() + 8 ) + 56);
       }
       else
       {
	        w_bytes_per_set = (unsigned int)(3 * ( 4 * posAtom.size() + 8 ));
       }
       w_dcd_frame = make_vector( 1, w_bytes_per_set / 4 );
	int afterheader = lseek( dcd, (off_t)(0), SEEK_CUR );

	char *outdcdhead=(char *)malloc((size_t)afterheader);
	lseek( dcd, (off_t)(0), SEEK_SET );
	read( dcd, outdcdhead, (size_t)(afterheader) );
	 *((unsigned int *)(&outdcdhead[8])) = (unsigned int)(ff);
	 *((unsigned int *)(&outdcdhead[afterheader-8])) = (unsigned int)(posAtom.size());

        write( wdcd, outdcdhead, (size_t)(afterheader) );


	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;

	int N=posAtom.size();

	/**DEBUG**/
	extern PointVector pdbpos;
	DoubleVector2 Dest;
	if(fit)
	{
		Dest.resize(CAs.size());
		for(int i=0;i<CAs.size();i++) Dest[i].resize(3);
	}
	/**DEBUG**/


	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		int i,atomno;

		if(fit)
		if(pdbpos.size()==0)
		{
			pdbpos.resize(CAs.size());
			for(int i=0;i<CAs.size();i++) 
			{
				pdbpos[i].x=CAs[i][0];
				pdbpos[i].y=CAs[i][1];
				pdbpos[i].z=CAs[i][2];
			}
		}
		if(fit) getFitted(fit,posAtom,pdbpos,CAs,Dest);

		if ( whave_cell )
                {
                        *((double *)(&w_dcd_frame[2] +  0)) = cell_a;
                        *((double *)(&w_dcd_frame[2] +  2)) = cell_gamma;
                        *((double *)(&w_dcd_frame[2] +  4)) = cell_b;
                        *((double *)(&w_dcd_frame[2] +  6)) = cell_beta;
                        *((double *)(&w_dcd_frame[2] +  8)) = cell_alpha;
                        *((double *)(&w_dcd_frame[2] + 10)) = cell_c;
                }


		for ( i=0 ; i < posAtom.size() ; i++ )
		{
	      		const int CELL_OFFSET = wcell_offset;
			int atomno=i;
			int offset=posAtom.size()+2;
	      		w_dcd_frame[ CELL_OFFSET + 2+atomno]=fit?Dest[posAtom[i]][0]:CAs[posAtom[i]][0];
	      		w_dcd_frame[ CELL_OFFSET + 2+atomno+offset]=fit?Dest[posAtom[i]][1]:CAs[posAtom[i]][1];
	      		w_dcd_frame[ CELL_OFFSET + 2+atomno+2*offset]=fit?Dest[posAtom[i]][2]:CAs[posAtom[i]][2];
		}
 		int bytes=write( wdcd, (void *)(&w_dcd_frame[1]), (size_t)(w_bytes_per_set) );
		++frame;
		icount++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
	free(w_dcd_frame);
	close(wdcd);
}

void	Dcd::getEnergy(IntVector &F,DoubleVector2 &E,IntVector &atomPos,vector<string> &list)
{
	extern string parfile;
	ParFile pf(parfile);
	F.resize(nframes());
	E.resize(list.size());
	for(int i=0;i<E.size();i++) E[i].resize(F.size());

	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	extern vector<atom> table;
	IntVector	bondFlag;
	bondFlag.resize(bondstruct.size());
	for(int j=0;j<bondstruct.size();j++)
	{
		int id1=bondstruct[j].first-1;
		int id2=bondstruct[j].second-1;
		int ifound1=0,ifound2=0;
		for(int k=0;k<atomPos.size();k++)
		{
			if(atomPos[k]==id1) ifound1=1;
			if(atomPos[k]==id2) ifound2=1;
			if(ifound1 && ifound2) break;
		}
		if(!ifound1 || !ifound2) bondFlag[j]=0;else bondFlag[j]=1;
	}
	IntVector	angleFlag;
	angleFlag.resize(anglestruct.size());
	for(int j=0;j<anglestruct.size();j++)
	{
		int id1=anglestruct[j].first-1;
		int id2=anglestruct[j].second-1;
		int id3=anglestruct[j].third-1;
		int ifound1=0,ifound2=0,ifound3=0;
		for(int k=0;k<atomPos.size();k++)
		{
			if(atomPos[k]==id1) ifound1=1;
			if(atomPos[k]==id2) ifound2=1;
			if(atomPos[k]==id3) ifound3=1;
			if(ifound1 && ifound2 && ifound3) break;
		}
		if(!ifound1 || !ifound2 ||!ifound3) angleFlag[j]=0;else angleFlag[j]=1;
	}

	IntVector	diheFlag;
	diheFlag.resize(dihedralstruct.size());
	for(int j=0;j<dihedralstruct.size();j++)
	{
		int id1=dihedralstruct[j].first-1;
		int id2=dihedralstruct[j].second-1;
		int id3=dihedralstruct[j].third-1;
		int id4=dihedralstruct[j].fourth-1;
		int ifound1=0,ifound2=0,ifound3=0,ifound4=0;
		for(int k=0;k<atomPos.size();k++)
		{
			if(atomPos[k]==id1) ifound1=1;
			if(atomPos[k]==id2) ifound2=1;
			if(atomPos[k]==id3) ifound3=1;
			if(atomPos[k]==id4) ifound4=1;
			if(ifound1 && ifound2 && ifound3 && ifound4) break;
		}
		if(!ifound1 || !ifound2 ||!ifound3 ||!ifound4) diheFlag[j]=0;else diheFlag[j]=1;
	}

	IntVector	impFlag;
	impFlag.resize(improperstruct.size());
	for(int j=0;j<improperstruct.size();j++)
	{
		int id1=improperstruct[j].first-1;
		int id2=improperstruct[j].second-1;
		int id3=improperstruct[j].third-1;
		int id4=improperstruct[j].fourth-1;
		int ifound1=0,ifound2=0,ifound3=0,ifound4=0;
		for(int k=0;k<atomPos.size();k++)
		{
			if(atomPos[k]==id1) ifound1=1;
			if(atomPos[k]==id2) ifound2=1;
			if(atomPos[k]==id3) ifound3=1;
			if(atomPos[k]==id4) ifound4=1;
			if(ifound1 && ifound2 && ifound3 && ifound4) break;
		}
		if(!ifound1 || !ifound2 ||!ifound3 ||!ifound4) impFlag[j]=0;else impFlag[j]=1;
	}
	IntVector2 vdwlist;
	int iiflag=0;
	for(int i=0;i<list.size();i++) if(list[i]=="vdw" || list[i]=="elect") {iiflag=1;break;}
	if(iiflag)
	for(int i=0;i<atomPos.size();i++)
	{
		for(int j=i+1;j<atomPos.size();j++)
		{
			if(i==j) continue;
			int id1=table[atomPos[i]].atom_id;
			int id2=table[atomPos[j]].atom_id;
			int iflag=0;
			for(int k=0;k<bondstruct.size();k++)
			{
				int i1=0;int i2=0;
				if(bondstruct[k].first==id1 || bondstruct[k].second==id1) i1=1;
				if(bondstruct[k].first==id2 || bondstruct[k].second==id2) i2=1;
				if(i1 && i2) {iflag=1;break;}
			}
			if(iflag) continue;
			for(int k=0;k<anglestruct.size();k++)
			{
				int i1=0,i2=0;
				if(anglestruct[k].first==id1 || anglestruct[k].second==id1 || 
				   anglestruct[k].third==id1) i1=1;
				if(anglestruct[k].first==id2 || anglestruct[k].second==id2 || 
				   anglestruct[k].third==id2) i2=1;
				if(i1 && i2) {iflag=1;break;}
			}
			if(iflag) continue;
			/*
			for(int k=0;k<dihedralstruct.size();k++)
			{
				int i1=0,i2=0;
				if(dihedralstruct[k].first==id1 || dihedralstruct[k].second==id1 || 
				   dihedralstruct[k].third==id1 || dihedralstruct[k].fourth==id1) i1=1; 
				if(dihedralstruct[k].first==id2 || dihedralstruct[k].second==id2 || 
				   dihedralstruct[k].third==id2 || dihedralstruct[k].fourth==id2) i2=1; 
				if(i1 && i2) {iflag=1;break;}
			}
			if(iflag) continue;
			*/
				
			int s=vdwlist.size();
			vdwlist.resize(s+1);
			vdwlist[s].resize(2);
			vdwlist[s][0]=id1;
			vdwlist[s][1]=id2;
		}
	}

	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<list.size();i++)
		{
			Data totalEnergy=0.0;
			if(list[i]=="bond")
			{
				for(int j=0;j<bondstruct.size();j++)
				{
					if(bondFlag[j]==0) continue;
					int id1=bondstruct[j].first-1;
					int id2=bondstruct[j].second-1;
					string a=table[id1].atom_type;
					string b=table[id2].atom_type;
					Data e=pf.getBondEnergy(a,b,getDistance(id1,id2));
					totalEnergy+=e;
				}	
			}
			else
			if(list[i]=="angle")
			{
				for(int j=0;j<anglestruct.size();j++)
				{
					if(angleFlag[j]==0) continue;
					int id1=anglestruct[j].first-1;
					int id2=anglestruct[j].second-1;
					int id3=anglestruct[j].third-1;
					string a=table[id1].atom_type;
					string b=table[id2].atom_type;
					string c=table[id3].atom_type;
					Data e=pf.getAngleEnergy(a,b,c,getAngle(id1,id2,id3),
						getDistance(id1,id3));
					totalEnergy+=e;
				}
			}
			else
			if(list[i]=="dihed")
			{
				for(int j=0;j<dihedralstruct.size();j++)
				{
					if(diheFlag[j]==0) continue;
					int id1=dihedralstruct[j].first-1;
					int id2=dihedralstruct[j].second-1;
					int id3=dihedralstruct[j].third-1;
					int id4=dihedralstruct[j].fourth-1;
					string a=table[id1].atom_type;
					string b=table[id2].atom_type;
					string c=table[id3].atom_type;
					string d=table[id4].atom_type;
					Data e;
					e=pf.getDihedralEnergy(a,b,c,d,dihedral(id1,id2,id3,id4));
							     pf.getImproperEnergy(a,b,c,d,dihedral(id1,id2,id3,id4));
					totalEnergy+=e;
				}
			}
			else
			if(list[i]=="imprp")
			{
				for(int j=0;j<improperstruct.size();j++)
				{
					if(impFlag[j]==0) continue;
					int id1=improperstruct[j].first-1;
					int id2=improperstruct[j].second-1;
					int id3=improperstruct[j].third-1;
					int id4=improperstruct[j].fourth-1;
					string a=table[id1].atom_type;
					string b=table[id2].atom_type;
					string c=table[id3].atom_type;
					string d=table[id4].atom_type;
					Data e;
					e=pf.getImproperEnergy(a,b,c,d,dihedral(id1,id2,id3,id4));
					totalEnergy+=e;
				}
			}
			else
			if(list[i]=="elect")
			{
				for(int j=0;j<vdwlist.size();j++)
				{
					int id1=vdwlist[j][0]-1;
					int id2=vdwlist[j][1]-1;
					Data d=getDistance(id1,id2);
					if(d>12) continue;
					Data e=332.0636*table[id1].dnum1 * table[id2].dnum1/(diel * d);
					totalEnergy+=e;
				}
			}
			else
			if(list[i]=="vdw")
			{
				for(int j=0;j<vdwlist.size();j++)
				{
					int id1=vdwlist[j][0]-1;
					int id2=vdwlist[j][1]-1;
					string a=table[id1].atom_type;
					string b=table[id2].atom_type;
					Data d=getDistance(id1,id2);
					if(d>12) continue;
					Data e;
					e=pf.getVdwEnergy(a,b,d);
					totalEnergy+=e;
				}
			}
			E[i][icount]=totalEnergy;
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}


void	Dcd::getFitted(int fitflag,IntVector &AtomPos,PointVector &pdbpos,
				DoubleVector2 &Orig,DoubleVector2 &Dest)
{
	int N=0;
	IntVector okFlag;
	okFlag.resize(AtomPos.size());
	
	for(int i=0;i<AtomPos.size();i++)
	{
		okFlag[i]=0;
		switch(fitflag)
		{
			case HFIT:
				if(table[AtomPos[i]].atom_name[0]!='H') okFlag[i]=1;
				break;
			case BBFIT:
				if(isBackbone3(table[AtomPos[i]])) okFlag[i]=1;
				break;
			case CAFIT:
				if(table[AtomPos[i]].atom_name=="CA") okFlag[i]=1;
				break;
		}
	}
	for(int i=0;i<AtomPos.size();i++) N+=okFlag[i];
	Data *v1=new Data[N * 3];
	Data *v2=new Data[N * 3];
	Data *v3=new Data[N * 3];
	Data *mtx=new Data[16];
	Data centerB[3]={0,0,0};
	
	int icount=0;
	for(int i=0;i<AtomPos.size();i++)
	{
		if(okFlag[i]==0) continue;
		v2[icount*3+0]=pdbpos[AtomPos[i]].x;
		v2[icount*3+1]=pdbpos[AtomPos[i]].y;
		v2[icount*3+2]=pdbpos[AtomPos[i]].z;
		centerB[0]+=v2[icount*3+0];
		centerB[1]+=v2[icount*3+1];
		centerB[2]+=v2[icount*3+2];
		icount++;
	}
	centerB[0]/=N;centerB[1]/=N;centerB[2]/=N;
	Data distance=0;
	Data centerA[3]={0,0,0};
	icount=0;
	for(int i=0;i<AtomPos.size();i++)
	{
		if(okFlag[i]==0) continue;
		v1[icount*3+0]=Orig[AtomPos[i]][0];
		v1[icount*3+1]=Orig[AtomPos[i]][1];
		v1[icount*3+2]=Orig[AtomPos[i]][2];
		centerA[0]+=v1[icount*3+0];
		centerA[1]+=v1[icount*3+1];
		centerA[2]+=v1[icount*3+2];
	}
	centerA[0]/=N;centerA[1]/=N;centerA[2]/=N;
	
	for(int i=0;i<N;i++)
	{
		if(okFlag[i]==0) continue;
		v1[i*3+0]=v1[i*3+0]+(centerB[0]-centerA[0]);
		v1[i*3+1]=v1[i*3+1]+(centerB[1]-centerA[1]);
		v1[i*3+2]=v1[i*3+2]+(centerB[2]-centerA[2]);
	}
	double old_distance=1e+10;
	double diff=0.0;
	for(int iters=1;iters<=10;iters++)
	{
		distance=rmsd(v1, v2, N,mtx,v3);
		diff=fabs(distance-old_distance);
		if(diff<0.1) break;
		old_distance=distance;
	}
	delete[] v1;
	delete[] v2;
	if (mtx[3]==0 && mtx[7]==0 && mtx[11]==0 && mtx[15]==1) 
	{
  		const Data dx = mtx[12];
		const Data dy = mtx[13];
    		const Data dz = mtx[14];
       		for (int i=0; i<AtomPos.size(); i++)
		{
			Data pos[3];
			pos[0]=Orig[AtomPos[i]][0];
			pos[1]=Orig[AtomPos[i]][1];
			pos[2]=Orig[AtomPos[i]][2];
         		Dest[AtomPos[i]][0] = pos[0]*mtx[0] + pos[1]*mtx[4] + pos[2]*mtx[8] + dx;
        		Dest[AtomPos[i]][1] = pos[0]*mtx[1] + pos[1]*mtx[5] + pos[2]*mtx[9] + dy;
         		Dest[AtomPos[i]][2] = pos[0]*mtx[2] + pos[1]*mtx[6] + pos[2]*mtx[10] + dz;
		}
	}
	else
	{
		for(int i=0;i<AtomPos.size();i++)
		{
			Data pos[3];
			pos[0]=Orig[AtomPos[i]][0];
			pos[1]=Orig[AtomPos[i]][1];
			pos[2]=Orig[AtomPos[i]][2];
			Data tmp[3];
			Data itmp3 = 1.0f / (pos[0]*mtx[3] + pos[1]*mtx[7] +
                           pos[2]*mtx[11] + mtx[15]);
     			tmp[0] = itmp3*pos[0];
			tmp[1] = itmp3*pos[1];
			tmp[2] = itmp3*pos[2];
     			Dest[AtomPos[i]][0]=tmp[0]*mtx[0] + tmp[1]*mtx[4] + 
					    tmp[2]*mtx[ 8] + itmp3*mtx[12];
     			Dest[AtomPos[i]][1]=tmp[0]*mtx[1] + tmp[1]*mtx[5] + 
					     tmp[2]*mtx[ 9] + itmp3*mtx[13];
     			Dest[AtomPos[i]][2]=tmp[0]*mtx[2] + tmp[1]*mtx[6] + 
					    tmp[2]*mtx[10] + itmp3*mtx[14];
		}
	}
	/*
	for(int i=0;i<N;i++) 
	{
		Dest[AtomPos[i]][0]=v3[i*3+0];
		Dest[AtomPos[i]][1]=v3[i*3+1];
		Dest[AtomPos[i]][2]=v3[i*3+2];
	}
	*/
	delete[] v3;
	delete[] mtx;
}

void	Dcd::getDisu(IntVector &F,IntVector &Index1,IntVector &Index2,IntVector &Index3,
				IntVector &Index4,DoubleVector2 &D,DoubleVector2 &A)
{
	F.resize(nframes());
	D.resize(Index1.size());
	A.resize(Index1.size());
	for(int i=0;i<D.size();i++)
	{
		D[i].resize(F.size());
		A[i].resize(F.size());
	}
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		for(int i=0;i<Index1.size();i++)
		{
			int DISTFIRST=Index1[i]-1;
			int DISTSECOND=Index2[i]-1;
			int DISTTHIRD=Index3[i]-1;
			int DISTFOURTH=Index4[i]-1;
			D[i][icount]=getDistance(DISTFIRST,DISTFOURTH);	
			A[i][icount]=dihedral(DISTFIRST,DISTSECOND,DISTTHIRD,DISTFOURTH);
		}
		F[icount]=frame;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}

void	Dcd::getSasa(IntVector &F,IntVector &posAtom,Data srad,DoubleVector &D)
{
	F.resize(nframes());
	D.resize(F.size());
	int FIRST = frameStart;
	int istart =lseek(dcd,(off_t)(0),SEEK_CUR);
	int frame= FIRST;
	int  firstframe = FIRST;
	int itotalframes=0;
	int icount=0;

	const int npts=100;
	DoubleVector spherepts;
	spherepts.resize(3*npts);

	Data maxrad=-1.0;
	extern vector<atom> table;
	for(int i=0;i<posAtom.size();i++)
	{
		if(table[posAtom[i]].radious>maxrad)
			maxrad=table[posAtom[i]].radious;
	}
	while(frame<=frameEnd)
	{
	 if((frame-firstframe)%frameStep==0)
	 {
		if(!readCAs()) break;
		IntVector2 pair;
		pair.resize(posAtom.size());
		for(int i=0;i<posAtom.size();i++)
		{
			for(int j=i+1;j<posAtom.size();j++)
			{
				Data d=getDistance(posAtom[i],posAtom[j]);
				if(d<=2.0*(maxrad+srad))
				{
					pair[i].push_back(j);
					pair[j].push_back(i);
				}
			}
		}
		Data RAND_MAX_INV=1.0/RAND_MAX;
		srand48(38572111);
		for(int i=0;i<npts;i++)
		{
			Data u1=drand48();
			Data u2=drand48();
			Data z=2*u1*RAND_MAX_INV-1.0;
			Data phi=(2.0*M_PI*u2*RAND_MAX_INV);
			Data R=sqrt(1.0-z*z);
			spherepts[3*i  ] = R*cosf(phi);
			spherepts[3*i+1] = R*sinf(phi);
			spherepts[3*i+2] = z;
		}
		const Data prefac = (float) (4 * M_PI / npts);
		Data totarea = 0.0f;
		for(int i=0;i<posAtom.size();i++)
		{
			Data rad=table[posAtom[i]].radious +srad;
			DoubleVector surfpos;surfpos.resize(3);
			int surfpts=npts;
			for(int j=0;j<npts;j++)
			{
				surfpos[0]=CAs[posAtom[i]][0]+
					rad*spherepts[3*j+0];
				surfpos[1]=CAs[posAtom[i]][1]+
					rad*spherepts[3*j+1];
				surfpos[2]=CAs[posAtom[i]][2]+
					rad*spherepts[3*j+2];
				int on=1;
				for(int k=0;k<pair[i].size();k++)
				{
					int ind=posAtom[pair[i][k]];
					Data radsq = table[ind].radious+srad; 
					radsq *= radsq;
					Data dx = surfpos[0]-CAs[ind][0];
					Data dy = surfpos[1]-CAs[ind][1];
					Data dz = surfpos[2]-CAs[ind][2];
					if (dx*dx + dy*dy + dz*dz <= radsq) 
					{
						on = 0;
						break;
					}
				}
				if(!on) surfpts--;
			}
			float atomarea = prefac * rad * rad * surfpts;
			totarea += atomarea;
		}
		F[icount]=frame;
		D[icount]=totarea;
	 	icount++;
	 	frame++;
	 }
	 else if(!skipFrame(frame)) break;
	}
	lseek(dcd,(off_t)istart,SEEK_SET);
}
