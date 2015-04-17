# include <global.h>
# include <math.h>
# include <command_rdf.h>

CommandRdf::CommandRdf(int id,string criteria)
	:Command("RDF")
{
	rdf_id = id;
	rdf_criteria = criteria;
	radious = 0.0;
	radious_step=0.25;
	radious_max = 10.0;	
	volume = -1.0;
	xst_filename = "";
	xst_freq=500.0;
}

void	CommandRdf::setRadious(Data r1,Data r2,Data r3)
{
	radious =r1;
	radious_step=r2;
	radious_max =r3;	
}

void	CommandRdf::setVolume(Data v)
{
	volume = v;
}

void	CommandRdf::setXstFilename(string s)
{
	xst_filename = s;
}


void	CommandRdf::setXstFreq(Data x)
{
	xst_freq = x;
}

void	CommandRdf::Run()
{
	//endiaferomaste monon gia nera
	vector<int> list;
	for(int i=0;i<table.size();i++)
	{
		if(table[i].res_name=="TIP3") continue;
		if(rdf_criteria!=table[i].atom_name)  continue;
		int s=list.size();
		list.resize(s+1);
		list[s]=table[i].atom_id;
	}
	Data dold=1e+100;
	Data V1,V2;
	Data rx;
	vector<Data> volume_per_frame;
	volume_per_frame.resize(dcd->getframes());

	if(volume<0)
	{
		FILE *fp=fopen(xst_filename.c_str(),"r");
		if(!fp) psf_error("YOU MUST SUPPLY A VOLUME OR A XST FILE");
		char line[1024];
		char word[1024];
		int lastframe=0;
		Data lastvolume;
		while(!feof(fp))
		{
			int iframe;
			Data r1,r2,r3;
			if(fgets(line,1024,fp)==NULL) break;
			line[strlen(line)-1]=0;
			if(line[0]=='#') continue;
			int index=0;
			getword(line,word,index);
			iframe=atoi(word)/xst_freq;
			if(iframe) iframe--;
			getword(line,word,index);	
			r1=atof(word);

			getword(line,word,index);	
			getword(line,word,index);	
			getword(line,word,index);	
			
			getword(line,word,index);	
			r2=atof(word);
			
			getword(line,word,index);	
			getword(line,word,index);	
			getword(line,word,index);	

			getword(line,word,index);	
			r3=atof(word);
			for(int i=lastframe;i<=iframe;i++)
			{
				volume_per_frame[i]=r1 * r2 * r3;
			}
			lastframe=iframe;
			lastvolume=r1*r2*r3;
		}
		
		for(int i=lastframe;i<volume_per_frame.size();i++)
			volume_per_frame[i]=lastvolume;
		fclose(fp);
	}
	else
	{
		for(int i=0;i<volume_per_frame.size();i++) volume_per_frame[i]=volume;
	}

	char filename[1024];
	sprintf(filename,"rdf_%d_%s_%s_%d.dat",table[rdf_id-1].res_id,
			table[rdf_id-1].res_name.c_str(),table[rdf_id-1].atom_name.c_str(),rdf_id);
	FILE *fp=fopen(filename,"w");
	for(Data r=radious;r<=radious_max;r+=radious_step)
	{
		V2=4.0 *M_PI*r*r*r;
		rx=list.size()*1.0/V2;
		rx=0.048;
		Data DV=V2-V1;
		int icount=0;
		Data dcount=dcd->countAtoms(rdf_id,r,list,volume_per_frame);
		if(dcount>=dold)
		{
		fprintf(fp,"%6.4lf %12.8lf %12.8lf \n",r,dcount/(DV),(dcount-dold)/(DV));
		printf("%6.4lf %12.8lf %12.8lf \n",r,dcount/(DV),(dcount-dold)/(DV));
		}
		else
		{
			DV=V2;
		}
		dold=dcount;
		V1=V2;
	}
	fclose(fp);
}


CommandRdf::~CommandRdf()
{
}
