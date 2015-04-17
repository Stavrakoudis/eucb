# include <parfile.h>
# include <math.h>
/*
typedef	struct
{
	string	atomtype1;
	string	atomtype2;
	string	atomtype3;
	string	atomtype4;
	Data	kpsi;
	Data	psi0;
}ParImproperStruct;
*/


ParFile::ParFile(string s)
{
	FILE *fp=fopen(s.c_str(),"r");
	if(!fp) psf_error("Unable to open "+s+" for reading");
	char line[1024];	
	while(!feof(fp))
	{
		fgets(line,1023,fp);
		vector<string> list;
		line[strlen(line)-1]=0;
		getAllWords(line,list);
		if(strlen(line)==0) continue;
		if(list[0][0]=='!') continue;//comment line
		if(list[0]=="BONDS")
		{
			while(1)	
			{
				list.resize(0);
				fgets(line,1023,fp);
				line[strlen(line)-1]=0;
				if(strlen(line)==0) break;
				getAllWords(line,list);
				if(list[0][0]=='!') continue;
				ParBondStruct st;
				st.atomtype1=list[0];
				if(list[0].find_first_of('!')!=string::npos) continue;
				st.atomtype2=list[1];
				st.kb=atof(list[2].c_str());
				st.b0=atof(list[3].c_str());
				bondstruct.push_back(st);
			}
		}
		
		else
		if(list[0]=="ANGLES")
		{
			while(1)	
			{
				list.resize(0);
				fgets(line,1023,fp);
				line[strlen(line)-1]=0;
				if(strlen(line)==0) break;
				getAllWords(line,list);
				if(list[0][0]=='!') continue;
				ParAngleStruct st;
				st.atomtype1=list[0];
				if(list[0].find_first_of('!')!=string::npos) continue;
				st.atomtype2=list[1];
				st.atomtype3=list[2];
				st.ktheta=atof(list[3].c_str());
				st.theta0=atof(list[4].c_str());
				if(list[5][0]!='!' && isdigit(list[5][0])) 
				{
					st.kub=atof(list[5].c_str());
					st.s0=atof(list[6].c_str());
				}
				else
				{
					st.kub=0.0;
					st.s0=0.0;
				}
				anglestruct.push_back(st);
			}
		}
		else
		if(list[0]=="DIHEDRALS")
		{
			while(1)	
			{
				list.resize(0);
				fgets(line,1023,fp);
				line[strlen(line)-1]=0;
				if(strlen(line)==0) break;
				getAllWords(line,list);
				if(list[0][0]=='!') continue;
				ParDihedralStruct st;
				st.atomtype1=list[0];
				if(list[0].find_first_of('!')!=string::npos) continue;
				st.atomtype2=list[1];
				st.atomtype3=list[2];
				st.atomtype4=list[3];
				st.kchi=atof(list[4].c_str());
				st.n=atoi(list[5].c_str());
				st.delta=atof(list[6].c_str());	
				dihedralstruct.push_back(st);
			}
		}
		else
		if(list[0]=="IMPROPER")
		{
			while(1)	
			{
				list.resize(0);
				fgets(line,1023,fp);
				line[strlen(line)-1]=0;
				if(strlen(line)==0) break;
				getAllWords(line,list);
				if(list[0][0]=='!') continue;
				ParImproperStruct st;
				st.atomtype1=list[0];
				if(list[0].find_first_of('!')!=string::npos) continue;
				st.atomtype2=list[1];
				st.atomtype3=list[2];
				st.atomtype4=list[3];
				st.kpsi=atof(list[4].c_str());
				st.psi0=atof(list[6].c_str());	
				improperstruct.push_back(st);
			}
		}
		else
		if(list[0]=="NONBONDED")
		{
			while(1)	
			{
				list.resize(0);
				char *d=fgets(line,1023,fp);
				if(!d) break;
				line[strlen(line)-1]=0;
				if(strlen(line)==0) continue;
				if(list[0]=="HBOND") break;
				getAllWords(line,list);
				if(list[0][0]=='!') continue;
				ParVdwStruct st;
				st.atomtype1=list[0];
				if(list[0].find_first_of('!')!=string::npos) continue;
				st.epsilon=atof(list[2].c_str());
				st.rmin=atof(list[3].c_str());
				if(st.atomtype1=="HBOND" || st.atomtype1=="END") break;
				vdwstruct.push_back(st);
			}
		}
		
	}
	fclose(fp);
}

Data	ParFile::getBondEnergy(string a,string b,Data B)
{
	for(int i=0;i<bondstruct.size();i++)
	{
		if((bondstruct[i].atomtype1==a && bondstruct[i].atomtype2==b) 
		  || (bondstruct[i].atomtype2==a && bondstruct[i].atomtype1==b)) 
		{
			Data e=bondstruct[i].kb*((B-bondstruct[i].b0)*(B-bondstruct[i].b0));
			return e;
		}
	}
	return 0.0;
}

Data	ParFile::getAngleEnergy(string a,string b,string c,Data theta,Data d)
{
	for(int i=0;i<anglestruct.size();i++)
	{
		if((anglestruct[i].atomtype1==a && anglestruct[i].atomtype2==b &&
		    anglestruct[i].atomtype3==c) ||
                   (anglestruct[i].atomtype1==c && anglestruct[i].atomtype2==b &&
		    anglestruct[i].atomtype3==a))
		{
			Data e=0.0;
			e=anglestruct[i].ktheta*(pow(deg2rad(theta)-deg2rad(anglestruct[i].theta0),2.0));
			if(fabs(anglestruct[i].kub)>1e-7 && fabs(anglestruct[i].s0)>1e-7)
			{
				e+=anglestruct[i].kub*(pow(d-anglestruct[i].s0,2.0));
			}
			return e;
		}
	}
	return 0.0;
}

Data	ParFile::getDihedralEnergy(string a,string b,string c,string d,Data chi)
{
	Data e=0.0;
	int first=0;
	for(int i=0;i<dihedralstruct.size();i++)
	{
		if(((dihedralstruct[i].atomtype1==a||dihedralstruct[i].atomtype1=="X") 
		     && dihedralstruct[i].atomtype2==b &&
		    dihedralstruct[i].atomtype3==c && 
		   (dihedralstruct[i].atomtype4==d||dihedralstruct[i].atomtype4=="X")) ||

		   ((dihedralstruct[i].atomtype1==d||dihedralstruct[i].atomtype1=="X") 
			 && dihedralstruct[i].atomtype2==c &&
		    dihedralstruct[i].atomtype3==b && 
			(dihedralstruct[i].atomtype4==a||dihedralstruct[i].atomtype4=="X")))
		{
			if(first==0) first=1;
			if(dihedralstruct[i].n)
			{
				e+=dihedralstruct[i].kchi*
				    ((1.0+cos(dihedralstruct[i].n*deg2rad(chi)-deg2rad(dihedralstruct[i].delta))));
			}
			else
			{
				 Data diff = deg2rad(chi)-deg2rad(dihedralstruct[i].delta);
      				if (diff < -M_PI)           diff += 2*M_PI;
      				else if (diff > M_PI)       diff -= 2*M_PI;

      				e+= dihedralstruct[i].kchi*diff*diff;
	
			}
		}
		else
		{
			if(first) return e;
		}
	}
			return e;
	return 0.0;
}

Data	ParFile::getImproperEnergy(string a,string b,string c,string d,Data psi)
{
	Data e=0.0;
	int first=0;
	for(int i=0;i<improperstruct.size();i++)
	{
		if(((improperstruct[i].atomtype1==a) 
		     && (improperstruct[i].atomtype2==b||improperstruct[i].atomtype2=="X") &&
		    (improperstruct[i].atomtype3==c || improperstruct[i].atomtype3=="X") && 
		   (improperstruct[i].atomtype4==d)) 
			||
		   ((improperstruct[i].atomtype1==d) 
			 && (improperstruct[i].atomtype2==c||improperstruct[i].atomtype2=="X") &&
		    	(improperstruct[i].atomtype3==b || improperstruct[i].atomtype3=="X") && 
			(improperstruct[i].atomtype4==a)))
		{
			if(first==0) first=1;
			 Data diff = deg2rad(psi)-deg2rad(improperstruct[i].psi0);
      			if (diff < -M_PI)           diff += 2*M_PI;
      			else if (diff > M_PI)       diff -= 2*M_PI;
			e+= improperstruct[i].kpsi*diff*diff;
		}
		else
		{
			if(first) return e;
		}
	}
			return e;
	return 0.0;
}

Data	ParFile::getVdwEnergy(string a,string b,Data d)
{
	int index1=-1,index2=-1;
	for(int i=0;i<vdwstruct.size();i++)
	{
		if(vdwstruct[i].atomtype1==a) index1=i;
		if(vdwstruct[i].atomtype1==b) index2=i;
		if(index1!=-1 && index2!=-1) break;
	}
	if(index1==-1 || index2==-1) return 0.0;
	Data e=0.0;
	Data eps=sqrt(vdwstruct[index1].epsilon * vdwstruct[index2].epsilon);
	Data sig=(vdwstruct[index1].rmin+vdwstruct[index2].rmin);
	e=eps * (pow(sig/d,12)-2*pow(sig/d,6));
	return e;
}

ParFile::~ParFile()
{
}
