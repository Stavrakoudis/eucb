# include <psf.h>
# include <util.h>
# include <stdlib.h>


vector<BondStruct> bondstruct;
vector<AngleStruct> anglestruct;
vector<DihedralStruct> dihedralstruct;
vector<DihedralStruct> improperstruct;

Data	getradious(string res_name,string atom_name)
{
	for(int i=0;i<rad.size();i++)
	{
		if(rad[i].res_name == res_name
		 ||(rad[i].res_name == "ASX" &&res_name=="ASP")	
		 ||(rad[i].res_name == "ASX" &&res_name=="ASN")
		 ||(rad[i].res_name=="GLX" &&res_name=="GLU")
		 ||(rad[i].res_name=="GLX" &&res_name=="GLN")
		 ||(rad[i].res_name=="HIS" &&res_name=="HSE")
		 ||(rad[i].res_name=="HIS" &&res_name=="HSD")
		 ||(rad[i].res_name=="HIS" &&res_name=="HSP")
		   )
		{
			if(rad[i].atom_name==atom_name 
			||(rad[i].atom_name=="AD1" && atom_name=="ND1")
			||(rad[i].atom_name=="AD1" && atom_name=="OD1")
			||(rad[i].atom_name=="AD2" && atom_name=="ND2")
			||(rad[i].atom_name=="AD2" && atom_name=="OD2")
			||(rad[i].atom_name=="AE1" && atom_name=="OE1")
			||(rad[i].atom_name=="AE1" && atom_name=="NE1")
			||(rad[i].atom_name=="AE2" && atom_name=="OE2")
			||(rad[i].atom_name=="AE2" && atom_name=="NE2")
			)
			{
				return rad[i].rad;
			}
		}
	}
	return -1.0;
}

string	printAtomLow(atom a)
{
	string s="";
	s=s+a.res_name[0];
	for(int i=1;i<a.res_name.size();i++)
	{
		char x=a.res_name[i];
		x=x+'a'-'A';
		s=s+x;
	}
	return s;
}

string	printAtomGreek(atom a,int isgroup)
{
	string s="";
	s=a.atom_name[0];
	s=s+"$^{";
	if(a.atom_name.size()>1)
	{
		char x=a.atom_name[1];
		if(x=='H') s=s+"\\eta";
		else
		if(x=='G') s=s+"\\gamma";
		else
		if(x=='A') s=s+"\\alpha";
		else
		if(x=='B') s=s+"\\beta";
		else
		if(x=='D') s=s+"\\delta";
		else
		if(x=='E') s=s+"\\epsilon";
		else
		if(x=='Z') s=s+"\\zeta";
		else
		if(x=='T') s=s+"\\tau";
		else
		s=s+x;
	}
	if(a.atom_name.size()>2) s=s+a.atom_name[2];
	if(a.atom_name.size()>3) s=s+a.atom_name[3];
	if(isgroup) s=s+"2";
	s=s+"}$";
	return s;
}

string 	printAtomLatex(atom a,int f)
{
	string s="";
	char xx[100];
	sprintf(xx,"%d",a.res_id);
	if(f)
	s=printAtomLow(a)+"$_{"+xx+a.chain_id+"}$:"+printAtomGreek(a,0);
	else
	s=printAtomLow(a)+"$_{"+xx+a.chain_id+"}$";
	sprintf(xx,"%28s",s.c_str());
	s=xx;
	return s;
}

string 	printAtomLatexGroup(atom a,int f)
{
	string s="";
	char xx[100];
	sprintf(xx,"%d",a.res_id);
	if(f)
	s=printAtomLow(a)+"$_{"+xx+a.chain_id+"}$:"+printAtomGreek(a,1);
	else
	s=printAtomLow(a)+"$_{"+xx+a.chain_id+"}$";
	sprintf(xx,"%28s",s.c_str());
	s=xx;
	return s;
}

string	printAtom1(atom a)
{
	string st=a.chain_id+"_"+printNumber(a.res_id)+"_"+a.res_name+"_"+a.atom_name;
	return st;
}

string	printAtom3(atom a)
{
	return printString(a.chain_id,5)+" "+printNumber(a.res_id,5)+" "+
		printString(a.res_name,4)+" "+
		printString(a.atom_name,4);
}

string 	printAtom4(atom a)
{
	string st= a.chain_id+":"+printNumber(a.res_id)+":"+
		a.res_name+":"+a.atom_name;
	return printString(st,-15);
	/*	
	char s[1024];
	sprintf(s,"%s:%d:%s:%s",a.chain_id.c_str(),a.res_id,
		a.res_name.c_str(),a.atom_name.c_str());
	char s1[100];
	sprintf(s1,"%-15s",s);
	return s1;
	*/
}

string	printAtom5(atom a)
{
	string st=a.chain_id+":"+printNumber(a.res_id)+":"+a.res_name;
	return st;
}

string 	printAtom6(atom a)
{
	string st=a.res_name+":"+printNumber(a.res_id)+":"+a.atom_name;
	return st;
}

string	printAtomHeader()
{
	return printString("CHAIN",5)+" "+printString("RESID",5)+" "+
		printString("RES",4)+" "+
		printString("ATOM",4);
}

string	printAtom2(atom a)
{
	string st=a.chain_id+"_"+printNumber(a.res_id)+"_"+a.res_name;
	return st;
}


/*	Store on the vector result all the chains in the psf file.	
 * */
void	enumerateChains(vector<atom> &table,vector<string> &result)
{
	result.resize(0);
	for(int i=0;i<table.size();i++)
	{
		string s=table[i].chain_id;
		int found=0;
		for(int j=0;j<result.size();j++) 
			if(result[j]==s) {found=1;break;}
		if(!found) result.push_back(s);
	}
}


int	countConnected(vector<atom> &table,int id)
{
	return table[id-1].bond.size();
}

/*	Return 1 if the atom with id id1 is connected with the atom with id id2.
 *	Otherwise return 0.
 * */
int 	isConnected(vector<atom> &table,int id1,int id2)
{
	for(int i=0;i<table[id1-1].bond.size();i++)
		if(table[id1-1].bond[i]==id2) return 1;
	return 0;
}

void	write_psf(string fname,vector<atom> &table,IntVector &posAtom)
{
	FILE *fp;
	char line[1024];
	fp=fopen(fname.c_str(),"w");
	if(!fp) 
	{
		string s1="PSF FILE "+fname+" NOT READY FOR OPENING";
		psf_error(s1);
	}
	strcpy(line,"PSF");
	fprintf(fp,"%s\n",line);
	int natoms=posAtom.size();
	sprintf(line,"%8d !NATOM ",natoms);
	fprintf(fp,"%s\n",line);
	int totalbonds=0;
	for(int i=0;i<natoms;i++)
	{
		
		sprintf(line,"%8d %-2s %3d %6s  %-4s %-4s %10.6f %13.4f %11d",
				i+1,
				table[posAtom[i]].chain_id.c_str(),
				table[posAtom[i]].res_id,
				table[posAtom[i]].res_name.c_str(),
				table[posAtom[i]].atom_name.c_str(),
				table[posAtom[i]].atom_type.c_str(),
				table[posAtom[i]].dnum1,
				table[posAtom[i]].dnum2,
				table[posAtom[i]].lastnum);
		fprintf(fp,"%s\n",line);
		
	}
	for(int i=0;i<bondstruct.size();i++)
	{
		if(isin(posAtom,bondstruct[i].first-1)!=-1 && 
		   isin(posAtom,bondstruct[i].second-1)!=-1) totalbonds++;
	}
	sprintf(line,"\n%8d !NBOND",totalbonds);
	fprintf(fp,"%s\n",line);
	int pair=0;
	for(int i=0;i<bondstruct.size();i++)
	{
		int id1=isin(posAtom,bondstruct[i].first-1);
		int id2=isin(posAtom,bondstruct[i].second-1);
		if(id1!=-1 && id2!=-1) 
			//fprintf(fp,"\t%d\t%d ",id1+1,id2+1);
			fprintf(fp," %7d %7d",id1+1,id2+1);
		else	continue;
		pair++;
		if(pair%5==0) fprintf(fp,"\n");
	}
	fclose(fp);
}

/*	The function parses the file psffile and stores the psf information into
 *	the vector table.
 * */
int	parse_psf(char *psffile,vector<atom> &table)
{
	makeRad();
	FILE *fp;
	char line[1024];
	fp=fopen(psffile,"r");
	if(!fp) 
	{
		string s2=psffile;
		string s1="PSF FILE "+s2+" NOT READY FOR OPENING";
		psf_error(s1);
	}
	fgets(line,1023,fp);
	line[strlen(line)-1]=0;
	vector<string> list;
	getAllWords(line,list);
	if(list[0][0]!='P' || list[0][1]!='S' || list[0][2]!='F')
	{
		fclose(fp);
		string s2=psffile;
		string s1="FILE "+s2+" IS NOT A VALID PSF FILE ";
		psf_error(s1);
	}
	while(1)
	{
		list.resize(0);
		if(!fgets(line,1023,fp)) break;
		line[strlen(line)-1]=0;
		getAllWords(line,list);
		int natoms;
		if(list[1]=="!NATOM")
		{
			natoms = atoi(list[0].c_str());
			table.resize(natoms);
			for(int i=0;i<natoms;i++)
			{
				list.resize(0);
				if(!fgets(line,1023,fp)) break;
				line[strlen(line)-1]=0;
				getAllWords(line,list);
				table[i].atom_id=atoi(list[0].c_str());
                                table[i].chain_id=list[1];
                                table[i].res_id=atoi(list[2].c_str());
                                table[i].res_name=list[3];
                                table[i].atom_name=list[4];
                                table[i].atom_type=list[5];
                                table[i].dnum1=atof(list[6].c_str());
                                table[i].dnum2=atof(list[7].c_str());
                                table[i].lastnum=atoi(list[8].c_str());
				table[i].ishydro=0;
				table[i].radious=getradious(
					table[i].res_name,table[i].atom_name);
			}
		}
		else
		if(list[1]=="!NBOND:")
		{
			BondStruct t1;
			int first,second;
			int nbonds = atoi(list[0].c_str());
			for(int j=0;j<2*nbonds;j++)
			{
				int k;
				fscanf(fp,"%d",&k);
				if(j%2==0) 
				{
					first = k;
				}
				else
				{
					second = k;
					
					table[first-1].bond.push_back(second);
					table[second-1].bond.push_back(first);

					t1.first=first;
					t1.second=second;
					bondstruct.push_back(t1);
				}
			}
		}
		else
		if(list[1]=="!NTHETA:")
		{
			AngleStruct t1;
			int first,second,third;
			int nbonds = atoi(list[0].c_str());
			for(int j=0;j<3*nbonds;j++)
			{
				int k;
				fscanf(fp,"%d",&k);
				if(j%3==0) 
				{
					first = k;
				}
				else
				if(j%3==1)
				{
					second = k;
				}
				else
				if(j%3==2)
				{
					third=k;
					t1.first=first;
					t1.second=second;
					t1.third=third;
					anglestruct.push_back(t1);
				}
			}
		}
		else
		if(list[1]=="!NPHI:")
		{
			DihedralStruct t1;
			int first,second,third,fourth;
			int nbonds = atoi(list[0].c_str());
			for(int j=0;j<4*nbonds;j++)
			{
				int k;
				fscanf(fp,"%d",&k);
				if(j%4==0) 
				{
					first = k;
				}
				else
				if(j%4==1)
				{
					second = k;
				}
				else
				if(j%4==2)
				{
					third=k;
				}
				else
				if(j%4==3)
				{
					fourth=k;
					t1.first=first;
					t1.second=second;
					t1.third=third;
					t1.fourth=fourth;
					dihedralstruct.push_back(t1);
				}
			}
		}
		else
		if(list[1]=="!NIMPHI:")
		{
			DihedralStruct t1;
			int first,second,third,fourth;
			int nbonds = atoi(list[0].c_str());
			for(int j=0;j<4*nbonds;j++)
			{
				int k;
				fscanf(fp,"%d",&k);
				if(j%4==0) 
				{
					first = k;
				}
				else
				if(j%4==1)
				{
					second = k;
				}
				else
				if(j%4==2)
				{
					third=k;
				}
				else
				if(j%4==3)
				{
					fourth=k;
					t1.first=first;
					t1.second=second;
					t1.third=third;
					t1.fourth=fourth;
					improperstruct.push_back(t1);
				}
			}
		}
	}
	fclose(fp);
	return 1;
}
int	isHydro(atom a)
{
	return a.atom_name[0]=='H';
}

int	isBackbone(atom a)
{
	return isBackbone4(a);
}

int	isBackbone4(atom a)
{
	return isBackbone3(a) || a.atom_name == "O";
}

int 	isBackbone3(atom a)
{
	return a.atom_name == "N" || a.atom_name=="CA" || a.atom_name=="C";
}

int	isSidechain(atom a)
{
	return a.atom_name[0]!='H' && !isBackbone(a);
}

int	isWater(atom a)
{
	return a.res_name == "TIP3";
}

int	isPositive(atom a)
{
	if(a.res_name=="LYS" && a.atom_name=="NZ") return 1;
	else
	if(a.res_name=="ARG")
	{
		if(a.atom_name == "NE" || a.atom_name == "NH1" || a.atom_name=="NH2") return 1;
	}
	else
	if(a.res_name == "HSP")
	{
		if(a.atom_name == "ND1" || a.atom_name=="NE2") return 1;
	}
	return 0;
}

int	isNegative(atom a)
{
	if(a.res_name == "ASP")
	{
		if(a.atom_name=="OD1" || a.atom_name=="OD2") return 1;
	}	
	else
	if(a.res_name == "GLU")
	{
		if(a.atom_name == "OE1" || a.atom_name == "OE2") return 1;
	}
	return 0;
}

int 	isHydroPhobic(atom a)
{
	if(a.res_name == "ALA" || a.res_name=="ILE" || a.res_name=="LEU" 
		|| a.res_name=="PRO" || a.res_name=="VAL" || a.res_name=="TYR" 
		|| a.res_name=="TRP" || a.res_name=="PHE")
	{
		if(a.atom_name[0]!='H') return 1;
	}
	return 0;
}
