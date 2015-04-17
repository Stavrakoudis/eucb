# include <team.h>
# define PEP_FLAG	0x1
# define NUM_FLAG 	0x2
# define ARG1_FLAG	0x4
# define ARG2_FLAG	0x8
# define WILDCARD	'_'
# include <string>
using namespace std;

/*	This class is used to create the 
 *	search sequences of the program.
 * */
		

Team	*makeHydrophobic()
{
	vector<string> arg1;
	vector<string> arg2;
	arg1.resize(35);
	arg2.resize(35);
	arg1[0]="PHE";arg2[0]="CG";
	arg1[1]="PHE";arg2[1]="CD1";
	arg1[2]="PHE";arg2[2]="CD2";
	arg1[3]="PHE";arg2[3]="CZ";

	arg1[4]="TYR";arg2[4]="CG";
	arg1[5]="TYR";arg2[5]="CD1";
	arg1[6]="TYR";arg2[6]="CD2";
	arg1[7]="TYR";arg2[7]="CZ";

	arg1[8]="TRP";arg2[8]="CG";
	arg1[9]="TRP";arg2[9]="CD1";
	arg1[10]="TRP";arg2[10]="NE1";
	arg1[11]="TRP";arg2[11]="CE2";
	arg1[12]="TRP";arg2[12]="CD2";
	arg1[13]="TRP";arg2[13]="CE3";
	arg1[14]="TRP";arg2[14]="CZ3";
	arg1[15]="TRP";arg2[15]="CZ2";
	arg1[16]="TRP";arg2[16]="CH2";

	arg1[17]="ALA";arg2[17]="CB";

	arg1[18]="LEU";arg2[18]="CB";
	arg1[19]="LEU";arg2[19]="CG";
	arg1[20]="LEU";arg2[20]="CD1";
	arg1[21]="LEU";arg2[21]="CD2";

	arg1[22]="ILE";arg2[22]="CB";
	arg1[23]="ILE";arg2[23]="CG1";
	arg1[24]="ILE";arg2[25]="CG2";
	arg1[25]="ILE";arg2[25]="CD";

	arg1[26]="PRO";arg2[26]="CB";
	arg1[27]="PRO";arg2[27]="CG";
	arg1[28]="PRO";arg2[28]="CD";


	arg1[29]="VAL";arg2[29]="CB";
	arg1[30]="VAL";arg2[30]="CG1";
	arg1[31]="VAL";arg2[31]="CG2";

	arg1[32]="MET";arg2[32]="CB";
	arg1[33]="MET";arg2[33]="SG";
	arg1[34]="MET";arg2[34]="CD";

	return new Team(arg1,arg2);
}

Team	*makeAromatic()
{
	vector<string> arg1;
	vector<string> arg2;
	arg1.resize(32);
	arg2.resize(32);
	arg1[0]="PHE";arg2[0]="CG";
	arg1[1]="PHE";arg2[1]="CD1";
	arg1[2]="PHE";arg2[2]="CD2";
	arg1[3]="PHE";arg2[3]="CZ";

	arg1[4]="TYR";arg2[4]="CG";
	arg1[5]="TYR";arg2[5]="CD1";
	arg1[6]="TYR";arg2[6]="CD2";
	arg1[7]="TYR";arg2[7]="CZ";

	arg1[8]="TRP";arg2[8]="CG";
	arg1[9]="TRP";arg2[9]="CD1";
	arg1[10]="TRP";arg2[10]="NE1";
	arg1[11]="TRP";arg2[11]="CE2";
	arg1[12]="TRP";arg2[12]="CD2";
	arg1[13]="TRP";arg2[13]="CE3";
	arg1[14]="TRP";arg2[14]="CZ3";
	arg1[15]="TRP";arg2[15]="CZ2";
	arg1[16]="TRP";arg2[16]="CH2";

	arg1[17]="HSP";arg2[17]="ND1";
	arg1[18]="HSP";arg2[18]="CG";
	arg1[19]="HSP";arg2[19]="CE1";
	arg1[20]="HSP";arg2[20]="NE2";
	arg1[21]="HSP";arg2[21]="CD2";

	arg1[22]="HSD";arg2[22]="ND1";
	arg1[23]="HSD";arg2[23]="CG";
	arg1[24]="HSD";arg2[24]="CE1";
	arg1[25]="HSD";arg2[25]="NE2";
	arg1[26]="HSD";arg2[26]="CD2";

	arg1[27]="HSE";arg2[27]="ND1";
	arg1[28]="HSE";arg2[28]="CG";
	arg1[29]="HSE";arg2[29]="CE1";
	arg1[30]="HSE";arg2[30]="NE2";
	arg1[31]="HSE";arg2[31]="CD2";
	return new Team(arg1,arg2);

}

Team	*makeAliphatic()
{
	vector<string> arg1;
	vector<string> arg2;
	arg1.resize(18);
	arg2.resize(18);
	arg1[0]="ALA";arg2[0]="CB";

	arg1[1]="LEU";arg2[1]="CB";
	arg1[2]="LEU";arg2[2]="CG";
	arg1[3]="LEU";arg2[3]="CD1";
	arg1[4]="LEU";arg2[4]="CD2";

	arg1[5]="ILE";arg2[5]="CB";
	arg1[6]="ILE";arg2[6]="CG1";
	arg1[7]="ILE";arg2[7]="CG2";
	arg1[8]="ILE";arg2[8]="CD";

	arg1[9]="PRO";arg2[9]="CB";
	arg1[10]="PRO";arg2[10]="CG";
	arg1[11]="PRO";arg2[11]="CD";


	arg1[12]="VAL";arg2[12]="CB";
	arg1[13]="VAL";arg2[13]="CG1";
	arg1[14]="VAL";arg2[14]="CG2";

	arg1[15]="MET";arg2[15]="CB";
	arg1[16]="MET";arg2[16]="SG";
	arg1[17]="MET";arg2[17]="CD";

	return new Team(arg1,arg2);

}

Team	*makePositive()
{
	vector<string> arg1;
	vector<string> arg2;
	arg1.resize(6);
	arg2.resize(6);
	arg1[0]="ARG";arg2[0]="NE";
	arg1[1]="ARG";arg2[1]="NH1";
	arg1[2]="ARG";arg2[2]="NH2";
	arg1[3]="LYS";arg2[3]="NZ";
	arg1[4]="HSP";arg2[4]="ND1";
	arg1[5]="HSP";arg2[5]="NE2";
	return new Team(arg1,arg2);
}

Team	*makeNegative()
{
	vector<string> arg1;
	vector<string> arg2;
	arg1.resize(4);
	arg2.resize(4);
	arg1[0]="ASP";arg2[0]="OD1";
	arg1[1]="ASP";arg2[1]="OD2";
	arg1[2]="GLU";arg2[2]="OE1";
	arg1[3]="GLU";arg2[3]="OE2";
	return new Team(arg1,arg2);
}

void	Team::enumerateAtoms(vector<atom> &list,vector<int> &res)
{
	res.resize(0);
	string chain;
	for(int i=0;i<list.size();i++)
		if(find(list[i],chain)) res.push_back(i);
}

void 	Team::replace(vector<string> &list)
{
	string chain;
	vector<int> reslist;
	int flag = 0;
	for(int j=0;j<list.size();j++)
	{
		BreakItem(list[j],chain,reslist);
		if(chain!="") flag|=PEP_FLAG;
		if(reslist.size()) flag|=NUM_FLAG;
		if(flag & NUM_FLAG)
		{
			for(int j=0;j<reslist.size();j++)
			{
				int s = table.size();
				table.resize(s+1);
				table[s].mask = flag;
				table[s].chain_id=chain;
				table[s].res_name="";
				table[s].atom_name="";
				table[s].res_id=reslist[j];
			}
		}
		else
		{
		int s = table.size();
		table.resize(s+1);
		table[s].mask = flag;
		table[s].chain_id=chain;
		table[s].res_name="";
		table[s].atom_name="";
		table[s].res_id=0;
		}
	}
}

Team::Team(vector<string> list)
{
	string chain;
	vector<int> reslist;
	int flag = 0;
	for(int j=0;j<list.size();j++)
	{
		BreakItem(list[j],chain,reslist);
		if(chain!="") flag|=PEP_FLAG;
		if(reslist.size()) flag|=NUM_FLAG;
		if(flag & NUM_FLAG)
		{
			for(int j=0;j<reslist.size();j++)
			{
				int s = table.size();
				table.resize(s+1);
				table[s].mask = flag;
				table[s].chain_id=chain;
				table[s].res_name="";
				table[s].atom_name="";
				table[s].res_id=reslist[j];
			}
		}
		else
		{
		int s = table.size();
		table.resize(s+1);
		table[s].mask = flag;
		table[s].chain_id=chain;
		table[s].res_name="";
		table[s].atom_name="";
		table[s].res_id=0;
		}
	}
}

Team::Team(char *file,vector<string> list)
{
	vector<string> arg1_list;
	vector<string> arg2_list;
	string chain;
	vector<int> reslist;
	FILE *fp=fopen(file,"r");
	if(!fp) return;
	while(!feof(fp))
	{
		char line[100];
		char word[100];
		if(fgets(line,80,fp)==NULL) break;
		line[strlen(line)-1]=0;
		int index=0;
		getword(line,word,index);
		int s;
		s=arg1_list.size();
		arg1_list.resize(s+1);
		arg1_list[s]=word;
		getword(line,word,index);
		s=arg2_list.size();
		arg2_list.resize(s+1);
		arg2_list[s]=word;
	}
	fclose(fp);

	for(int i=0;i<arg1_list.size();i++)
	{
		int flag = ARG1_FLAG | ARG2_FLAG;
		for(int j=0;j<list.size();j++)
		{
			BreakItem(list[j],chain,reslist);
			if(chain!="") flag|=PEP_FLAG;
			if(reslist.size()) flag|=NUM_FLAG;
			if(flag & NUM_FLAG)
			{
				for(int j=0;j<reslist.size();j++)
				{
					int s = table.size();
					table.resize(s+1);
					table[s].mask = flag;
					table[s].chain_id=chain;
					table[s].res_name=arg1_list[i];
					table[s].atom_name,arg2_list[i];
					table[s].res_id=reslist[j];
				}
			}
			else
			{
			int s = table.size();
			table.resize(s+1);
			table[s].mask = flag;
			table[s].chain_id=chain;
			table[s].res_name=arg1_list[i];
			table[s].atom_name,arg2_list[i];
			table[s].res_id=0;
			}
		}
		if(!list.size())
		{
			int s = table.size();
			table.resize(s+1);
			table[s].mask = flag;
			table[s].chain_id=chain;
			table[s].res_name=arg1_list[i];
			table[s].atom_name=arg2_list[i];
			table[s].res_id=0;
		}
	}
}

Team::Team(char *file,char *str1,char *str2)
{
	vector<int> num_list;
	vector<string> pep_list;
	vector<string> arg1_list;
	vector<string> arg2_list;
	string s1=str1,s2=str2;
	string t="";
	int index=0;
	while(index!=s1.size())
	{
		if(s1[index]!=',') t=t+s1[index];
		else
		{
			int s=pep_list.size();
			pep_list.resize(s+1);
			pep_list[s]=t;
			t="";
		}
		index++;
	}
	if(t.size())
	{
			int s=pep_list.size();
			pep_list.resize(s+1);
			pep_list[s]=t;
	}
	t="";
	index=0;
	int minus_flag=0;
	while(index!=s2.size())
	{
		if(s2[index]!=',' && s2[index]!='-') t=t+s2[index];
		else
		{
			if(minus_flag==0)
			{
				int s=num_list.size();
				num_list.resize(s+1);
				num_list[s]=atoi(t.c_str());
			}
			else
			{
				int start=num_list[num_list.size()-1];
				int end = atoi(t.c_str());
				for(int i=start+1;i<=end;i++)
				{
					int s=num_list.size();
					num_list.resize(s+1);
					num_list[s]=i;
				}
				minus_flag=0;
				
			}
			t="";
			if(s2[index]=='-') minus_flag=1;
		}
		index++;
	}
	if(t.size())
	{
			if(minus_flag==0)
			{
				int s=num_list.size();
				num_list.resize(s+1);
				num_list[s]=atoi(t.c_str());
			}
			else
			{
				int start=num_list[num_list.size()-1];
				int end = atoi(t.c_str());
				for(int i=start+1;i<=end;i++)
				{
					int s=num_list.size();
					num_list.resize(s+1);
					num_list[s]=i;
				}
				minus_flag=0;
				
			}
	}
	FILE *fp=fopen(file,"r");
	if(!fp) return;
	while(!feof(fp))
	{
		char line[100];
		char word[100];
		if(fgets(line,80,fp)==NULL) break;
		line[strlen(line)-1]=0;
		int index=0;
		getword(line,word,index);
		int s;
		s=arg1_list.size();
		arg1_list.resize(s+1);
		arg1_list[s]=word;
		getword(line,word,index);
		s=arg2_list.size();
		arg2_list.resize(s+1);
		arg2_list[s]=word;
	}
	fclose(fp);
	int flag=ARG1_FLAG | ARG2_FLAG;
	if(pep_list.size()) flag|=PEP_FLAG;
	if(num_list.size()) flag|=NUM_FLAG;
	for(int i=0;i<arg1_list.size();i++)
	{
		if(pep_list.size() && !num_list.size())
		{
			for(int j=0;j<pep_list.size();j++)
			{
				int s=table.size();
				table.resize(s+1);
				table[s].mask = flag;
				table[s].chain_id=pep_list[j];
				table[s].res_id=0;
				table[s].res_name=arg1_list[i];
				table[s].atom_name=arg2_list[i];
			}
		}
		if(!pep_list.size() && num_list.size())
		{
			for(int j=0;j<num_list.size();j++)
			{
				int s=table.size();
				table.resize(s+1);
				table[s].mask = flag;
				table[s].chain_id="";
				table[s].res_id=num_list[j];
				table[s].res_name=arg1_list[i];
				table[s].atom_name=arg2_list[i];
			}
		}
		else
		if(pep_list.size() && num_list.size())
		{
			for(int j=0;j<pep_list.size();j++)
			{
				for(int k=0;k<num_list.size();k++)
				{
				int s=table.size();
				table.resize(s+1);
				table[s].mask = flag;
				table[s].chain_id=pep_list[j];
				table[s].res_id=num_list[k];
				table[s].res_name=arg1_list[i];
				table[s].atom_name=arg2_list[i];
				}
			}
		}
		if(!pep_list.size() && !num_list.size())
		{
			int s=table.size();
			table.resize(s+1);
			table[s].mask = flag;
			table[s].chain_id="";
			table[s].res_id=0;
			table[s].res_name=arg1_list[i];
			table[s].atom_name=arg2_list[i];
		}
	}
}

static int isWildCardIn(string s)
{
	for(int i=0;i<s.size();i++) if(s[i]==WILDCARD) return i;
	return -1;
}

void	Team::setResList(vector<string> &s)
{
	if(res_list.size()!=0) return ;
	res_list.resize(s.size());
	for(int i=0;i<s.size();i++) res_list[i]=s[i];
}

/*	Return 1 if the atom is in the team and 0 otherwise.
 *	The function returns the chain of the atom in the second
 *	argument.
 * */
int Team::find(atom atom,string &chain)
{
	int ifound=0;
	
	for(int i=0;i<res_list.size();i++)
	{
		if(atom.res_name==res_list[i])
		{
			ifound=1;
			break;
		}
	}
	if(res_list.size() && !ifound) return 0;

	if(atom_type.size())
	{
		ifound=0;
		for(int i=0;i<atom_type.size();i++)
		{
			if(atom_type[i]=="noh" && atom.atom_name[0]!='H') {ifound=1;break;}
			if(atom_type[i]=="ca" && atom.atom_name=="CA") {ifound=1;break;}
			if(atom_type[i]=="sidechain" && isSidechain(atom)) {ifound=1;break;}
			if(atom_type[i]=="backbone" && isBackbone(atom)) {ifound=1;break;}
			if(atom_type[i]=="backbone3" && isBackbone3(atom)) {ifound=1;break;}
			if(atom_type[i]=="backbone4" && isBackbone4(atom)) {ifound=1;break;}
		}
		if(!ifound) return 0;
	}

	if(res_name.size())
	{
		ifound = 0;
		for(int i=0;i<res_name.size();i++)
		{
				if(atom.res_name==res_name[i]) {ifound=1;break;}
		}
		if(!ifound) return 0;
	}

	if(table.size()==0){
	       	return 1;
	}

	for(int i=0;i<table.size();i++)
	{
		int icount=0;
		if(table[i].mask & PEP_FLAG) 
		{
			if(table[i].chain_id==atom.chain_id || table[i].chain_id=="*")
				icount=1;
			else
				icount=-1;
			int wildPos=isWildCardIn(table[i].chain_id);
			if(wildPos!=-1)
			{
				if(!strncmp(table[i].chain_id.c_str(),atom.chain_id.c_str(),wildPos)) icount=1;
			}
		}
		if(icount==-1) continue;
		if(table[i].mask & NUM_FLAG)
		{
			if(table[i].res_id==atom.res_id) 
			{
				icount=1;
			}
			else
				icount=-1;
		}
		if(icount==-1) continue;
		if(table[i].mask & ARG1_FLAG) 
		{
			if(table[i].res_name==atom.res_name || table[i].res_name=="*")
				icount=1;
			else
				icount=-1;
		}
		if(icount==-1) continue;
		if(table[i].mask & ARG2_FLAG) 
		{
			if(table[i].atom_name==atom.atom_name || 
			  (table[i].atom_name=="C*" && atom.atom_name[0]=='C' && atom.atom_name!="CA")
			  || (table[i].atom_name=="*"))
				icount=1;
			else
				icount=-1;
		}
		if(icount==1) 
		{
			chain=table[i].chain_id;
			int res_found=1;
			int atom_found=1;
			if(res_name.size()!=0)
			{
				res_found=0;
				for(int l=0;l<res_name.size();l++)
					if(res_name[l]==atom.res_name) 
				{
					res_found=1;
					break;
				}
			}
			if(atom_name.size()!=0)
			{
				atom_found=0;
				for(int l=0;l<atom_name.size();l++)
				{
					if(atom_name[l]==atom.atom_name)
					{
						atom_found=1;
						break;
					}
				}
			}
			if(!res_found || !atom_found) return 0;
			return 1;
		}
	}
	return 0;
}

Team::Team(vector<string> arg1,vector<string> arg2)
{
	table.resize(arg1.size());
	for(int i=0;i<table.size();i++)
	{
		table[i].mask = ARG1_FLAG | ARG2_FLAG;
		table[i].res_name=arg1[i];
		table[i].atom_name=arg2[i];
	}
}

/*	Return all the chains in the team.
 * */
void	Team::enumerateChains(vector<string> &result)
{
	result.resize(0);
	for(int i=0;i<table.size();i++)
	{
		string s=table[i].chain_id;
		if(s.size()==0) continue;
		int found=0;
		for(int j=0;j<result.size();j++)
		{
			if(result[j]==s) {found=1;break;}
		}
		if(!found) result.push_back(s);
	}
}

/*	Return all the resid-s in the team.
 * */
void	Team::enumerateResId(string s,vector<int> &res)
{
	res.resize(0);
	for(int i=0;i<table.size();i++)
	{
		string str=table[i].chain_id;
		if(table[i].res_id==0) continue;
		if(s==str)
		{
			int found=0;
			for(int j=0;j<res.size();j++)
				if(res[j]==table[i].res_id) {found=1;break;}
			if(!found) 
			{
				res.push_back(table[i].res_id);
			}
		}
	}
	if(res.size() == 0)
	{
		extern vector<atom> table;
	for(int i=0;i<table.size();i++)
	{
		string str=table[i].chain_id;
		if(s==str)
		{
			int found=0;
			for(int j=0;j<res.size();j++)
				if(res[j]==table[i].res_id) {found=1;break;}
			if(!found) 
			{
				res.push_back(table[i].res_id);
			}
		}
	}
	}
}

/*	Set the sequence string.
 * */
void	Team::setSeqString(vector<string> list)
{
	vector<string> arg1_list;
	vector<string> arg2_list;
	arg1_list.resize(table.size());
	arg2_list.resize(table.size());
	string chain;
	vector<int> reslist;
	for(int i=0;i<table.size();i++)
	{
		arg1_list[i]=table[i].res_name;
		arg2_list[i]=table[i].atom_name;
	}
	table.resize(0);

	for(int i=0;i<arg1_list.size();i++)
	{
		int flag = ARG1_FLAG | ARG2_FLAG;
		for(int j=0;j<list.size();j++)
		{
			BreakItem(list[j],chain,reslist);
			if(chain!="") flag|=PEP_FLAG;
			if(reslist.size()) flag|=NUM_FLAG;
			if(flag & NUM_FLAG)
			{
				for(int j=0;j<reslist.size();j++)
				{
					int s = table.size();
					table.resize(s+1);
					table[s].mask = flag;
					table[s].chain_id=chain;
					table[s].res_name=arg1_list[i];
					table[s].atom_name=arg2_list[i];
					table[s].res_id=reslist[j];
				}
			}
			else
			{
			int s = table.size();
			table.resize(s+1);
			table[s].mask = flag;
			table[s].chain_id=chain;
			table[s].res_name=arg1_list[i];
			table[s].atom_name=arg2_list[i];
			table[s].res_id=0;
			}
		}
		if(!list.size())
		{
			int s = table.size();
			table.resize(s+1);
			table[s].mask = flag;
			table[s].chain_id=chain;
			table[s].res_name=arg1_list[i];
			table[s].atom_name=arg2_list[i];
			table[s].res_id=0;
		}
	}
}

/*	Print the components of the team to the stdout.
 * */
void Team::print()
{
	for(int i=0;i<table.size();i++)
	{
		printf("TEAM ELEMENT[%4d]= ",i);
		if(table[i].mask & PEP_FLAG) printf(" %s ",table[i].chain_id.c_str());
		if(table[i].mask & NUM_FLAG) printf(" %6d ",table[i].res_id);
		if(table[i].mask & ARG1_FLAG) printf(" %s ",table[i].res_name.c_str());
		if(table[i].mask & ARG2_FLAG) printf(" %s ",table[i].atom_name.c_str());
		printf("\n");
	}
}

void	Team::setResName(vector<string> &s)
{
	res_name.resize(s.size());
	for(int i=0;i<s.size();i++) res_name[i]=s[i];
}

void	Team::setAtomName(vector<string> &s)
{
	atom_name.resize(s.size());
	for(int i=0;i<s.size();i++) atom_name[i]=s[i];
}

void	Team::setAtomType(vector<string> &s)
{
	atom_type.resize(s.size());
	for(int i=0;i<s.size();i++) atom_type[i]=s[i];
}

int	Team::getAtomNameSize()
{
	return atom_name.size();
}

Team::~Team()
{
}
