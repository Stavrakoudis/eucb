# include <global.h>
# include <command_anal.h>

typedef struct
{
	string 	chain;
	int    	firstRes;
	vector<string> resName;
	vector<int>    countResName;
	int	lastRes;
	int	resCount;
}analStruct;

typedef struct
{
	string atom1_chainid;
	string atom1_resname;
	int atom1_resid;
	string atom2_chainid;
	string atom2_resname;
 	int atom2_resid;
}conStruct;

void	sortAnalStruct(analStruct &anal)
{
	for(int i=0;i<anal.resName.size();i++)
	{
		for(int j=0;j<anal.resName.size()-1;j++)
		{
			if(anal.resName[j+1]<anal.resName[j])
			{
				string tstr=anal.resName[j+1];
				anal.resName[j+1]=anal.resName[j];
				anal.resName[j]=tstr;
				int itemp;
				itemp=anal.countResName[j+1];
				anal.countResName[j+1]=anal.countResName[j];
				anal.countResName[j]=itemp;
			}
		}
	}
}

CommandAnal::CommandAnal()
	:Command("ANAL")
{
}

/*	This function performs the psf analyzing.
 * */
void CommandAnal::Run()
{
	Team *psfanal_team = team;
	vector<analStruct> AnalTable;
	vector<conStruct>  ConTable;
	int last_resid=0;
	for(int i=0;i<table.size();i++)
	{
		string st;
		if(psfanal_team!=NULL && !psfanal_team->find(table[i],st)) continue;

		if(table[i].res_name=="CYS" && table[i].atom_name=="SG")
		{
			for(int j=0;j<table[i].bond.size();j++)
			{
				int id=table[i].bond[j];
				if(table[id-1].res_name=="CYS" && table[id-1].atom_name=="SG")
				{

					int ifound=0;
					for(int k=0;k<ConTable.size();k++)
					{
						if(ConTable[k].atom2_chainid==table[i].chain_id
							&& ConTable[k].atom2_resname==table[i].res_name
							&& ConTable[k].atom2_resid==table[i].res_id
							&& ConTable[k].atom1_chainid==table[id-1].chain_id
							&& ConTable[k].atom1_resname==table[id-1].res_name
							&& ConTable[k].atom1_resid==table[id-1].res_id) 
						{
							ifound=1;
							break;
						}
					}

					if(ifound) continue;

					int s=ConTable.size();
					ConTable.resize(s+1);
					ConTable[s].atom1_chainid=table[i].chain_id;
					ConTable[s].atom1_resname=table[i].res_name;
					ConTable[s].atom1_resid=table[i].res_id;
					ConTable[s].atom2_chainid=table[id-1].chain_id;
					ConTable[s].atom2_resname=table[id-1].res_name;
					ConTable[s].atom2_resid=table[id-1].res_id;
				}
			}
		}

		string str=table[i].chain_id;
		int found=-1;
		for(int j=0;j<AnalTable.size();j++)
		{
			if(AnalTable[j].chain==str) {found=j;break;}
		}
		if(found==-1)
		{
			int s=AnalTable.size();
			AnalTable.resize(s+1);
			AnalTable[s].chain=str;
			AnalTable[s].firstRes=table[i].res_id;
			AnalTable[s].lastRes=table[i].res_id;
			AnalTable[s].resCount=1;
			AnalTable[s].resName.push_back(table[i].res_name);
			AnalTable[s].countResName.push_back(1);
			last_resid=table[i].res_id;
		}
		else
		{
			int id=table[i].res_id;
			if(id!=last_resid) 
			{
				AnalTable[found].resCount++;
				int ifound=-1;
				for(int j=0;j<AnalTable[found].resName.size();j++)
				{
					if(AnalTable[found].resName[j]==table[i].res_name)
					{
						ifound=j;
						break;
					}
				}
				if(ifound==-1)
				{
					AnalTable[found].resName.push_back(table[i].res_name);
					AnalTable[found].countResName.push_back(1);
				}
				else
				{
					AnalTable[found].countResName[ifound]++;
				}
				last_resid=id;
			}
			if(id>AnalTable[found].lastRes)
				AnalTable[found].lastRes=id;
			if(id<AnalTable[found].firstRes)
				AnalTable[found].firstRes=id;
		}
	}
	Log("Found "+printNumber((int)AnalTable.size())+" chains ");
	printf("ENUMERATION OF CHAINS IN PSF \n");
	printf("================================================================================= \n");
	for(int i=0;i<AnalTable.size();i++)
	{
		printf("%d) CHAIN: %s RESIDUES:%d  RANGE: %d-%d\n",i+1,
				AnalTable[i].chain.c_str(),
				AnalTable[i].resCount,
				AnalTable[i].firstRes,AnalTable[i].lastRes);
		printf("\t   RESNAME    COUNT    \n");
		printf("================================\n");
		sortAnalStruct(AnalTable[i]);
		for(int j=0;j<AnalTable[i].resName.size();j++)
		{
			printf(" %4d) %10s %8d\n",j+1,AnalTable[i].resName[j].c_str(),AnalTable[i].countResName[j]);
		}
	}
	printf("================================================================================= \n");

	printf("Disulfide Bonds \n");
	printf("================================================================================= \n");
	for(int i=0;i<ConTable.size();i++)
	{
		printf("%s:%s:%d->%s:%s:%d\n",
					ConTable[i].atom1_chainid.c_str(),
					ConTable[i].atom1_resname.c_str(),
					ConTable[i].atom1_resid,
					ConTable[i].atom2_chainid.c_str(),
					ConTable[i].atom2_resname.c_str(),
					ConTable[i].atom2_resid);
	}
	printf("================================================================================= \n");
}
