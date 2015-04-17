# include <pdb.h>
# include <psf.h>
# include <math.h>

PointVector pdbpos;
vector<PdbAtom> pdb_table;
string pdb_name;

Data	getPdbDistance(PdbAtom a,PdbAtom b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}

int readPdb(string filename)
{
	FILE *fp=fopen(filename.c_str(),"r");
	if(!fp) 
	{
		return 0;
	}
	pdb_name = filename;
	vector<string> list;
	extern vector<atom> table;
	int icount = 0;
	while(1)
	{
		char str[1024];
		if(!fgets(str,1023,fp)) break;
		str[strlen(str)-1]=0;
		list.resize(0);
		getAllWords(str,list);
		if(list[0]=="ATOM")
		{
			int atom_id=list[1][0]=='*'?table[icount].atom_id:atoi(list[1].c_str());
			string atom_name=list[2];
			string res_name=list[3];
			int special = (res_name=="TIP3" || res_name=="SOD"|| res_name=="CLA");
		        if(res_name[0]=='T' && res_name[1]=='I' && res_name[2]=='P' ||
				res_name[3]=='3')
			{
				res_name="TIP3";
				special=1;
			}
			string chain_id=special?list[10]:list[11];
			int    res_id = special?atoi(list[4].c_str()):atoi(list[5].c_str());
			string pdb_chain = list[4];
			if(atom_id!=table[icount].atom_id) 
			{
				printf("different atom id \n");
				fclose(fp);
				return 0;
			}
			if(atom_name!=table[icount].atom_name)
			{
				printf("different atom name \n");
				fclose(fp);
				return 0;
			}
			if(res_name!=table[icount].res_name) 
			{
				
				printf("different res name %s versus %s in line %d of pdb file \n",
					res_name.c_str(),table[icount].res_name.c_str(),icount);
				fclose(fp);
				return 0;
			}
			if(chain_id!=table[icount].chain_id) 
			{
				printf("different chain id  atom = %d first = %s second = %s\n",icount,table[icount].chain_id.c_str(),chain_id.c_str());
				fclose(fp);
				return 0;
			}
			if(res_id!=table[icount].res_id) 
			{
				printf("different red id \n");
				fclose(fp);
				return 0;
			}
			Point pt;
			pt.x=special?atof(list[5].c_str()):atof(list[6].c_str());
			pt.y=special?atof(list[6].c_str()):atof(list[7].c_str());
			pt.z=special?atof(list[7].c_str()):atof(list[8].c_str());
			pdbpos.push_back(pt);
			pdb_table.resize(icount+1);

			pdb_table[icount].atom_id=atom_id;
			pdb_table[icount].atom_name = atom_name;
			pdb_table[icount].res_name = res_name;
			pdb_table[icount].pdb_chain = special?"":pdb_chain;
			pdb_table[icount].res_id = res_id;
			pdb_table[icount].x = pt.x;
			pdb_table[icount].y = pt.y;
			pdb_table[icount].z = pt.z;
			pdb_table[icount].val1 = special?atof(list[8].c_str()):atof(list[9].c_str());
			pdb_table[icount].val2 = special?atof(list[9].c_str()):atof(list[10].c_str());
			pdb_table[icount].chain_id=chain_id;
			pdb_table[icount].pdb_letter = special?list[11]:list[12];

			icount++;
		}
	}
	fclose(fp);
	return 1;
}
