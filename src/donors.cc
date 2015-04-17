# include <donors.h>
# include <psf.h>
# include <team.h>
# include <helpfiles.h>
/*	Donor: all the possible donors.
 *	Hydro: all the possible hydrogens.
 *	Acceptor: all the possible acceptors.
 * */
vector<int> Donor;
vector<int> Hydro;
vector<int> Acceptor;
vector<IntVector> DonorCon;
vector<IntVector> AcceptorCon;


void	parseDonorFile(vector<atom> &table)
{
	for(int i=0;i<table.size();i++)
	{
		//water acceptor
		if(table[i].res_name=="TIP3" && table[i].atom_name=="OH2")
		{
			 Acceptor.push_back(table[i].atom_id);
		}
		else //protein acceptor
		if(table[i].atom_name[0]=='O') // thanos
		{
			if(table[i].dnum1<0)
			{
				int ifound=0;
				for(int j=0;j<table[i].bond.size();j++)
				{
					int id=table[i].bond[j];
					if(table[id-1].atom_name[0]=='C')
					{
						ifound=1;
						break;
					}
				}
				if(ifound) Acceptor.push_back(table[i].atom_id);
			}
		}
		//water donor
		if(table[i].res_name=="TIP3" && table[i].atom_name[0]=='O')
		{
			int foundH=0;
			for(int j=0;j<table[i].bond.size();j++)
			{
				 int id=table[i].bond[j];
				 if(table[id-1].atom_name[0]=='H') foundH++;
				 if(foundH==2) break;
			}
			if(foundH==2) 
			{
				Donor.push_back(table[i].atom_id);
			}
		}
		else //protein donor
		if(table[i].atom_name[0]=='O' || table[i].atom_name[0]=='N'
		|| table[i].atom_name[0]=='S')
		{
		//	if(table[i].dnum1>0)
			{
				int foundC=0;
				int foundH=0;
				for(int j=0;j<table[i].bond.size();j++)
				{
				 int id=table[i].bond[j];
				 if(table[id-1].atom_name[0]=='C') foundC=1;
				 if(table[id-1].atom_name[0]=='H') foundH=1;
				 if(foundC && foundH) break;
				}
				if(foundC && foundH) 
				{
					Donor.push_back(table[i].atom_id);
				}
			}
		}
	}
	DonorCon.resize(Donor.size());
	AcceptorCon.resize(Acceptor.size());
	for(int i=0;i<Donor.size();i++)
	{
		int donor_id=Donor[i];
		for(int j=0;j<table[donor_id-1].bond.size();j++)
		{
			int id=table[donor_id-1].bond[j];
			if(table[id-1].atom_name[0]=='H') DonorCon[i].push_back(id);
			if(table[id-1].atom_name[0]=='H')	
			{
				int ifound=0;
				for(int k=0;k<Hydro.size();k++)
				{
					if(Hydro[k]==id) {ifound=1;break;}
				}
				if(!ifound) Hydro.push_back(id);
			}
		}	
	}
	for(int i=0;i<Acceptor.size();i++)
	{
		int acc_id=Acceptor[i];
		for(int j=0;j<table[acc_id-1].bond.size();j++)
		{
			int id=table[acc_id-1].bond[j];
			AcceptorCon[i].push_back(id);
		}
	}
	/*
	vector<DonorT> donor_table2;

	char line[1025];
	char arg1_word[1024];
	char arg2_word[1024];
	char arg3_word[1024];
	int count = 0;
	char word[1024];
	for(int icount=0;icount<donorsrc_table.size();icount++)
	{
		strcpy(line,donorsrc_table[icount].c_str());
		char word[1024];
		int index=0;
		getword(line,word,index);
		if(!strcmp(word,"RESI"))
		{

			getword(line,arg1_word,index);
		}
		else
		if(!strcmp(word,"ACCEPTOR"))
		{
			getword(line,arg2_word,index);
			getword(line,arg3_word,index);
			count++;
			donor_table2.resize(count);
			donor_table2[count-1].flag=ISACCEPTOR;
			donor_table2[count-1].res_name=arg1_word;
			donor_table2[count-1].atom_name=arg2_word;
			donor_table2[count-1].atom_type=arg3_word;
		}
		else
		if(!strcmp(word,"DONOR"))
		{
			getword(line,arg2_word,index);
			getword(line,arg3_word,index);

			count++;
			donor_table2.resize(count);
			donor_table2[count-1].flag=ISHYDRO;
			donor_table2[count-1].res_name=arg1_word;
			donor_table2[count-1].atom_name=arg2_word;
			donor_table2[count-1].atom_type=arg3_word;

			count++;
			donor_table2.resize(count);
			donor_table2[count-1].flag=ISDONOR;
			donor_table2[count-1].res_name=arg1_word;
			donor_table2[count-1].atom_name=arg2_word;
			donor_table2[count-1].atom_type=arg3_word;
		}
	}
	Donor.resize(0);
	Hydro.resize(0);
	Acceptor.resize(0);
	for(int i=0;i<donor_table2.size();i++)
	{
		for(int j=0;j<table.size();j++)
		{
			if(table[j].res_name!=donor_table2[i].res_name) continue;
			if(donor_table2[i].flag == ISHYDRO)
			{
				if(table[j].atom_name==donor_table2[i].atom_name)
				{
					int s=Hydro.size();

					int found=0;
					for(int k=0;k<s;k++)
						if(Hydro[k]==table[j].atom_id) 
						{
							found=1;
							break;
						}
					if(found) continue;
					table[j].ishydro=1;
					Hydro.push_back(table[j].atom_id);
				}
			}
			
			if(donor_table2[i].flag == ISDONOR)
			{
				if(table[j].atom_name==donor_table2[i].atom_type)
				{
					int s=Donor.size();
					int found=0;
					for(int k=0;k<s;k++)
						if(Donor[k]==table[j].atom_id) 
						{
							found=1;
							break;
						}
					if(found) continue;
					Donor.resize(s+1);
					Donor[s]=table[j].atom_id;
					DonorCon.resize(s+1);
					DonorCon[s].resize(0);
				}
			}
			if(donor_table2[i].flag == ISACCEPTOR)
			{
				if(table[j].atom_name==donor_table2[i].atom_name)
				{
					int s=Acceptor.size();
					int found=0;
					for(int k=0;k<s;k++)
						if(Acceptor[k]==table[j].atom_id) 
						{
							found=1;
							break;
						}
					if(found) continue;
					Acceptor.resize(s+1);
					Acceptor[s]=table[j].atom_id;
					AcceptorCon.resize(s+1);
					
					AcceptorCon[s].resize(table[j].bond.size());
					for(int l=0;l<table[j].bond.size();l++)
					{
						int id =table[j].bond[l];
						AcceptorCon[s][l]=table[id-1].atom_id;
					}
				}
			}
		}
	}
	for(int i=0;i<Donor.size();i++)
	{
		int s=table[Donor[i]-1].bond.size();
		for(int j=0;j<s;j++)
		{
			int id = table[Donor[i]-1].bond[j];
			if(table[id-1].ishydro) 
			{
				DonorCon[i].push_back(id);
			}
		}
	}
	*/
}


/*	Return 1 if the atom with atom_id id is a Water Donor.
 * */
int	isWaterDonor(int id,vector<atom> &table)
{	
	return table[id-1].res_name=="TIP3";
}

/*	Return 1 if the atom with atom_id id is a Protein Donor.
 * */
int	isProteinDonor(int id,vector<atom> &table)
{
	return !isWaterDonor(id,table);
}

/*	Return 1 if the atom with atom_id id is a Backbone Donor.
 * */
int	isBackBoneDonor(int id,vector<int> &con,vector<atom> &table)
{
	int ifound=0;
	if(table[id-1].atom_name!="N") return 0;
	for(int i=0;i<con.size();i++)
	{
		if(table[con[i]-1].atom_name=="HN")
		{
			ifound=1;
			break;
		}
	}
	return ifound;
}

/*	Return 1 if the atom with atom_id id is a Sidechain Donor.
 * */
int	isSideChainDonor(int id,vector<int> &con,vector<atom> &table)
{
	if(!isProteinDonor(id,table)) return 0;
	if(isBackBoneDonor(id,con,table)) return 0;
	return 1;
}

/*	Return 1 if the atom with atom_id id is a Water acceptor.
 * */
int	isWaterAcceptor(int id,vector<atom> &table)
{
	return table[id-1].res_name=="TIP3";
}

/*	Return 1 if the atom with atom_id id is a Protein acceptor.
 * */
int	isProteinAcceptor(int id,vector<atom> &table)
{
	return !isWaterAcceptor(id,table);
}

/*	Return 1 if the atom with atom_id id is a BackBone acceptor.
 * */
int	isBackBoneAcceptor(int id,vector<int> &con,vector<atom> &table)
{
	if(!isProteinAcceptor(id,table)) return 0;
	if(table[id-1].atom_name!="O") return 0;
	int ifound=0;
	for(int i=0;i<con.size();i++)
	{
		if(table[con[i]-1].atom_name!="C")
		{
			ifound=1;
			break;
		}
	}	
	return ifound;
}

/*	Return 1 if the atom with atom_id id is a sidechain acceptor.
 * */
int	isSideChainAcceptor(int id,vector<int> &con,vector<atom> &table)
{
	if(!isProteinAcceptor(id,table)) return 0;
	if(isBackBoneAcceptor(id,con,table)) return 0;
	return 1;
}

/*	Parse the donorsrc file located in the HOME directory of the user.
 * */
void	parse_donorsrc(vector<atom> &table)
{
	parseDonorFile(table);
}

void SelectDonorAcceptor(vector<int> &CopyDonor,vector<IntVector> &CopyDonorCon,
	vector<int> &CopyAcceptor,vector<IntVector> &CopyAcceptorCon,
	Team *donorTeam,Team *acceptorTeam)
{
	extern vector<atom> table;
	extern string donor_selection;
	extern string acceptor_selection;
	for(int i=0;i<Donor.size();i++)
	{
		int donor_id=Donor[i];
		string s1;
		if(donorTeam==NULL || donorTeam->find(table[donor_id-1],s1))
		{
			if(donor_selection=="backbone")
			{
				if(!isBackBoneDonor(donor_id,DonorCon[i],table))  continue;
			}
			else
			if(donor_selection=="sidechain")
			{
				if(!isSideChainDonor(donor_id,DonorCon[i],table)) continue;
			}
			else
			if(donor_selection=="water")
			{
				if(!isWaterDonor(donor_id,table)) continue;
			}
			else
			if(donor_selection=="protein")
			{
				if(!isProteinDonor(donor_id,table)) continue;
			}
			int s=CopyDonor.size();
			CopyDonor.resize(s+1);
			CopyDonor[s]=donor_id;
			CopyDonorCon.resize(s+1);
			CopyDonorCon[s].resize(DonorCon[i].size());
			for(int j=0;j<DonorCon[i].size();j++) CopyDonorCon[s][j]=DonorCon[i][j];
			
		}
	}
	for(int i=0;i<Acceptor.size();i++)
	{
		int acceptor_id=Acceptor[i];
		string s1;
		if(acceptorTeam==NULL || acceptorTeam->find(table[acceptor_id-1],s1))
		{
			if(acceptor_selection=="backbone")
			{
				if(!isBackBoneAcceptor(acceptor_id,AcceptorCon[i],table)) continue;
			}
			else
			if(acceptor_selection=="sidechain")
			{
				if(!isSideChainAcceptor(acceptor_id,AcceptorCon[i],table)) continue;
			}
			else
			if(acceptor_selection=="water")
			{
				if(table[acceptor_id-1].res_name!="TIP3")	continue;
			//	if(!isWaterAcceptor(acceptor_id,table)) continue;
			}
			else
			if(acceptor_selection=="protein")
			{
				if(!isProteinAcceptor(acceptor_id,table)) continue;
			}
			int s=CopyAcceptor.size();
			CopyAcceptor.resize(s+1);
			CopyAcceptor[s]=acceptor_id;
			CopyAcceptorCon.resize(s+1);
			CopyAcceptorCon[s].resize(AcceptorCon[i].size());
			for(int j=0;j<AcceptorCon[i].size();j++)
			{
				CopyAcceptorCon[s][j]=AcceptorCon[i][j];
			}
			
		}
	}
}
