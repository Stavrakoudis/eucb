# include <stdlib.h>
# include <string.h>
# include <string.h>
# include <util.h>

typedef struct
{
	string res_name;
	string atom_name;
	float  rad;
}RadStruct;

vector<RadStruct> rad;

int main(int argc,char **argv)
{
	FILE *fp=fopen(argv[1],"r");
	if(!fp) return 1;
	printf("# include <string>\n");
	printf("# include <vector>\n");
	printf("using namespace std;\n");
	printf("typedef struct\n{\n");
	printf("\tstring res_name;\n");
	printf("\tstring atom_name;\n");
	printf("\tdouble rad;\n");
	printf("}RadStruct;\n");
	while(1)
	{
		char s[1024];
		vector<string> str;
		if(!fgets(s,1023,fp)) break;
		if(strlen(s)==0) break;
		getAllWords(s,str);	
		if(str[0]=="RESIDUE" && str[1]=="ATOM")
		{
			string res_name=str[2];
			int count=atoi(str[3].c_str());	
			for(int i=0;i<count;i++)
			{
				char s1[1024];
				vector<string> str1;
				if(!fgets(s1,1023,fp)) break;
				getAllWords(s1,str1);	
				RadStruct st;
				st.res_name=res_name;
				st.atom_name=str1[1];
				st.rad=atof(str1[2].c_str());
				rad.push_back(st);
				fprintf(stderr,"%s\t%s\t%lf\t%s\n",
					res_name.c_str(),str1[1].c_str(),st.rad,
					str1[3].c_str());
			}
		}
	}
	printf("vector<RadStruct> rad;\n");
	printf("void makeRad()\n{\n");
	printf("rad.resize(%d);\n",rad.size());
	for(int i=0;i<rad.size();i++) 
	{
		printf("rad[%d].res_name=\"%s\";rad[%d].atom_name=\"%s\";rad[%d].rad=%f;\n",i,rad[i].res_name.c_str(),i,rad[i].atom_name.c_str(),i,
		rad[i].rad);
	}
	printf("}\n");
	fclose(fp);
	return 0;
}
