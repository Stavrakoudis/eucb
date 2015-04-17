# include <string.h>
# include <stdlib.h>
# include <stdio.h>
# include <vector>
# include <string>
using namespace std;

int main()
{
	FILE *fout=fopen("helpfiles.cc","w");
	fprintf(fout,"#include <vector>\n");
	fprintf(fout,"#include <string>\n");
	fprintf(fout,"using namespace std;\n");
	fprintf(fout,"vector<string> donorsrc_table;\n");
	fprintf(fout,"vector<string> eucbhelp_table;\n");
	fprintf(fout,"void makeHelpTables()\n");
	fprintf(fout,"{\n");
	FILE *fp=fopen("donorsrc","r");
	char line[1024];
	while(1)
	{
	 if(!fgets(line,1023,fp)) break;
	 line[strlen(line)-1]=0;
	 fprintf(fout,"donorsrc_table.push_back(\"%s\");\n",line);
	}
	fclose(fp);

	fp=fopen("eucbhelp","r");
	while(1)
	{
	 if(!fgets(line,1023,fp)) break;
	 line[strlen(line)-1]=0;
	 fprintf(fout,"eucbhelp_table.push_back(\"%s\");\n",line);
	}
	fclose(fp);
	fprintf(fout,"}\n");
	fclose(fout);
	return 0;
}
