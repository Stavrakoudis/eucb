
# ifndef __RAD__H
# define __RAD__H
# include <string>
# include <vector>
using namespace std;
typedef struct
{
	string res_name;
	string atom_name;
	double rad;
}RadStruct;

extern vector<RadStruct> rad;

extern void makeRad();
# endif
