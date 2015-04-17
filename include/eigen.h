# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <util.h>
# include <vector>
# include <string>
using namespace std;

typedef vector<Data> DoubleVector;
typedef vector<DoubleVector> DoubleVector2;

int	eigen(int n,int m,DoubleVector2 &data,DoubleVector2 &symmat,
	DoubleVector &evals);
