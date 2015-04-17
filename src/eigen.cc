# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# include <util.h>
# include <vector>
# include <string>
using namespace std;

typedef vector<Data> DoubleVector;
typedef vector<DoubleVector> DoubleVector2;


#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )


void erhand(string err_msg)
{
    fprintf(stderr,"Run-time error:\n");
    fprintf(stderr,"%s\n", err_msg.c_str());
    fprintf(stderr,"Exiting to system.\n");
    exit(1);
}

void covcol(DoubleVector2 &data, int n, int m, DoubleVector2 &symmat)
{
 DoubleVector mean;
	mean.resize(m);
	
int i, j, j1, j2;

/* Allocate storage for mean vector */


/* Determine mean of column vectors of input data matrix */

for (j = 1; j <= m; j++)
    {
    mean[j-1] = 0.0;
    for (i = 1; i <= n; i++)
        {
        mean[j-1] += data[i-1][j-1];
        }
    mean[j-1] /= (float)n;
    }


/* Center the column vectors. */

for (i = 1; i <= n; i++)
    {
    for (j = 1; j <= m; j++)
        {
        data[i-1][j-1] -= mean[j-1];
        }
    }

/* Calculate the m * m covariance matrix. */
for (j1 = 1; j1 <= m; j1++)
    {
    for (j2 = j1; j2 <= m; j2++)
        {
        symmat[j1-1][j2-1] = 0.0;
        for (i = 1; i <= n; i++)
            {
            symmat[j1-1][j2-1] += data[i-1][j1-1] * data[i-1][j2-1];
            }
        symmat[j2-1][j1-1] = symmat[j1-1][j2-1];
        }
    }

return;

}




void tred2(DoubleVector2 &a, int n, DoubleVector &d, DoubleVector &e)
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
{
int l, k, j, i;
float scale, hh, h, g, f;

for (i = n; i >= 2; i--)
    {
    l = i - 1;
    h = scale = 0.0;
    if (l > 1)
       {
       for (k = 1; k <= l; k++)
           scale += fabs(a[i-1][k-1]);
       if (scale == 0.0)
          e[i-1] = a[i-1][l-1];
       else
          {
          for (k = 1; k <= l; k++)
              {
              a[i-1][k-1] /= scale;
              h += a[i-1][k-1] * a[i-1][k-1];

     }
          f = a[i-1][l-1];
          g = f>0 ? -sqrt(h) : sqrt(h);
          e[i-1] = scale * g;
          h -= f * g;
          a[i-1][l-1] = f - g;
          f = 0.0;
          for (j = 1; j <= l; j++)
              {
              a[j-1][i-1] = a[i-1][j-1]/h;
              g = 0.0;
              for (k = 1; k <= j; k++)
                  g += a[j-1][k-1] * a[i-1][k-1];
              for (k = j+1; k <= l; k++)
                  g += a[k-1][j-1] * a[i-1][k-1];
              e[j-1] = g / h;
              f += e[j-1] * a[i-1][j-1];
              }
          hh = f / (h + h);
          for (j = 1; j <= l; j++)
              {
              f = a[i-1][j-1];
              e[j-1] = g = e[j-1] - hh * f;
              for (k = 1; k <= j; k++)
                  a[j-1][k-1] -= (f * e[k-1] + g * a[i-1][k-1]);
              }
         }
    }
    else
        e[i-1] = a[i-1][l-1];
    d[i-1] = h;
    }
d[1-1] = 0.0;
e[1-1] = 0.0;
for (i = 1; i <= n; i++)
    {
    l = i - 1;
    if (d[i-1])
       {
       for (j = 1; j <= l; j++)
           {
           g = 0.0;
           for (k = 1; k <= l; k++)
               g += a[i-1][k-1] * a[k-1][j-1];
           for (k = 1; k <= l; k++)
               a[k-1][j-1] -= g * a[k-1][i-1];
           }
       }
       d[i-1] = a[i-1][i-1];
       a[i-1][i-1] = 1.0;
       for (j = 1; j <= l; j++)
           a[j-1][i-1] = a[i-1][j-1] = 0.0;
    }
}


int tqli(DoubleVector &d, DoubleVector &e, int n, DoubleVector2 &z)
{
int m, l, iter, i, k;
float s, r, p, g, f, dd, c, b;

for (i = 2; i <= n; i++)
    e[i-1-1] = e[i-1];
e[n-1] = 0.0;
for (l = 1; l <= n; l++)
    {
    iter = 0;
    do
      {
      for (m = l; m <= n-1; m++)
          {
          dd = fabs(d[m-1]) + fabs(d[m+1-1]);
          if (fabs(e[m-1]) + dd == dd) break;
          }
          if (m != l)
             {
             if (iter++ == 200) return 0;//erhand("No convergence in TLQI.");
             g = (d[l+1-1] - d[l-1]) / (2.0 * e[l-1]);
             r = sqrt((g * g) + 1.0);
             g = d[m-1] - d[l-1] + e[l-1] / (g + SIGN(r, g));
             s = c = 1.0;
             p = 0.0;
             for (i = m-1; i >= l; i--)
                 {
                 f = s * e[i-1];
                 b = c * e[i-1];

    if (fabs(f) >= fabs(g))
                    {
                    c = g / f;
                    r = sqrt((c * c) + 1.0);
                    e[i+1-1] = f * r;
                    c *= (s = 1.0/r);
                    }
                 else
                    {
                    s = f / g;
                    r = sqrt((s * s) + 1.0);
                    e[i+1-1] = g * r;
                    s *= (c = 1.0/r);
                    }
                 g = d[i+1-1] - p;
                 r = (d[i-1] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i+1-1] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++)
                     {
                     f = z[k-1][i+1-1];
                     z[k-1][i+1-1] = s * z[k-1][i-1] + c * f;
                     z[k-1][i-1] = c * z[k-1][i-1] - s * f;
                     }
                 }
                 d[l-1] = d[l-1] - p;
                 e[l-1] = g;
                 e[m-1] = 0.0;
             }
          }  while (m != l);
      }
	return 1;
 }



void corcol(DoubleVector2 &data, int n, int m, DoubleVector2 &symmat)
{
float eps = 0.005;
float x;
int i, j, j1, j2;

/* Allocate storage for mean and std. dev. vectors */

DoubleVector mean;mean.resize(m);
DoubleVector stddev; stddev.resize(m);

/* Determine mean of column vectors of input data matrix */

for (j = 1; j <= m; j++)
    {
    mean[j-1] = 0.0;
    for (i = 1; i <= n; i++)
        {
        mean[j-1] += data[i-1][j-1];
        }
    mean[j-1] /= (float)n;
    }

for (j = 1; j <= m; j++)
    {

	  stddev[j-1] = 0.0;
    for (i = 1; i <= n; i++)
        {
        stddev[j-1] += (   ( data[i-1][j-1] - mean[j-1] ) *
                         ( data[i-1][j-1] - mean[j-1] )  );
        }
        stddev[j-1] /= (float)n;
        stddev[j-1] = sqrt(stddev[j-1]);
        /* The following in an inelegant but usual way to handle
        near-zero std. dev. values, which below would cause a zero-
        divide. */
        if (stddev[j-1] <= eps) stddev[j-1] = 1.0;
    }

/* Center and reduce the column vectors. */

for (i = 1; i <= n; i++)
    {
    for (j = 1; j <= m; j++)
        {
        data[i-1][j-1] -= mean[j-1];
        x = sqrt((float)n);
        x *= stddev[j-1];
        data[i-1][j-1] /= x;
        }
    }
/* Calculate the m * m correlation matrix. */
for (j1 = 1; j1 <= m-1; j1++)
    {
    symmat[j1-1][j1-1] = 1.0;
    for (j2 = j1+1; j2 <= m; j2++)
        {
        symmat[j1-1][j2-1] = 0.0;
        for (i = 1; i <= n; i++)
            {
            symmat[j1-1][j2-1] += ( data[i-1][j1-1] * data[i-1][j2-1]);
            }
        symmat[j2-1][j1-1] = symmat[j1-1][j2-1];
        }
    }
    symmat[m-1][m-1] = 1.0;

return;

}


int	eigen(int n,int m,DoubleVector2 &data,DoubleVector2 &symmat,
	DoubleVector &evals)
{
	int i,j;
	
	DoubleVector interm;
	DoubleVector2 symmat2;

	symmat.resize(m);
	symmat2.resize(m);
	for(i=0;i<m;i++) 
	{
		symmat[i].resize(m);
		symmat2[i].resize(m);
	}
	evals.resize(m);
	interm.resize(m);

          //covcol(data, n, m, symmat);
    for (i = 1; i <= m; i++) {
     for (j = 1; j <= m; j++) {
		symmat[i-1][j-1]=data[i-1][j-1];
      symmat2[i-1][j-1] = symmat[i-1][j-1]; /* Needed below for col. projections */
                              }
                             }
    tred2(symmat, m, evals, interm);  /* Triangular decomposition */
    if(!tqli(evals, interm, m, symmat)) 
	return 0;   /* Reduction of sym. trid. matrix */
	
}

