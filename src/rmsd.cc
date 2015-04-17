#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
# include <util.h>
typedef struct
{
  Data m[4][4];
} MATRIX;

#define vdiff2(a,b) ( ((a)[0]-(b)[0]) * ((a)[0]-(b)[0]) +	\
  ((a)[1]-(b)[1]) * ((a)[1]-(b)[1]) + \
  ((a)[2]-(b)[2]) * ((a)[2]-(b)[2]) )

static Data alignedrmsd(Data *v1, Data *v2, int N);
static void centroid(Data *ret, Data *v, int N);
static int getalignmtx(Data *v1, Data *v2, int N, MATRIX *mtx);
static void crossproduct(Data *ans, Data *pt1, Data *pt2);
static void mtx_root(MATRIX *mtx);
static int almostequal(MATRIX *a, MATRIX *b);
void mulpt(MATRIX *mtx, Data *pt);
static void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y);
static void mtx_identity(MATRIX *mtx);
static void mtx_trans(MATRIX *mtx, Data x, Data y, Data z);
static int mtx_invert(Data *mtx, int N);
static Data absmaxv(Data *v, int N);

/*
  calculate rmsd between two structures
  Params: v1 - first set of points
          v2 - second set of points
          N - number of points
          mtx - return for transfrom matrix used to align structures
  Returns: rmsd score
  Notes: mtx can be null. Transform will be rigid. Inputs must
         be previously aligned for sequence alignment
 */
Data rmsd(Data *v1, Data *v2, int N, Data *mtx,Data *v3)
{
  Data cent1[3];
  Data cent2[3];
  MATRIX tmtx;
  MATRIX tempmtx;
  MATRIX move1;
  MATRIX move2;
  int i,j;
  Data answer;
  Data *temp1 = 0;
  Data *temp2 = 0;
  int err;

  assert(N > 3);

  temp1 = (Data *)malloc(N * 3 * sizeof(Data));
  temp2 = (Data *)malloc(N * 3 * sizeof(Data));
  if(!temp1 || !temp2)
    goto error_exit;

  centroid(cent1, v1, N);
  centroid(cent2, v2, N);
  for(i=0;i<N;i++)
  {
    temp1[i*3+0] = v1[i*3+0] - cent1[0];
    temp1[i*3+1] = v1[i*3+1] - cent1[1];
    temp1[i*3+2] = v1[i*3+2] - cent1[2];
    
    temp2[i*3+0] = v2[i*3+0] - cent2[0];
    temp2[i*3+1] = v2[i*3+1] - cent2[1];
    temp2[i*3+2] = v2[i*3+2] - cent2[2];
  }

  err = getalignmtx(temp1, temp2, N, &tmtx);
  if(err == -1)
    goto error_exit;
 
  mtx_trans(&move1, -cent2[0], -cent2[1], -cent2[2]);
  mtx_mul(&tempmtx, &move1, &tmtx);
  mtx_trans(&move2, cent1[0], cent1[1], cent1[2]);
  mtx_mul(&tmtx, &tempmtx, &move2);
  memcpy(temp2, v2, N * sizeof(Data) * 3);
  for(i=0;i<N;i++)
  {
    mulpt(&tmtx, temp2 + i * 3);
  }
	for(i=0;i<N;i++)
	{
		for(j=0;j<3;j++)
		{
			v3[i*3+j]=temp2[i*3+j];
		}
	}
  answer = alignedrmsd(v1, temp2, N);
  free(temp1);
  free(temp2);
  if(mtx)
    memcpy(mtx, &tmtx.m, 16 * sizeof(Data));

  return answer;
 error_exit:
  free(temp1);
  free(temp2);
  if(mtx)
  {
    for(i=0;i<16;i++)
      mtx[i] = 0;
  }
  return sqrt(-1.0);
}

/*
  calculate rmsd between two aligned structures (trivial)
  Params: v1 - first structure
          v2 - second structure
          N - number of points
  Returns: rmsd
 */
static Data alignedrmsd(Data *v1, Data *v2, int N)
{
  Data answer =0;
  int i;

  for(i=0;i<N;i++)
    answer += vdiff2(v1 + i *3, v2 + i * 3);
  return sqrt(answer/N);
}

/*
  compute the centroid
 */
static void centroid(Data *ret, Data *v, int N)
{
  int i;

  ret[0] = 0;
  ret[1] = 0;
  ret[2] = 0;
  for(i=0;i<N;i++)
  {
    ret[0] += v[i*3+0];
    ret[1] += v[i*3+1];
    ret[2] += v[i*3+2];
  }
  ret[0] /= N;
  ret[1] /= N;
  ret[2] /= N;
}

/*
  get the matrix needed to align two structures
  Params: v1 - reference structure
          v2 - structure to align
          N - number of points
          mtx - return for rigid body alignment matrix
  Notes: only calculates rotation part of matrix.
         assumes input has been aligned to centroids 
 */
static int getalignmtx(Data *v1, Data *v2, int N, MATRIX *mtx)
{
  MATRIX A = { {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}} };
  MATRIX At;
  MATRIX Ainv;
  MATRIX temp;
  Data tv[3];
  Data tw[3];
  Data tv2[3];
  Data tw2[3];
  int k, i, j;
  int flag = 0;
  Data correction;

  correction = absmaxv(v1, N * 3) * absmaxv(v2, N * 3);
  
  for(k=0;k<N;k++)
    for(i=0;i<3;i++)
      for(j=0;j<3;j++)
	A.m[i][j] += (v1[k*3+i] * v2[k*3+j])/correction;
  
  while(flag < 3)
  {
    for(i=0;i<4;i++)
      for(j=0;j<4;j++)
        At.m[i][j] = A.m[j][i];

    memcpy(&Ainv, &A, sizeof(MATRIX));
    /* this will happen if all points are in a plane */
    if( mtx_invert((Data *) &Ainv, 4) == -1)
    {
      if(flag == 0)
      {
        crossproduct(tv, v1, v1+3);
        crossproduct(tw, v2, v2+3);
      }
      else
      {
        crossproduct(tv2, tv, v1);
        crossproduct(tw2, tw, v2);
        memcpy(tv, tv2, 3 * sizeof(Data));
        memcpy(tw, tw2, 3 * sizeof(Data));
      }
      for(i=0;i<3;i++)
        for(j=0;j<3;j++)
	  A.m[i][j] += tv[i] * tw[j];
      
      flag++;
    }
    else
      flag = 5;
  }
  if(flag != 5)
    return -1;

  mtx_mul(&temp, &At, &A);
  mtx_root(&temp);
  mtx_mul(mtx, &temp, &Ainv); 
  return 0;
}

/*
  get the crossproduct of two vectors.
  Params: ans - return pinter for answer.
          pt1 - first vector
		  pt2 - second vector.
  Notes: crossproduct is at right angles to the two vectors.
*/
static void crossproduct(Data *ans, Data *pt1, Data *pt2)
{
  ans[0] = pt1[1] * pt2[2] - pt1[2] * pt2[1];
  ans[1] = pt1[0] * pt2[2] - pt1[2] * pt2[0];
  ans[2] = pt1[0] * pt2[1] - pt1[1] * pt2[0];
}

/*
  Denman-Beavers square root iteration
 */
static void mtx_root(MATRIX *mtx)
{
  MATRIX Y = *mtx;
  MATRIX Z;
  MATRIX Y1;
  MATRIX Z1;
  MATRIX invY;
  MATRIX invZ;
  MATRIX Y2;
  int iter = 0;
  int i, ii;

  mtx_identity(&Z);

  do
  {
    invY = Y;
    invZ = Z;
    if( mtx_invert((Data *) &invY, 4) == -1)
      return;
    if( mtx_invert((Data *) &invZ, 4) == -1)
      return;
    for(i=0;i<4;i++)
      for(ii=0;ii<4;ii++)
      {
        Y1.m[i][ii] = 0.5 * (Y.m[i][ii] + invZ.m[i][ii]);
	Z1.m[i][ii] = 0.5 * (Z.m[i][ii] + invY.m[i][ii]);
      }
    Y = Y1;
    Z = Z1;

    mtx_mul(&Y2, &Y, &Y);
  }
  while(!almostequal(&Y2, mtx) && iter++ < 20 );

  *mtx = Y;
}

/*
  Check two matrices for near-enough equality
  Params: a - first matrix
          b - second matrix
  Returns: 1 if almost equal, else 0, epsilon 0.0001f.
 */
static int almostequal(MATRIX *a, MATRIX *b)
{
  int i, ii;
  Data epsilon = 0.001f;

  for(i=0;i<4;i++)
    for(ii=0;ii<4;ii++)
      if(fabs(a->m[i][ii] - b->m[i][ii]) > epsilon) 
        return 0;
  return 1;
}  

/*
  multiply a point by a matrix.
  Params: mtx - matrix
          pt - the point (transformed)
*/
void mulpt(MATRIX *mtx, Data *pt)
{
  Data ans[4] = {0};
  int i;
  int ii;

  for(i=0;i<4;i++)
  {
    for(ii=0;ii<3;ii++)
	{
	  ans[i] += pt[ii] * mtx->m[ii][i];
	}
	ans[i] += mtx->m[3][i];
  }
  pt[0] = ans[0];
  pt[1] = ans[1];
  pt[2] = ans[2];
} 

/*
  multiply two matrices.
  Params: ans - return pointer for answer.
          x - first matrix
		  y - second matrix.
  Notes: ans may not be equal to x or y.
*/
static void mtx_mul(MATRIX *ans, MATRIX *x, MATRIX *y)
{
  int i;
  int ii;
  int iii;

  for(i=0;i<4;i++)
    for(ii=0;ii<4;ii++)
	{
	  ans->m[i][ii] = 0;
	  for(iii=0;iii<4;iii++)
	    ans->m[i][ii] += x->m[i][iii] * y->m[iii][ii];
    }
}


/*
  create an identity matrix.
  Params: mtx - return pointer.
*/
static void mtx_identity(MATRIX *mtx)
{
  int i;
  int ii;

  for(i=0;i<4;i++)
    for(ii=0;ii<4;ii++)
	{
	  if(i==ii)
	    mtx->m[i][ii] = 1.0f;
	  else
	    mtx->m[i][ii] = 0;
	}
}

/*
  create a translation matrix.
  Params: mtx - return pointer for matrix.
          x - x translation.
		  y - y translation.
		  z - z translation
*/
static void mtx_trans(MATRIX *mtx, Data x, Data y, Data z)
{
  mtx->m[0][0] = 1;
  mtx->m[0][1] = 0;
  mtx->m[0][2] = 0;
  mtx->m[0][3] = 0;

  mtx->m[1][0] = 0;
  mtx->m[1][1] = 1;
  mtx->m[1][2] = 0;
  mtx->m[1][3] = 0;
  
  mtx->m[2][0] = 0;
  mtx->m[2][1] = 0;
  mtx->m[2][2] = 1;
  mtx->m[2][3] = 0;
  
  mtx->m[3][0] = x;
  mtx->m[3][1] = y;
  mtx->m[3][2] = z;
  mtx->m[3][3] = 1; 
}

/*
   matrix invert routine
  Params: mtx - the matrix in raw format, in/out
          N - width and height
  Returns: 0 on success, -1 on fail
 */
static int mtx_invert(Data *mtx, int N)
{
  int indxc[100]; /* these 100s are the only restriction on matrix size */
  int indxr[100];
  int ipiv[100];
  int i, j, k;
  int irow, icol;
  Data big;
  Data pinv;
  int l, ll;
  Data dum;
  Data temp;
  
  assert(N <= 100);

  for(i=0;i<N;i++)
    ipiv[i] = 0;

  for(i=0;i<N;i++)
  {
    big = 0.0;

    /* find biggest element */
    for(j=0;j<N;j++)
      if(ipiv[j] != 1)
        for(k=0;k<N;k++)
          if(ipiv[k] == 0)
            if(fabs(mtx[j*N+k]) >= big)
	    {
	       big = fabs(mtx[j*N+k]);
               irow = j;
               icol = k;
	    }       
	      
    ipiv[icol]=1;

    if(irow != icol)
      for(l=0;l<N;l++)
      {
        temp = mtx[irow * N + l];
        mtx[irow * N + l] = mtx[icol * N + l];
	mtx[icol * N + l] = temp;
      }

    indxr[i] = irow;
    indxc[i] = icol;

       
    /* if biggest element is zero matrix is singular, bail */
    if(mtx[icol* N + icol] == 0)
       goto error_exit;
		  
    pinv = 1.0/mtx[icol * N + icol];
           
    mtx[icol * N + icol] = 1.0;
 
    for(l=0;l<N;l++)
      mtx[icol * N + l] *= pinv;
                   
    for(ll=0;ll<N;ll++)
      if(ll != icol)
      {
        dum = mtx[ll * N + icol];
        mtx[ll * N + icol] = 0.0;
        for(l=0;l<N;l++) 
          mtx[ll * N + l] -= mtx[icol * N + l]*dum;
      }
  }                
   

  /* unscramble matrix */
  for (l=N-1;l>=0;l--) 
  {
    if (indxr[l] != indxc[l])
    for (k=0;k<N;k++)
    {
      temp = mtx[k * N + indxr[l]];
      mtx[k * N + indxr[l]] = mtx[k * N + indxc[l]];
      mtx[k * N + indxc[l]] = temp;
    }
  } 

  return 0;

 error_exit:
  return -1;
}

/*
  get the asolute maximum of an array
 */
static Data absmaxv(Data *v, int N)
{
  Data answer;
  int i;

  for(i=0;i<N;i++)
    if(answer < fabs(v[i]))
      answer = fabs(v[i]);
  return answer;
}
