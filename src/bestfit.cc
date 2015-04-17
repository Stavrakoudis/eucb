typedef float rotate_t[3][3];
typedef float xlate_t[3];
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
#define DSQR(a) ((double)(a)*(double)(a))

#define SQRTABS(a) ((a) > 0 ? sqrt(a) : sqrt(-a))
#define SQR(a) ((a)*(a))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define FABS(a) ((a) > 0 ? (a) : (-a))

double vsumsqr(const float v[], unsigned int n)
  {
  double sum;
  unsigned int i;
# ifdef _UNROLL
  unsigned int j;
#endif

  sum = DSQR(*v); n--; v++;
# ifdef _UNROLL
  j = n & 7;
  for (i = 0; i != j; i++, v++) sum += DSQR(*v);
  for ( ; i != n; i += 8, v += 8)
    sum += DSQR(v[0]) + DSQR(v[1]) + DSQR(v[2]) + DSQR(v[3]) +
           DSQR(v[4]) + DSQR(v[5]) + DSQR(v[6]) + DSQR(v[7]);
# else
  for (i = 0; i != n; i++, v++) sum += DSQR(*v);
# endif
  return sum;
  }



double vdotby3(const float a[], const float b[], unsigned int n)
  {
  double sum;
  unsigned int i;
# ifdef _UNROLL
  unsigned int j;
#endif

  sum = (double)(*a) * (double)(*b); n--; a += 3; b += 3;
# ifdef _UNROLL
  j = n & 7;
  for (i = 0; i != j; i++, a += 3, b += 3)
    sum += (double)(*a) * (double)(*b);
  for ( ; i != n; i += 8, a += 24, b += 24)
    sum += ((double)a[0] * (double)b[0] +
            (double)a[3] * (double)b[3] +
            (double)a[6] * (double)b[6] +
            (double)a[9] * (double)b[9] +
            (double)a[12] * (double)b[12] +
            (double)a[15] * (double)b[15] +
            (double)a[18] * (double)b[18] +
            (double)a[21] * (double)b[21]);
# else
  for (i = 0; i != n; i++, a += 3, b += 3)
    sum += (double)(*a) * (double)(*b);
# endif
  return sum;
  }

static void hinds_matrix(rotate_t u, xlate_t t,
                   int sw, double e[3], double rr[6], double r[3][3],
                   double sx[3], double sy[3], int npts)
  {
  double ss[6], d, h, p, a[3][3], b[3][3];
  int i, j, l, m, m1, m2, m3;

  if (sw == 4) {
    /* Case of three identical roots */
    for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++)
        if (i == j) a[i][j] = 1.0; else a[i][j] = 0.0;
    }
  else {
    if (sw == 1) {
      /* Case of three distinct roots */
      for (l=0; l < 2; l++) {
        d = e[l];
        ss[0] = (d-rr[2]) * (d-rr[5]) - SQR(rr[4]);
        ss[1] = (d-rr[5]) * rr[1] + rr[3]*rr[4];
        ss[2] = (d-rr[0]) * (d-rr[5]) - SQR(rr[3]);
        ss[3] = (d-rr[2]) * rr[3] + rr[1]*rr[4];
        ss[4] = (d-rr[0]) * rr[4] + rr[1]*rr[3];
        ss[5] = (d-rr[0]) * (d-rr[2]) - SQR(rr[1]);
        if (FABS(ss[0]) >= FABS(ss[2])) {
          if (FABS(ss[0]) >= FABS(ss[5])) {
            a[0][l] = ss[0]; a[1][l] = ss[1]; a[2][l] = ss[3];
            }
          else {
            a[0][l] = ss[3]; a[1][l] = ss[4]; a[2][l] = ss[5];
            }
          }
        else if (FABS(ss[2]) >= FABS(ss[5])) {
          a[0][l] = ss[1]; a[1][l] = ss[2]; a[2][l] = ss[4];
          }
        else {
          a[0][l] = ss[3]; a[1][l] = ss[4]; a[2][l] = ss[5];
          }
        d = sqrt(SQR(a[0][l]) + SQR(a[1][l]) + SQR(a[2][l]));
        a[0][l] /= d; a[1][l] /= d; a[2][l] /= d;
        }
      m1 = 2; m2 = 0; m3 = 1;
      }
    else {
      /* Cases of two distinct roots */
      if (sw == 2) {
        m = 0; m1 = 2; m2 = 0; m3 = 1;
        }
      else {
        m = 2; m1 = 0; m2 = 1; m3 = 2;
        }
      h = e[2];
      a[0][1] = 1.0; a[1][1] = 1.0; a[2][1] = 1.0;
      if (FABS(rr[0]-h) > FABS(rr[2]-h)) {
        if (FABS(rr[0]-h) > FABS(rr[5]-h)) {
          a[0][m] = rr[0]-h; a[1][m] = rr[1]; a[2][m] = rr[3];
          p = -(rr[0] + rr[1] + rr[3] - h);
          a[0][1] = p/a[0][m];
          }
        else {
          a[0][m] = rr[3]; a[1][m] = rr[4]; a[2][m] = rr[5]-h;
          p = -(rr[3] + rr[4] + rr[5] - h);
          a[2][1] = p/a[2][m];
          }
        }
      else {
        if (FABS(rr[2]-h) > FABS(rr[5]-h)) {
          a[0][m] = rr[1]; a[1][m] = rr[2]-h; a[2][m] = rr[4];
          p = -(rr[1] + rr[2] + rr[4] - h);
          a[1][1] = p/a[1][m];
          }
        else {
          a[0][m] = rr[3]; a[1][m] = rr[4]; a[2][m] = rr[5]-h;
          p = -(rr[3] + rr[4] + rr[5] - h);
          a[2][1] = p/a[2][m];
          }
        }
      d = sqrt(SQR(a[0][m]) + SQR(a[1][m]) + SQR(a[2][m]));
      p = sqrt(SQR(a[0][1]) + SQR(a[1][1]) + SQR(a[2][1]));
      for (i = 0; i < 3; i++) {
        a[i][1] /= p;
        a[i][m] /= d;
        }
      }
    /* Common for either two or three distinct roots */
    a[0][m1] = a[1][m2]*a[2][m3] - a[1][m3]*a[2][m2];
    a[1][m1] = a[2][m2]*a[0][m3] - a[2][m3]*a[0][m2];
    a[2][m1] = a[0][m2]*a[1][m3] - a[0][m3]*a[1][m2];
    }
  
  for (l = 0; l < 2; l++) {
    d = 0;
    for (i = 0; i < 3; i++) {
      b[i][l] = r[0][i]*a[0][l] + r[1][i]*a[1][l] + r[2][i]*a[2][l];
      d += SQR(b[i][l]);
      }
    d = sqrt(d);
    if (d == 0) {
    			printf("\033[31m\033[1mRotation matrix determination failed. Abort.\033[0m\n");
    			exit(1);
    					}
    					
    for (i = 0; i < 3; i++) b[i][l] /= d;
    }
  b[0][2] = b[1][0]*b[2][1] - b[1][1]*b[2][0];
  b[1][2] = b[2][0]*b[0][1] - b[2][1]*b[0][0];
  b[2][2] = b[0][0]*b[1][1] - b[0][1]*b[1][0];
  
  /* Calculate rotation matrix */
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      u[i][j] = b[i][0]*a[j][0] + b[i][1]*a[j][1] + b[i][2]*a[j][2];
  /* Calculate translation vector */
  for (i = 0; i < 3; i++)
    t[i] = (sy[i] - u[i][0]*sx[0] - u[i][1]*sx[1] - u[i][2]*sx[2])/npts;
  }

static int solve(double e[3],
                   double e0, double det, double spur, double cof)
  {
  double d, h, g, cth, sth, sqrth;
#define epsilon 1.0e-37
#define sqrt3 1.7320508075688772

  det *= det;
  d = SQR(spur);
  h = d - cof;
  g = spur * (cof*1.5 - d) - det*0.5;
  if (h > d*epsilon) {
    sqrth = SQRTABS(h);
    d = -g/(h*sqrth);
    if (d > 1-epsilon) {    /* Two identical roots */
      e[0] = spur + 2*sqrth;
      e[2] = e[1] = spur - sqrth;
      return 2;
      }
    else if (d < -1+epsilon) {    /* Two identical roots */
      e[1] = e[0] = spur+sqrth;
      e[2] = MAX(0, spur - 2*sqrth);
      return 3;
      }
    else {    /* Three distinct roots */
      d = acos(d)/3.0;
      cth = sqrth*cos(d); sth = sqrth*sqrt3*sin(d);
      e[0] = spur + 2*cth;
      e[1] = spur - cth + sth;
      e[2] = MAX(0, spur - cth - sth);
      return 1;
      }
    }
else {    /* Three identical roots */
    e[0] = e[1] = e[2] = spur;
    return 4;
    }

  } /* solve */


double vsumby3(const float v[], unsigned int n)
  {
  double sum;
  unsigned int i;
# ifdef _UNROLL
  unsigned int j;
#endif

  sum = *v; n--; v += 3;
# ifdef _UNROLL
  j = n & 7;
  for (i = 0; i != j; i++, v += 3) sum += *v;
  for ( ; i != n; i += 8, v += 24)
    sum += (double)v[0] + (double)v[3] +
           (double)v[6] + (double)v[9] +
           (double)v[12] + (double)v[15] +
           (double)v[18] + (double)v[21];
#else
  for (i = 0; i != n; i++, v += 3) sum += *v;
#endif
  return sum;
  }

float bestfit(rotate_t u, xlate_t t,
                float *a, float *b, int npts)
  {
  double sx[3], sy[3], sxy[3][3], sx2, sy2;
  double r[3][3], rr[6];
  double d, spur, det, cof, e0, e[3], inpts;
  int i, j, m, sw;
  float *aa=0, *bb=0;
 sx[0] = vsumby3(&a[0+0], npts);
  sx[1] = vsumby3(&a[0+1], npts);
  sx[2] = vsumby3(&a[0+2], npts);
  sy[0] = vsumby3(&b[0+0], npts);
  sy[1] = vsumby3(&b[0+1], npts);
  sy[2] = vsumby3(&b[0+2], npts);
  sx2 = vsumsqr(&a[0], 3*npts);
  sy2 = vsumsqr(&b[0], 3*npts);
  sxy[0][0] = vdotby3(&a[0], &b[0], npts);
  sxy[0][1] = vdotby3(&a[0], &b[1], npts);
  sxy[0][2] = vdotby3(&a[0], &b[2], npts);
  sxy[1][0] = vdotby3(&a[1], &b[0], npts);
  sxy[1][1] = vdotby3(&a[1], &b[1], npts);
  sxy[1][2] = vdotby3(&a[1], &b[2], npts);
  sxy[2][0] = vdotby3(&a[2], &b[0], npts);
  sxy[2][1] = vdotby3(&a[2], &b[1], npts);
  sxy[2][2] = vdotby3(&a[2], &b[2], npts);
inpts = 1.0/npts;
  e0 = (sx2 - (SQR(sx[0])+SQR(sx[1])+SQR(sx[2]))*inpts +
        sy2 - (SQR(sy[0])+SQR(sy[1])+SQR(sy[2]))*inpts)*inpts;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      r[i][j] = (sxy[i][j] - sx[i]*sy[j]*inpts)*inpts;
      }
    }

  for (m = j = 0; j < 3; j++) {
    for (i = 0; i <= j; i++,m++) {
      rr[m] = r[i][0]*r[j][0] + r[i][1]*r[j][1] + r[i][2]*r[j][2];
      }
    }

  det = r[0][0] * (r[1][1]*r[2][2] - r[2][1]*r[1][2])
      - r[1][0] * (r[0][1]*r[2][2] - r[2][1]*r[0][2])
      + r[2][0] * (r[0][1]*r[1][2] - r[1][1]*r[0][2]);

  spur = (rr[0]+rr[2]+rr[5])/3.0;
  cof = (rr[2]*rr[5] - SQR(rr[4]) + rr[0]*rr[5] -
         SQR(rr[3]) + rr[0]*rr[2] - SQR(rr[1]))/3.0;


  sw = solve(e, e0, det, spur, cof);

        if ( sw != 1 )
                {
                printf("\033[31m\033[1mInternal problem with bestfit(). Please report.\033[0m\n");
                exit( 1 );
                }

  hinds_matrix(u, t, sw, e, rr, r, sx, sy, npts);
  d = SQRTABS(e[2]);
  if (det < 0) d = -d;
  d += SQRTABS(e[1]) + SQRTABS(e[0]);
  d = e0-2*d;

  return SQRTABS(d);


}
