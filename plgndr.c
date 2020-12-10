#include <math.h>
#include <stdio.h>
#include "plgndr.h"
#include <stdlib.h>
double plgndr(int l, int m, double x)
{
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;
  
  if (m < 0 || m > l || fabs(x) > 1.0) {
    fprintf(stderr, "Bad arguments in routine plgndr()\n");
    fprintf(stderr, "l: %d m: %d x: %g\n", l, m, x);
    exit(1);
  }
  pmm=1.0;
  if (m > 0) {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=m+2;ll<=l;ll++) {
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;
      }
      return pll;
    }
  }
}
double sphericalY(int l, int m, double theta)
{
  int i;
  double lmm = 1.0, lpm;

  for (i = 2; i <= l - m; i++)
    lmm *= (double)i;
  lpm = lmm;
  for (i = l - m + 1; i <= l + m; i++)
    lpm *= (double)i;
  return sqrt((2.*l + 1.)*lmm/(4.*M_PI*lpm))*plgndr(l, m, cos(theta));
}
