/* numerical stuff */

#include <math.h>
#include <stdio.h>
#include "numerical.h"
#include "nrutil.h"

/* Numerical Recipes Stuff */

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY  1.0e-20
#define MAX(a,b) ((a) > (b)  ? (a) : (b))
#define SIGN(a,b)  ((b)  > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d)  (a)=(b);(b)=(c);(c)=(d);
#define FUNC(x) ((*func)(x))
#define FUNC2(x,p) ((*func)(x,p))
 
#define EPS 1.0e-5
#define JMAX 50


/* this is from numerical recipes; brackets a minimum with ax, bx, cx */
void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)()){
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if(*fb > *fa){
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx = *bx + GOLD * (*bx-*ax);
	*fc=(*func)(*cx);
	while(*fb > *fc){
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u = *bx-((*bx-*cx)*q-(*bx-*ax)*r) /
		    (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim = *bx + GLIMIT * (*cx-*bx);
		if((*bx-u) * (u-*cx) > 0.0){
			fu=(*func)(u);
			if(fu < *fc){
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				return;
   	 			}
			else if(fu > *fb){
				*cx=u;
				*fc=fu;
				return;
				}
			u = *cx + GOLD * (*cx-*bx);
			fu = (*func)(u);
			}
		else if((*cx-u)*(u-ulim) > 0.0){
			fu=(*func)(u);
			if(fu < *fc){
				SHFT(*bx,*cx,u,*cx + GOLD * (*cx - *bx));
				SHFT(*fb,*fc,fu,(*func)(u));
			}
		}
		else if((u-ulim)*(ulim-*cx) >= 0.0){
			u=ulim;
			fu=(*func)(u);
		}
		else{
			u = *cx+GOLD*(*cx-*bx);
			fu=(*func)(u);
			}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
		}
	return;
}

/* this is from numerical recipes; brackets a minimum with ax, bx, cx */
void mnbrak2(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double, void * ), double beta, void * p){
	double ulim,u,r,q,fu,dum;
	

	*fa=(*func)(*ax, p);
	*fb=(*func)(*bx, p);
	if(*fb > *fa){
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx = *bx + GOLD * (*bx-*ax);
	*fc=(*func)(*cx, p);
	while(*fb > *fc){
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u = *bx-((*bx-*cx)*q-(*bx-*ax)*r) /
		    (2.*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim = *bx + GLIMIT * (*cx-*bx);
		if((*bx-u) * (u-*cx) > 0.0){
			fu=(*func)(u, p);
			if(fu < *fc){
				*ax = *bx;
				*bx = u;
				*fa = *fb;
				*fb = fu;
				return;
   	 			}
			else if(fu > *fb){
				*cx=u;
				*fc=fu;
				return;
				}
			u = *cx + GOLD * (*cx-*bx);
			fu = (*func)(u, p);
			}
		else if((*cx-u)*(u-ulim) > 0.0){
			fu=(*func)(u, p);
			if(fu < *fc){
				SHFT(*bx,*cx,u,*cx + GOLD * (*cx - *bx));
				SHFT(*fb,*fc,fu,(*func)(u, p));
			}
		}
		else if((u-ulim)*(ulim-*cx) >= 0.0){
			u=ulim;
			fu=(*func)(u, p);
		}
		else{
			u = *cx+GOLD*(*cx-*bx);
			fu=(*func)(u, p);
			}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
		}
	return;
}


/* Numerical Integration Routines */


#define FUNC(x) ((*func)(x))
double d_trapzd(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

/* specific version for GSL function */
double d_trapzd2(double (*func)(double, void *), double a, double b, int n, void * p)
{
	double x,tnm,sum,del;
	static double s;
	int it,j;

	if (n == 1) {
	  return (s=0.5*(b-a)*(FUNC2(a,p)+FUNC2(b,p)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC2(x,p);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}

double d_qtrap(double (*func)(double), double a, double b)
{
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double s,olds;
	int j;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=d_trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}

/* specific version for GSL function */
double d_qtrap2(double (*func)(double, void *), double a, double b, void * p){
  double s,olds;
  int j;

  olds = -1.0e30;
  for (j=1;j<=JMAX;j++) {
    s=d_trapzd2(func,a,b,j,p);
    if (fabs(s-olds) < EPS*fabs(olds)) return s;
    olds=s;
  }
  nrerror("Too many steps in routine qtrap");
  return 0.0;
}

double d_qsimp(double (*func)(double), double a, double b)
{
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	int j;
	double s,st,ost,os;

	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=d_trapzd(func,a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) return s;
		os=s;
		ost=st;
	}
	nrerror("Too many steps in routine qsimp");
	return 0.0;
}

#define JMAXP (JMAX+1)
#define K 5
double d_qromb(double (*func)(double), double a, double b)
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	double trapzd(double (*func)(double), double a, double b, int n);
	void nrerror(char error_text[]);
	double ss,dss;
	double s[JMAXP+1],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=d_trapzd(func,a,b,j);
		if (j >= K) {
			d_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=0.25*h[j];
	}
	nrerror("Too many steps in routine qromb");
	return 0.0;
}

double d_qromb2(double (*func)(double, void *), double a, double b, void * p){

  double ss,dss;
  double s[JMAXP+1],h[JMAXP+1];
  int j;

  h[1]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j]=midpnt2(func,a,b,j,p);
    if (j >= K) {
      d_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
      if (fabs(dss) < EPS*fabs(ss)) return ss;
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  nrerror("Too many steps in routine qromb");
  return 0.0;
}
#undef JMAXP
#undef K

#define NRANSI
#include "nrutil.h"
void d_polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=dvector(1,n);
	d=dvector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_dvector(d,1,n);
	free_dvector(c,1,n);
}
#undef NRANSI



/* from chp. 4 of numerical recipes */
double midpnt(double (*func)(double), double a, double b, int n){
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC(0.5*(a+b)));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC(x);
      x += ddel;
      sum += FUNC(x);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}

/* from chp. 4 of numerical recipes, specific for GSL functions */
double midpnt2(double (*func)(double, void *), double a, double b, int n, void * p){
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;
  
  if (n == 1) {
    return (s=(b-a)*FUNC2(0.5*(a+b),p));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += FUNC2(x,p);
      x += ddel;
      sum += FUNC2(x,p);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}




