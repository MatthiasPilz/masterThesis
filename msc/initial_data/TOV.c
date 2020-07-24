/* TOV equations solver
   dtim import Sebastianos TOV-solver from bam */
/* at this point first use only polytropic EOS and use the integration over r*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI  3.1415926535897932

#define PR 0
// PR 1 some output PR 2 all output

double EOS_GAMMA;
double EOS_GAMMAMO;
double EOS_KAPPA;

/* init simpelest possible EOS*/ 
void eos_init() {	
  EOS_GAMMA   = 2.;
  EOS_GAMMAMO = EOS_GAMMA - 1.;
  EOS_KAPPA   = 100.;
}

void eos_simple(double rho, double *p, double *dpdrho, double *epsl) {
	
  *p      = EOS_KAPPA * pow(rho,EOS_GAMMA);
  *dpdrho = EOS_KAPPA * EOS_GAMMA *  pow(rho,EOS_GAMMAMO);
  *epsl   = EOS_KAPPA/EOS_GAMMAMO * pow(rho,EOS_GAMMAMO);
	
	}


/* rest-mass density calculation */
double computeMb( double *r, double *m, double *rho, 
		  int n )
{
  int i;
  double sum = 0.;
  double dr, tmp;
  for (i=1; i<n; i++)  {
    dr    = r[i]-r[i-1];
    tmp   = 1./sqrt(1.-2.*m[i]/r[i]);
    sum  += r[i]*r[i]*rho[i]*tmp*dr;
  }
  return 4.*PI*sum;
}

#define nvar 5
#define forv for (v=0; v<nvar; v++)


int     interp_locate(double *x, int Nx, double xval)
{
  int ju,jm,jl;
  int ascnd;
  
  jl=-1;
  ju=Nx;
  
  if (xval <= x[0]) { 
    if (xval < x[0]) if (PR) printf("  pt to locate is outside (xval<xx).\n"); 
    return 0; 
  } else if (xval >= x[Nx-1]) { 
    if (xval > x[Nx-1])if (PR)  printf("  pt to locate is outside (xval>xx).\n"); 
    return Nx-1; 
  }
  
  ascnd = (x[Nx-1] >= x[0]);
  
  while (ju-jl > 1) {
    
    jm = (ju+jl) >> 1;
    
    if (xval >= x[jm] == ascnd)
      jl=jm;
    else
      ju=jm;
    
  }
  
  return jl;
}

void    interp_lag4(double *f, double *x, int Nx, double xv, 
                    double *fv_p, double *dfv_p, double *ddfv_p )
{
  /* Given the values in xv, it returns the interpolated values fv and 
  its 1st and 2nd derivatives dfv, ddfv of the fuction f(x) 
  Lagrangian 4 pts interpolation is used */
  
  if (Nx < 4) { printf(" too few points for interpolation"); exit(0); }
  
  int i = interp_locate(x,Nx,xv);
  
  if( i < 1 ){ 
   if(PR) printf(" too few points on the left => interpolation maybe be inaccurate! (v=%e)\n",xv);
    i = 1;
  }
  if( i > (Nx-3) ){ 
   if(PR) printf(" too few points on the right => interpolation maybe be inaccurate! (v=%e   -> %e %e)\n",xv, x[Nx-2],x[Nx-1]);
    i = Nx-3;
  }
  
  double ximo =  x[i-1];
  double xi   =  x[i];
  double xipo =  x[i+1]; 
  double xipt =  x[i+2]; 
  
  double C1   = (f[i] - f[i-1])/(xi - ximo);
  double C2   = (-f[i] + f[i+1])/(-xi + xipo);
  double C3   = (-f[i+1] + f[i+2])/(-xipo + xipt);
  double CC1  = (-C1 + C2)/(-ximo + xipo);
  double CC2  = (-C2 + C3)/(-xi + xipt);
  double CCC1 = (-CC1 + CC2)/(-ximo + xipt);
  
  *fv_p   = f[i-1] + (-ximo + xv)*(C1 + (-xi + xv)*(CC1 + CCC1*(-xipo + xv)));
  *dfv_p  = C1 - (CC1 - CCC1*(xi + xipo - 2.*xv))*(ximo - xv)
      + (-xi + xv)*(CC1 + CCC1*(-xipo + xv));
  *ddfv_p = 2.*(CC1 - CCC1*(xi + ximo + xipo - 3.*xv));
}

/* tov r rhs */
int tov_r_rhs(double dr, double *u, double *k)
{
  double r   = u[0];
  double rho = u[1];
  double m   = u[2];
  double phi = u[3];
  double I   = u[4];
  
  double p,epsl,dpdrho,e; 
  
  /* simple polytropic EOS*/ 
  eos_simple(rho, &p, &dpdrho, &epsl);
  e      = rho*(1.+epsl);
  
  if (r==0) r=1e-10;
  double tmp1   = m+4.*PI*r*r*r*p;
  double tmp2   = r*r*(1.-2.*m/r);
  double tmp    = (r==0.)?0.:tmp1/tmp2;
  
  double drhodr = -(e+p) * tmp / dpdrho;
  double dmdr   = 4.*PI*r*r*e;
  double dphidr = tmp;
  double f      = sqrt(1.-2.*m/r);
  double dIdr   = ( 1.-f )/( r*f );

  k[0] = 0;
  k[1] = drhodr;
  k[2] = dmdr;
  k[3] = dphidr;
  k[4] = dIdr;
  
  return 0;
}


/* tov r solver */
int tov_r(double rhoc, double R, int *npts, 
	  double **p_r, double **p_m,
	  double **p_rho, double **p_pre,double **p_phi, 
	  double **p_riso)
{

  
  int n,v,i;
  
  double **u = (double **) malloc((*npts)*sizeof(double*));
  for (n=0; n<*npts; n++)
    u[n] = (double *) malloc((nvar)*sizeof(double));
  double u1[nvar],u2[nvar],u3[nvar],k[nvar];
  double fact = 1./6.;
  
  double stp = R/(*npts); 
  double pc, dpdrhoc, epslc, ec;
  eos_init();
  eos_simple(rhoc, &pc, &dpdrhoc, &epslc);
  ec    = rhoc*(1.+epslc);
  
  // get into
  u[0][0] = 0;
  u[0][1] = rhoc;
  u[0][2] = 0;
  u[0][3] = 0;
  u[0][4] = 0;
  
  
  if(PR) {
  printf("tov_r: solve TOV star (in this cell):\n");
  printf("    drho = %.16e npts = %d\n",stp, *npts);
  printf("    rhoc = %.16e\n",rhoc);
  printf("    ec   = %.16e\n",ec);
  printf("    pc   = %.16e\n",pc);
  }
  double rhoo = u[0][1];
  int stop = 0;
  n=0;
  while (u[n][1]>0. && u[n][1]<=1.01*rhoo && stop==0) {

    stop += tov_r_rhs(stp, u[n], k);
    
    // u_1 = u + dt/2 k
    forv u1[v] = u[n][v] + 0.5*stp*k[v]; 
    
    // r = rhs(u_1)
    stop += tov_r_rhs(stp, u1, k);
  
    // u_2 = u + dt/2 k
    forv u2[v] = u[n][v] + 0.5*stp*k[v]; 

    // r = rhs(u_2)
    stop += tov_r_rhs(stp, u2, k);

    // u_3 = u + dt k
    forv u3[v] = u[n][v] + stp*k[v];

    // r = rhs(u_3)
    stop += tov_r_rhs(stp, u3, k);
  
    // u = 1/6 ( -2 u + 2 u_1 + 4 u_2 + 2 u_3 + dt k ) 
    forv u[n+1][v] = fact*( 2.*( - u[n][v] + u1[v] + u3[v] ) + 4.*u2[v] + stp*k[v] );
    
    u[n+1][0] += stp;
    
    if (PR>1) printf("%d %.16e %.16e %.16e %.16e %.16e    %e %d\n",n,u[n+1][0],u[n+1][1],u[n+1][2],u[n+1][3],u[n+1][4],   rhoo,nvar);

    if (n>=(*npts)-5) {
      u = (double **) realloc(u,(*npts*2)*sizeof(double*));
      for (i=*npts; i<*npts*2; i++)
        u[i] = (double *) malloc((nvar)*sizeof(double));
      *npts = (*npts)*2;
    }
    rhoo = u[n][1];
    
    n++;
  }
  
  double p,dpdrho,epsl,phi,C, M,Mb, IR,phiR,phiRa;
  
  *npts = n;
  R     = u[*npts-1][0];
  M     = u[*npts-1][2];
  phiR  = u[*npts-1][3];
  IR    = u[*npts-1][4];
  phiRa = 0.5*log(1.-2.*M/R);
  C     = 1/(2*R) * (sqrt(R*R-2*M*R)+R-M) * exp(-IR);
  
  *p_r   = (double*) malloc (*npts*sizeof(double));
  *p_m   = (double*) malloc (*npts*sizeof(double));
  *p_rho = (double*) malloc (*npts*sizeof(double));
  *p_pre = (double*) malloc (*npts*sizeof(double));
  *p_phi = (double*) malloc (*npts*sizeof(double));
  *p_riso= (double*) malloc (*npts*sizeof(double));
  
  
  for (n=0; n<*npts; n++) {
    
    (*p_r)[n]   = u[n][0]; 
    (*p_rho)[n] = u[n][1];
    (*p_m)[n]   = u[n][2];
    (*p_phi)[n] = (u[n][3]-phiR + phiRa);
    
    eos_simple((*p_rho)[n], &p,&dpdrho,&epsl);
    
    (*p_pre)[n] = p;
    (*p_riso)[n]= (*p_r)[n] * C * exp(u[n][4]);
    
  }

  Mb = computeMb( *p_r, *p_m, *p_rho, *npts-1 );
  if (PR) {
  printf("    R    = %.16e   (%.16e)\n",R,(*p_riso)[*npts-1]);
  printf("    M    = %.16e\n",M);
  printf("    Mb   = %.16e\n",Mb);
  }
  
  // free memory
  for (i=0; i<*npts; i++)
    free(u[i]);
  free(u);
  
  return 0;
}

void grhdc_p2c_point(
  double rho, double eps, double v1, double v2, double v3, double p, double psi4,
  double *D, double *tau, double *S1, double *S2, double *S3)
{
  double W, WW, h, gam, gWWrh;
  
  W     = 1/sqrt(1 - psi4*v1*v1 - psi4*v2*v2 - psi4*v3*v3);
  WW    = W*W;
  h     = 1 + eps + p/rho;
  gam   = pow(psi4,1.5);
  gWWrh = gam * rho *WW * h;
  
  *D   = gam * rho * W;
  *tau = gWWrh - gam * (  p + rho * W);
  *S1  = gWWrh * psi4 * v1;
  *S2  = gWWrh * psi4 * v2; 
  *S3  = gWWrh * psi4 * v3; 
  
}











