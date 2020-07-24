int tov_r(double rhoc, double R, int *npts, 
	  double **p_r, double **p_m,
	  double **p_rho, double **p_pre,double **p_phi, 
	  double **p_riso);

void grhdc_p2c_point(
  double rho, double eps, double v1, double v2, double v3, double p, double psi4,
  double *D, double *tau, double *S1, double *S2, double *S3);

void    interp_lag4(double *f, double *x, int Nx, double xv, 
                    double *fv_p, double *dfv_p, double *ddfv_p );

void eos_simple(double rho, double *p, double *dpdrho, double *epsl);
