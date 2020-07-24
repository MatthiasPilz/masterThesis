#include "ader.h"
#include "grhdc.h"

#define Power(x,y) ( y==1.5?(sqrt((double)(x))*((double)(x))):exp(((double) (y))*log((double) (x))) )
#define PR 1
#define DEBUG_GRHDC 0
#define BIG 1e+15

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))


enum {
    C2P_COLD,
    C2P_ATM,
    C2P_OK,
    C2P_FAIL,
    C2P_NR
};

int grhdc_sym_flags[2*9] = { 1, -1,  1, 1, 1, 1, -1, 1, 1,
                            -1,  1, -1, 1, 1, 1,  1, 1, 1 };

void grhdc_Var_Num ( PARA *par ) {

    par->var = 9;
    par->sym_flags = grhdc_sym_flags;

    printf ( "grhdc equation needs %d variables\n", par->var );

}

// ### ### ### ### ### ### ###

static void create_initial_data_grhdc ( PARA *par ) {

    char command[255];
    sprintf ( command, "../msc/initial_data/run_average %d %lf %lf %s %s", par->N, par->xMin, par->xRange, par->project, par->shape );
    if ( !system ( command ) )
        printf ( "# SUCCESS # - created initial data.\n" );

}

static void set_initial_data_grhdc ( double **Avg, PARA *par ) {

    ASSERT ( Avg );
    ASSERT ( par );

    int i = 0;
    int check;
    double temp[5];
    double unnecessary;
    char filename[255];

// first cell averages
    sprintf ( filename, "../math/initial_data/%s/%s/%d_ICA.txt", par->project, par->shape, par->N );
    FILE *data = fopen ( filename, "r" );   // open file for reading
    ASSERT ( data );

    OUTPUT_EVO ( "opened cell average data" )

    for ( i = 0; i < par->N; ++i ) {
        check = fscanf ( data, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &unnecessary, temp+0, temp+1, temp+2, temp+3, temp+4, temp+5, temp+6, temp+7, temp+8);
        ASSERT ( check );
        Avg[0][i] = temp[0];
        Avg[1][i] = temp[1];
        Avg[2][i] = temp[2];
        Avg[3][i] = temp[3];
        Avg[4][i] = temp[4];
        Avg[5][i] = temp[5];
        Avg[6][i] = temp[6];
        Avg[7][i] = temp[7];
        Avg[8][i] = temp[8];
    }
    fclose ( data );			// close file

}

void grhdc_Init_Data	( FV_DATA *V, PARA *par ) {

    ASSERT ( V );
    ASSERT ( par );

    create_initial_data_grhdc( par );
    set_initial_data_grhdc ( V->Avg, par );

    printf ( "# SUCCESS # - initialized grhdc FV_DATA.\n" );

}

// #### #### ####

/*****************************************************************************/
/* general routine for sound speed for p,rho,epsl */
double eos_cs2_rep(double rho, double epsl, double p, double dpdrho, double dpdepsl)
{
  if (rho==0.)
    return 0.;
  else
    return( ( dpdrho + p*dpdepsl/(rho*rho) )/( 1.0 + epsl + p/rho ) );
}

/*****************************************************************************/
/* ideal gas */
int eos_ideal(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl, double EOSGAMMA)
{

  *p      = (EOSGAMMA-1.)*(*rho)*(*epsl);
  *dpdepsl= (EOSGAMMA-1.)*(*rho);
  *dpdrho = (EOSGAMMA-1.)*(*epsl);
  *cs2    = eos_cs2_rep(*rho,*epsl, *p,*dpdrho,*dpdepsl);

  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1;
  return 0;
}

/*****************************************************************************/
/* polytropic EoS */
int eos_poly(double *p, double *cs2, double *dpdrho, double *dpdepsl, double *rho, double *epsl, double EOSGAMMA, double EOSKAPPA)
{

  double krhotog1 = EOSKAPPA * pow( *rho, EOSGAMMA-1.); // compute less pow !

  *p      = *rho * krhotog1;
  *dpdepsl= 0.;
  *dpdrho = (EOSGAMMA) * krhotog1;
  if (*epsl<=-1.) {
    // used when this eos is the cold part of ideal gas
    *epsl = krhotog1 /(EOSGAMMA-1.);
  }
  *cs2    = eos_cs2_rep(*rho,*epsl, *p,*dpdrho,*dpdepsl);

  if (!finite(*rho) || !finite(*epsl) || !finite(*cs2)) return 1;
  return 0;
}

// c2p using full EoS
int try_hot_part(double conD, double conS2, double conT,
                 double *rho, double *epsl, double *pres,
                 double init, double errormax,double countmax,
                 double patm, double rhoatm, double atmfac, double EOSGAMMA)
{
    int count;

    double temp,tempSQRTarg,tempSQRT, Wlor;
    double pold,peos,pmin;
    double f,df,drhodp,depsldp,kappa,chi,cs2,v2;

    double problem,error,lostdigits;

    // set minimal p
    pmin = sqrt(conS2)-conD-conT;
    pmin = DMAX(pmin,patm);

    // guess pold from previous step
    pold = DMAX(init,pmin);
    if (!finite(pold)) {
        if (PR) printf("starting guess for p is not finite");
        return C2P_FAIL;
    }

    // compute values with first guess
    temp        = conT+pold+conD;
    tempSQRTarg = temp*temp-conS2;
    tempSQRT    = sqrt(tempSQRTarg);
    Wlor        = temp/tempSQRT;
    *rho        = conD/Wlor;
    *epsl       = (tempSQRT - Wlor*pold)/conD - 1.;


    // check if cons vars have physical values
    if (tempSQRTarg<=0.0)
        return C2P_COLD;


    // start Newton Raphson and find pres
    for (count=0,error=2*errormax; count<=countmax; count++) {

        // EOS call
        problem = eos_ideal(&peos,&cs2,&chi,&kappa, rho,epsl, EOSGAMMA);
        //problem = EOS.comp("re","","","pc","re","",
        //                   *rho,*epsl, &peos,&cs2, &chi,&kappa);

        f       = pold - peos;
        drhodp  = conD*conS2/(tempSQRT *temp*temp);
        depsldp = pold*conS2/(conD*tempSQRTarg*tempSQRT);
        df      = 1. - chi*drhodp - kappa*depsldp;

        *pres   = DMAX( pmin, pold-f/df );
        error   = fabs( 1.-pold/(*pres) );
        pold    = *pres;

        temp        = conT+pold+conD;
        tempSQRTarg = temp*temp-conS2;
        tempSQRT    = sqrt(tempSQRTarg);
        Wlor        = temp/tempSQRT;
        *rho        = conD/Wlor;
        *epsl       = (tempSQRT - Wlor*pold)/conD - 1.;

        // go out if something
        if (*rho<atmfac*rhoatm) return C2P_ATM;
        if (problem) return C2P_ATM;
        if (df==0.) return C2P_NR;
        if (error<errormax) break;

        if ((DEBUG_GRHDC) && (count>2)) {
            printf("    %d   %+2.16e   %+2.16e %+2.16e %2.16e %+2.16e %2.16e\n", count,error,df,f,*pres,*epsl,*rho);
            printf("         => %+2.16e %+2.16e\n", chi*drhodp,kappa*depsldp);
            printf("         => %+2.16e %+2.16e\n",  pold,f/df);
            printf("         => %+2.16e %+2.16e %+2.16e\n", conT,pold,conD);
            printf("         => %+2.16e %+2.16e\n",  temp*temp,conS2);
            printf("         => %+2.16e %+2.16e\n",  tempSQRT,Wlor*pold);
            printf("         => %+2.16e %+2.16e\n",  (tempSQRT - Wlor*pold)/conD,1.);
        }

    }


    if (count>=countmax) {
        lostdigits = 0;
        lostdigits = DMAX(lostdigits , fabs(log10(conD)-log10(conT)));
        lostdigits = DMAX(lostdigits , fabs(log10(conD)-log10(*pres)));
        lostdigits = DMAX(lostdigits , fabs(log10(conT)-log10(*pres)));
        lostdigits = DMAX(lostdigits , fabs(log10(tempSQRT)-log10(Wlor*pold)));

        if ( lostdigits  <  log10(error)-log10(1e-15) || error > 1e-4 ) {
            if (PR) {
                printf("too many Newton-Raphson steps (hot)\n");
                printf("  %d  %e %e    %e %e %e\n",count,error, errormax, f/df, f,df);
                printf("  %+e %+e %+e (values  rho,epsl,p)\n",*rho,*epsl,*pres);
                printf("  %+e %+e %+e (D,S2,T)\n",conD,conS2,conT);
            }
            return C2P_NR;
        }
        return C2P_NR;
    }


    if (!finite(*pres) || !finite(*rho) || !finite(*epsl)) return C2P_COLD;
    if (*epsl<0.) return C2P_COLD;

    // test velocity
    v2 = conS2/((conT+conD+*pres)*(conT+conD+*pres));
    if (v2>0.999)
        return C2P_ATM;

    return C2P_OK;
}

// c2p using C2P_COLD part soley
int try_cold_part(double conD, double conS2, double *conT,
                  double *rho, double *epsl, double *pres,
                  double init, double errormax,double countmax,
                  double patm, double rhoatm, double atmfac, double EOSGAMMA, double EOSKAPPA)
{
    int count;

    double temp,tempSQRTarg,tempSQRT, Wlor;
    double rhomin,rhold;
    double g,dg,h,dhdrho,chi, v2;
    double problem,error;
    double cs2,dpdrho,dpdepsl;

    // check if cons vars have physical values
    if (conD <= 0.)
        return C2P_ATM;

    rhomin = rhoatm;
    rhold  = DMAX(init,rhomin);

    if (!finite(rhold)) {
        if (PR) printf("starting guess for rho is not finite");
        return C2P_FAIL;
    }

    // start Newton Raphson and find rho
    for (count=0,error=2*errormax; count<=countmax; count++) {

      // EOS call
      *epsl = -10.;
      problem = eos_poly(pres,&cs2,&dpdrho,&dpdepsl,rho,epsl, EOSGAMMA, EOSKAPPA);
      //EOS.comp("r","","","pe","r","", rhold, pres,epsl, &chi);
      h       = 1.0 + *epsl + *pres/rhold;
      dhdrho  = chi/rhold;

      if (conD*h==0.)       tempSQRTarg = 1.0 + conS2*BIG;
        else                tempSQRTarg = 1.0 + conS2/((conD*h)*(conD*h));
        Wlor    = sqrt(tempSQRTarg);

        g       = Wlor*rhold - conD;
        dg      = Wlor - (rhold*conS2*dhdrho)/(Wlor*conD*conD*h*h*h);

        if ((!finite(dg)) || (dg==0.0)) {
      if (PR) ;//printf("dg==0\n");
      return C2P_NR;

        }

        *rho    = DMAX( rhomin, rhold-g/dg );
        error   = fabs( 1.-rhold/(*rho) );
        rhold   = *rho;

        if (error<errormax) break;
    }

    if (count>=countmax) {
        if (PR) {
            printf("too many Newton-Raphson steps (C2P_COLD)\n");
            printf("  %d  %e %e    %e %e %e\n",count,error, errormax, g/dg, g,dg);
            printf("  %+e %+e %+e (values  rho,epsl,p)\n",*rho,*epsl,*pres);
            printf("  %+e %+e %+e (D,S2,T)\n",conD,conS2,*conT);
        }
        return C2P_NR;
    }


    // EOS call
    *epsl = -10.;
    problem = eos_poly(pres,&cs2,&dpdrho,&dpdepsl,rho,epsl, EOSGAMMA, EOSKAPPA);
    //EOS.comp("r","","","pe","r","", rhold, pres,epsl, &chi);
    h       = 1.0 + *epsl + *pres/rhold;
    dhdrho  = chi/rhold;

    if (conD*h==0.)     tempSQRTarg = 1.0 + conS2*BIG;
    else                tempSQRTarg = 1.0 + conS2/((conD*h)*(conD*h));
    Wlor    = sqrt(tempSQRTarg);

    // test for atm
    if (*rho<atmfac*rhoatm)
        return C2P_ATM;

    // update tau
    *conT = Wlor*Wlor*(*rho)*h - *pres - conD;

    // test velocity
    v2 = conS2/((*conT+conD+*pres)*(*conT+conD+*pres));
    if (v2>0.999)
        return C2P_ATM;

    return C2P_OK;
}

void grhdc_c2p(double *gD, double *gS, double *gtau, double *gpsi4, double *galpha, double *grho, double *gv, double *geps, double *gp, double EOSGAMMA) {

    // local variables
    double chi;
    double d;
    double depsldp;
    double df;
    double drhodp;
    double eps;
    double err;
    double f;
    double h;
    double kappa;
    double p;
    double peos;
    double pmin;
    double pold;
    double rho;
    double S1;
    double S2;
    double S3;
    double Sinv1;
    double Sinv2;
    double Sinv3;
    double SS;
    double tau;
    double temp;
    double tempSQRT;
    double tempSQRTarg;
    double W;


    const int itermax    = 20;
    const double errgoal = 1e-12;
    const double errmax  = 1e20;
    double gamma    = EOSGAMMA;
    double EOSkappa = 100.0;

    double rhoc     = 0.00128;
    double rhoatm   = 1e-8 * rhoc;
    double rhoatmf  = rhoatm * 100.0;
    double patm     = EOSkappa * pow(rhoatm,gamma);
    double patmf    = EOSkappa * pow(rhoatmf,gamma);
    double epsatm   = EOSkappa / (gamma-1) * pow(rhoatm,gamma-1.);
    double epsatmf  = EOSkappa / (gamma-1) * pow(rhoatmf,gamma-1.);
    double cc,vv,dp,SQRTdetg,ooSQRTdetg,SSSqrt,ooDpT;

    int result;
    int i;

    // loop over all points

    SQRTdetg = sqrt(*gpsi4)*(*gpsi4);
    ooSQRTdetg = 1./SQRTdetg;


    d   = ooSQRTdetg*(*gD);
    tau = ooSQRTdetg*(*gtau);
    S1  = ooSQRTdetg*(*gS);
    S2  = 0.0;
    S3  = 0.0;

    Sinv1 = S1/(*gpsi4);
    Sinv2 = S2/(*gpsi4);
    Sinv3 = S3/(*gpsi4);

    SS = fabs(S1*Sinv1 + S2*Sinv2 + S3*Sinv3);
    SSSqrt= sqrt(SS);

    result = C2P_OK;

    if (isnan(d)||isinf(d)||
         isnan(tau)||isinf(tau)||
          isnan(S1)||isinf(S1)||
          isnan(S2)||isinf(S2)||
          isnan(S3)||isinf(S3))
     {
        printf("d %e S %e tau %e psi4 %e alpha %e \n",*gD, *gS, *gtau, *gpsi4, *galpha);
        printf("d %e tau %e S1 %e S2 %e S3 %e \n",d, tau, S1, S2,S3);
        printf("D, S, or tau is nan\n");
        result = C2P_FAIL;
     }

        if (isnan(SS)||isinf(SS)||
            isnan(Sinv1)||isinf(Sinv1)||
            isnan(Sinv2)||isinf(Sinv2)||
            isnan(Sinv3)||isinf(Sinv3)) {
          printf("Sinv, or SS is nan, reset to zero, try to go on\n");
          S1 = Sinv1 = 0.;
          S2 = Sinv2 = 0.;
          S3 = Sinv3 = 0.;
          SS = 0.;
          result = C2P_FAIL;
            }

    if ((result!=C2P_FAIL)) {
         result = try_hot_part(d,SS,tau, &rho,&eps,&p, *gp, errgoal,itermax, patm, rhoatm, 100, gamma);
         //printf("HOT ");
    }

    if (result==C2P_COLD) {
          result = try_cold_part(d,SS,&tau, &rho,&eps,&p, rho, errgoal,itermax, patm, rhoatm, 100, gamma, EOSkappa);
          *gtau = SQRTdetg*tau;
          //printf("COLD ");
            }



        if ((result!=C2P_OK) && (result!=C2P_ATM)) {
//          printf("Problem in cons2prim: (%d)\n",result);
          result = C2P_ATM;
//          printf(" set ATM, try go on\n");
        }


            // set atm
            if (result==C2P_ATM) {
                p     = patm;
                rho   = rhoatm;
                eps   = epsatm;

                d     = rhoatm;
                S1    = Sinv1 = 0.;
                S2    = Sinv2 = 0.;
                S3    = Sinv3 = 0.;
                SS    = 0.;
                tau   = rhoatm*epsatm;

                *gD   = SQRTdetg*d;
                *gS   = SQRTdetg*S1;
                *gtau = SQRTdetg*tau;
            }


         // set primitives
       ooDpT  = 1.0/(d+p+tau);
       *grho  = rho;
       *geps  = eps;
       *gp    = p;
       *gv    = Sinv1 * ooDpT;

       //printf("%.15f %.15f %.15f\n",d,p,tau);
       //printf("%d %.15f %.15f %.15f %.15f %.15f %.15f\n",g->ngrid,grho[i],geps[i],gp[i],gv1[i],gv2[i],gv3[i]);
       if(isnan(*grho)) exit(0);
       if(isnan(*geps)) exit(0);
       if(isnan(*gp)) exit(0);
       if(isnan(*gv)) exit(0);

}

void grhdc_update_primitives_double_all(double **vec, PARA *par) {

    double D,S,tau,alpha,psi4,p,rho,v,eps;

    int i,j;
    for(i=0;i<par->N;i++) {
        D     = vec[0][i];
        S     = vec[1][i];
        tau   = vec[2][i];
        alpha = vec[3][i];
        psi4  = vec[4][i];
        p     = vec[8][i];

        grhdc_c2p(&D,&S,&tau,&psi4,&alpha,&rho,&v,&eps,&p,par->e_gamma);

        vec[0][i] = D;
        vec[1][i] = S;
        vec[2][i] = tau;
        vec[5][i] = rho;
        vec[6][i] = v;
        vec[7][i] = eps;
        vec[8][i] = p;
    }

}

void grhdc_update_primitives_i(gsl_vector ***vec, int length, PARA *par, int i) {

    double D,S,tau,alpha,psi4,p,rho,v,eps;

    int j;
    for(j=0;j<length;j++) {
        D     = gsl_vector_get (vec[0][i], j);
        S     = gsl_vector_get (vec[1][i], j);
        tau   = gsl_vector_get (vec[2][i], j);
        alpha = gsl_vector_get (vec[3][i], j);
        psi4  = gsl_vector_get (vec[4][i], j);
        p     = gsl_vector_get (vec[8][i], j);

        grhdc_c2p(&D,&S,&tau,&psi4,&alpha,&rho,&v,&eps,&p,par->e_gamma);

        gsl_vector_set (vec[0][i], j, D);
        gsl_vector_set (vec[1][i], j, S);
        gsl_vector_set (vec[2][i], j, tau);
        gsl_vector_set (vec[5][i], j, rho);
        gsl_vector_set (vec[6][i], j, v);
        gsl_vector_set (vec[7][i], j, eps);
        gsl_vector_set (vec[8][i], j, p);
    }

}

void grhdc_update_primitives_i_t(gsl_vector ****vec, int length, PARA *par, int i, int t) {

    double D,S,tau,alpha,psi4,p,rho,v,eps;

    int j;
    for(j=0;j<length;j++) {
        D     = gsl_vector_get (vec[0][i][t], j);
        S     = gsl_vector_get (vec[1][i][t], j);
        tau   = gsl_vector_get (vec[2][i][t], j);
        alpha = gsl_vector_get (vec[3][i][t], j);
        psi4  = gsl_vector_get (vec[4][i][t], j);
        p     = gsl_vector_get (vec[8][i][t], j);

        grhdc_c2p(&D,&S,&tau,&psi4,&alpha,&rho,&v,&eps,&p,par->e_gamma);

        gsl_vector_set (vec[0][i][t], j, D);
        gsl_vector_set (vec[1][i][t], j, S);
        gsl_vector_set (vec[2][i][t], j, tau);
        gsl_vector_set (vec[5][i][t], j, rho);
        gsl_vector_set (vec[6][i][t], j, v);
        gsl_vector_set (vec[7][i][t], j, eps);
        gsl_vector_set (vec[8][i][t], j, p);
    }

}

void grhdc_update_primitives_all(gsl_vector ***vec, int length, PARA *par) {

    int i;
    for(i=0;i<par->N;i++)
      grhdc_update_primitives_i(vec, length, par, i);

}

// ### ### ### ### ### ### ###

void grhdc_Output_Averages ( double **V, PARA *par, char *name ) {

    int i;
    char filename[255];

    double rho, v, p, eps;

    sprintf ( filename, "../data/%s/%s/%s/%d_comp_rho.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp1 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_v.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp2 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_p.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp3 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_D.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp4 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_S.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp5 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_tau.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp6 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_alpha.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp7 = fopen ( filename, "a" );
    sprintf ( filename, "../data/%s/%s/%s/%d_comp_psi4.txt", par->project, par->shape, par->scheme, par->N );
    FILE *fp8 = fopen ( filename, "a" );

    fprintf ( fp1, "\"Time=%lf\n", par->time );
    fprintf ( fp2, "\"Time=%lf\n", par->time );
    fprintf ( fp3, "\"Time=%lf\n", par->time );
    fprintf ( fp4, "\"Time=%lf\n", par->time );
    fprintf ( fp5, "\"Time=%lf\n", par->time );
    fprintf ( fp6, "\"Time=%lf\n", par->time );
    fprintf ( fp7, "\"Time=%lf\n", par->time );
    fprintf ( fp8, "\"Time=%lf\n", par->time );

    for ( i = 0; i < par->N; ++i ) {
        p=0;
//        fprintf ( fp1, "%.15f %.15f \n", -(par->xMin + (i+0.5)*par->dx),-V[5][i]+2.0);
        fprintf ( fp1, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[5][i]);
        fprintf ( fp2, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[6][i]);
        fprintf ( fp3, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[8][i]);
        fprintf ( fp4, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[0][i]);
        fprintf ( fp5, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[1][i]);
        fprintf ( fp6, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[2][i]);
        fprintf ( fp7, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[3][i]);
        fprintf ( fp8, "%.15f %.15f \n", par->xMin + (i+0.5)*par->dx,V[4][i]);
    }

    fclose ( fp1 ); fclose ( fp2 ); fclose ( fp3 ); fclose ( fp4 ); fclose ( fp5 ); fclose ( fp6 ); fclose ( fp7 ); fclose ( fp8 );

}

// ### ### ### ### ### ### ###

void grhdc_Output_Averages_single_time ( double **V, PARA *par, char *name ) {

    int i;
    char filename[255];

    double rho, v, p, eps;

    sprintf ( filename, "../data/%s/%s/%s/%d_%s_t%.5f.csv", par->project, par->shape, par->scheme, par->N, name, par->time );

    FILE *fp = fopen ( filename, "w" );

    for ( i = 0; i < par->N; ++i ) {
        fprintf ( fp, "%.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f %.15f\n", par->xMin + (i+0.5)*par->dx,
                  V[0][i],V[1][i],V[2][i],V[3][i],V[4][i],V[5][i],V[6][i],V[7][i],V[8][i] );
    }

    fclose ( fp );

}

// #### #### ####

void grhdc_Flux_Coeff_i ( ADER_DATA *A, PARA *par, int i ) {

    double rho, v, eps, p, D, S, tau, alpha, psi4, sqrtdetgamma, sqrtmdetg;

    int j;
    for ( j = 0; j < par->size; ++j ) {

      D     = gsl_vector_get (A->Vold[0][i], j);
      S     = gsl_vector_get (A->Vold[1][i], j);
      tau   = gsl_vector_get (A->Vold[2][i], j);
      alpha = gsl_vector_get (A->Vold[3][i], j);
      psi4  = gsl_vector_get (A->Vold[4][i], j);

      rho   = gsl_vector_get (A->Vold[5][i], j);
      v     = gsl_vector_get (A->Vold[6][i], j);
      eps   = gsl_vector_get (A->Vold[7][i], j);
      p     = gsl_vector_get (A->Vold[8][i], j);

      sqrtdetgamma = pow(psi4,1.5);
      sqrtmdetg    = sqrtdetgamma * alpha;

      gsl_vector_set (A->F[0][i], j, v * alpha * D);
      gsl_vector_set (A->F[1][i], j, sqrtmdetg*p + v*alpha*S);
      gsl_vector_set (A->F[2][i], j, v*(sqrtmdetg*p + alpha*tau));

    }

    gsl_vector_set_zero( A->F[3][i] ); // ALPHA     // DONT CHANGE GEOMETRY HERE
    gsl_vector_set_zero( A->F[4][i] ); // PSI4
    gsl_vector_set_zero( A->F[5][i] ); // RHO       // DONT CHANGE PRIMITIVES
    gsl_vector_set_zero( A->F[6][i] ); // V
    gsl_vector_set_zero( A->F[7][i] ); // EPS
    gsl_vector_set_zero( A->F[8][i] ); // P

    par->evalTempF += par->size;
}

void grhdc_Flux_Coeff_all ( ADER_DATA *A, PARA *par ) {

    double v, p, D, S, tau, alpha, psi4, sqrtdetgamma, sqrtmdetg;

    int i,j;
    for ( i = 0; i < par->N; ++i ) {

       grhdc_Flux_Coeff_i ( A, par, i );

    }

}

// #### #### ####

void grhdc_Source_Coeff_i_3D ( ADER_DATA *A, PARA *par, int i ) {

    double x;
    double rho, dv, v, p, eps, D, S, tau, alpha, psi4, dalpha, dpsi4, sqrtdetgamma, sqrtmdetg, W, W2hrho;
    double T00, T0i, Tii;

    gsl_vector *vdalpha = gsl_vector_calloc(par->size);
    gsl_vector *vdpsi4  = gsl_vector_calloc(par->size);

    int k,l,m;
    for(m=0; m<par->Mp1; m++)
    for(k=0; k<par->Mp1; k++) {
        gsl_vector_set(vdalpha, k + m*par->Mp1, 0.0);
        gsl_vector_set(vdpsi4,  k + m*par->Mp1, 0.0);
        for(l=0; l<par->Mp1; l++) {
            gsl_vector_set(vdalpha, k + m*par->Mp1, gsl_vector_get(vdalpha, k + m*par->Mp1) + gsl_matrix_get(A->D, l,k) * gsl_vector_get(A->Vold[3][i], l + m*par->Mp1));
            gsl_vector_set(vdpsi4 , k + m*par->Mp1, gsl_vector_get(vdpsi4 , k + m*par->Mp1) + gsl_matrix_get(A->D, l,k) * gsl_vector_get(A->Vold[4][i], l + m*par->Mp1));
        }
    }

    int j;
    for ( j = 0; j < par->size; ++j ) {

        x = par->xMin + (i + A->LAM[j%(par->Mp1)])*par->dx;

        D     = gsl_vector_get (A->Vold[0][i], j);
        S     = gsl_vector_get (A->Vold[1][i], j);
        tau   = gsl_vector_get (A->Vold[2][i], j);
        alpha = gsl_vector_get (A->Vold[3][i], j);
        psi4  = gsl_vector_get (A->Vold[4][i], j);

        rho   = gsl_vector_get (A->Vold[5][i], j);
        v     = gsl_vector_get (A->Vold[6][i], j);
        eps   = gsl_vector_get (A->Vold[7][i], j);
        p     = gsl_vector_get (A->Vold[8][i], j);

        dalpha= gsl_vector_get (vdalpha, j) / par->dx;
        dpsi4 = gsl_vector_get (vdpsi4 , j) / par->dx;

        sqrtdetgamma = pow(psi4,1.5);
        sqrtmdetg    = sqrtdetgamma * alpha;
        W            = pow((1 - (v*v*psi4)),-0.5);
        W2hrho       = W*W*(p + rho + eps*rho);

        T00 = pow(alpha,-2.)*(W2hrho + (-p));
        T0i = v*W2hrho*1/(alpha);
        Tii = (v*v*W2hrho + p*3./(psi4));

        gsl_vector_set (A->S[0][i], j, -(2.0/x)*alpha* D *v );
        gsl_vector_set (A->S[1][i], j, -(2.0/x)*alpha* S *v + sqrtmdetg*((-T00*alpha*dalpha) + 0.5*(Tii)*dpsi4) );
        gsl_vector_set (A->S[2][i], j, -(2.0/x)*alpha*tau*v - sqrtmdetg*(T0i*dalpha + (2.0/x)*p*v));

    }

    gsl_vector_set_zero( A->S[3][i] ); // ALPHA     // DONT CHANGE GEOMETRY HERE
    gsl_vector_set_zero( A->S[4][i] ); // PSI4
    gsl_vector_set_zero( A->S[5][i] ); // RHO       // DONT CHANGE PRIMITIVES
    gsl_vector_set_zero( A->S[6][i] ); // V
    gsl_vector_set_zero( A->S[7][i] ); // EPS
    gsl_vector_set_zero( A->S[8][i] ); // P

    par->evalTempS += par->size;

    gsl_vector_free(vdalpha);
    gsl_vector_free(vdpsi4);

}

void grhdc_Source_Coeff_i_1D ( ADER_DATA *A, PARA *par, int i ) {

    double x;
    double rho, dv, v, p, eps, D, S, tau, alpha, psi4, dalpha, dpsi4, sqrtdetgamma, sqrtmdetg, W, W2hrho;
    double T00, T0i, Tii;

    gsl_vector *vdalpha = gsl_vector_calloc(par->size);
    gsl_vector *vdpsi4  = gsl_vector_calloc(par->size);

    int k,l,m;
    for(m=0; m<par->Mp1; m++)
    for(k=0; k<par->Mp1; k++) {
        gsl_vector_set(vdalpha, k + m*par->Mp1, 0.0);
        gsl_vector_set(vdpsi4,  k + m*par->Mp1, 0.0);
        for(l=0; l<par->Mp1; l++) {
            gsl_vector_set(vdalpha, k + m*par->Mp1, gsl_vector_get(vdalpha, k + m*par->Mp1) + gsl_matrix_get(A->D, l,k) * gsl_vector_get(A->Vold[3][i], l + m*par->Mp1));
            gsl_vector_set(vdpsi4 , k + m*par->Mp1, gsl_vector_get(vdpsi4 , k + m*par->Mp1) + gsl_matrix_get(A->D, l,k) * gsl_vector_get(A->Vold[4][i], l + m*par->Mp1));
        }
    }

    int j;
    for ( j = 0; j < par->size; ++j ) {

        x = par->xMin + (i + A->LAM[j%(par->Mp1)])*par->dx;

        D     = gsl_vector_get (A->Vold[0][i], j);
        S     = gsl_vector_get (A->Vold[1][i], j);
        tau   = gsl_vector_get (A->Vold[2][i], j);
        alpha = gsl_vector_get (A->Vold[3][i], j);
        psi4  = gsl_vector_get (A->Vold[4][i], j);

        rho   = gsl_vector_get (A->Vold[5][i], j);
        v     = gsl_vector_get (A->Vold[6][i], j);
        eps   = gsl_vector_get (A->Vold[7][i], j);
        p     = gsl_vector_get (A->Vold[8][i], j);

        dalpha= gsl_vector_get (vdalpha, j) / par->dx;
        dpsi4 = gsl_vector_get (vdpsi4 , j) / par->dx;

        sqrtdetgamma = pow(psi4,1.5);
        sqrtmdetg    = sqrtdetgamma * alpha;
        W            = pow((1 - (v*v*psi4)),-0.5);
        W2hrho       = W*W*(p + rho + eps*rho);

        T00 = pow(alpha,-2.)*(W2hrho + (-p));
        T0i = v*W2hrho*1/(alpha);
        Tii = (v*v*W2hrho + p*3./(psi4));

        gsl_vector_set (A->S[0][i], j, 0.0 );
        gsl_vector_set (A->S[1][i], j, sqrtmdetg*((-T00*alpha*dalpha) + 0.5*(Tii)*dpsi4) );
        gsl_vector_set (A->S[2][i], j, -sqrtmdetg*T0i*dalpha);

    }

    gsl_vector_set_zero( A->S[3][i] ); // ALPHA     // DONT CHANGE GEOMETRY HERE
    gsl_vector_set_zero( A->S[4][i] ); // PSI4
    gsl_vector_set_zero( A->S[5][i] ); // RHO       // DONT CHANGE PRIMITIVES
    gsl_vector_set_zero( A->S[6][i] ); // V
    gsl_vector_set_zero( A->S[7][i] ); // EPS
    gsl_vector_set_zero( A->S[8][i] ); // P

    par->evalTempS += par->size;

    gsl_vector_free(vdalpha);
    gsl_vector_free(vdpsi4);

}

void grhdc_Source_Coeff_all ( ADER_DATA *A, PARA *par ) {
    /// not used!
//    int i;
//    for ( i = 0; i < par->N; ++i ) {
//        grhdc_Source_Coeff_i ( A, par, i );
//    }
}

// #### #### ####

double grhdc_Calc_Speed ( double **V, PARA *par, int i ) {

    double rho, v, p, eps, D, S, tau, alpha, psi4;
    double gam = par->e_gamma;
    double cc,cs,vsq;
    double lp,lm,l0,lmax;
    double GAMxx;

    D     = V[0][i];
    S     = V[1][i];
    tau   = V[2][i];
    alpha = V[3][i];
    psi4  = V[4][i];
    rho   = V[5][i];
    v     = V[6][i];
    eps   = V[7][i];
    p     = V[8][i];

    GAMxx = 1.0 / psi4;

//    cc = ((-gam) + gam*gam)*eps*1.0/((1.0 + gam*eps));  /// MARCUS - cc
    cc = ((gam*gam - gam)*eps)/(1.0 + gam*eps);   /// MATTHIAS - cc
    cs = sqrt(cc);
    vsq = v * v *psi4;

//    lp = fabs(1.0/((1.0 + (-cc*vsq)))*(v*(alpha + (-cc*alpha)) + (cs*alpha*sqrt((v*v*(-1.0 + cc + vsq + (-cc*vsq)) + (-vsq*(1.0/(psi4) + cc*1.0/(psi4))) + 1.0/(psi4) + cc*vsq*vsq*1.0/(psi4))))));    /// MARCUS - lp

//    lm = fabs(1.0/((1.0 + (-cc*vsq)))*(v*(alpha + (-cc*alpha)) + (-cs*alpha*sqrt((v*v*(-1.0 + cc + vsq + (-cc*vsq)) + (-vsq*(1.0/(psi4) + cc*1.0/(psi4))) + 1.0/(psi4) + cc*vsq*vsq*1.0/(psi4))))));    /// MARCUS - lm


    lp = fabs ( alpha * ( v * ( 1.0 - cc )
                    + ( cs*sqrt( (1.0-vsq)*( GAMxx*( 1.0 - vsq*cc ) - v*v*( 1.0 - cc ) ) ) ) )
                    / ( 1.0 - vsq*cc ) );    /// MATTHIAS - lp

    lm = fabs ( alpha * ( v * ( 1.0 - cc )
                    - ( cs*sqrt( (1.0-vsq)*( GAMxx*( 1.0 - vsq*cc ) - v*v*( 1.0 - cc ) ) ) ) )
                    / ( 1.0 - vsq*cc ) );    /// MATTHIAS - lm


    l0 = fabs(v*alpha);

    lmax = fmax(fmax(l0,lm), lp);

    return lmax;

}

// ### #### ###

void grhdc_Convergence_Output ( FV_DATA *V, PARA *par, char *name ) {

    ASSERT ( V );
    ASSERT ( par );
    ASSERT ( name );

    int i;
    char file1[255];
    char file2[255];

    sprintf ( file1, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_%s.txt", par->project, par->M_sub, par->shape, par->scheme, par->CFL, par->N, name );
    FILE *fp1 = fopen ( file1, "w" );

    sprintf ( file2, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/v_%s.txt", par->project, par->M_sub, par->shape, par->scheme, par->CFL, par->N, name );
    FILE *fp2 = fopen ( file2, "w" );

    /// #### ATTENTION --- changed to CSV format --- not compatible with convergence routine! TODO

    for ( i = 0; i < par->N-1; ++i ) {
        fprintf ( fp1, "%.15f, %.15f\n", par->xMin+(i+0.5)*par->dx, V->Avg[5][i] );
        fprintf ( fp2, "%.15f, %.15f\n", par->xMin+(i+0.5)*par->dx, V->Avg[6][i] );
    }

    fprintf ( fp1, "%.15f, %.15f\n", par->xMin+(par->N-0.5)*par->dx, V->Avg[5][(par->N-1)] );
    fprintf ( fp2, "%.15f, %.15f\n", par->xMin+(par->N-0.5)*par->dx, V->Avg[6][(par->N-1)] );

    fclose ( fp1 );
    fclose ( fp2 );

    printf ( "# SUCCESS # - output convergence file: ../%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/rho_%s.txt\n", par->project, par->M_sub, par->shape, par->scheme, par->CFL, par->N, name);
    printf ( "# SUCCESS # - output convergence file: ../%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d/v_%s.txt\n\n", par->project, par->M_sub, par->shape, par->scheme, par->CFL, par->N, name);

}

// #### #### ####

void grhdc_F_RK ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t ) {

    gsl_vector_memcpy   ( F[0], q[0][i][t] );
    gsl_vector_mul      ( F[0], F[0] );
    gsl_vector_scale    ( F[0], 0.5 );

    double rho, v, p, eps, D, S, tau, alpha, psi4, sqrtdetgamma, sqrtmdetg;

    int j;
    for ( j = 0; j < par->Mp1; ++j ) {

      D     = gsl_vector_get (q[0][i][t], j);
      S     = gsl_vector_get (q[1][i][t], j);
      tau   = gsl_vector_get (q[2][i][t], j);
      alpha = gsl_vector_get (q[3][i][t], j);
      psi4  = gsl_vector_get (q[4][i][t], j);

      p=0;
      grhdc_c2p(&D, &S, &tau, &psi4, &alpha, &rho, &v, &eps, &p, par->e_gamma);

      sqrtdetgamma = pow(psi4,1.5);
      sqrtmdetg    = sqrtdetgamma * alpha;

      gsl_vector_set (F[0], j, v * alpha * D);
      gsl_vector_set (F[1], j, sqrtmdetg*p + v*alpha*S);
      gsl_vector_set (F[2], j, v*(sqrtmdetg*p + alpha*tau));

    }

    gsl_vector_set_zero( F[3] ); // ALPHA     // DONT CHANGE GEOMETRY HERE
    gsl_vector_set_zero( F[4] ); // PSI4


    par->evalTempF += par->Mp1;

}

void grhdc_S_RK_1D ( gsl_vector **Src, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *Der, double *LAM ) {

    double x;
    double rho, v, p, eps, D, S, tau, alpha, psi4, dalpha, dpsi4, sqrtdetgamma, sqrtmdetg, W, W2hrho;
    double T00, T0i, Tii;

    gsl_vector *vdalpha = gsl_vector_calloc(par->Mp1);
    gsl_vector *vdpsi4  = gsl_vector_calloc(par->Mp1);

    int k,l;
    for(k=0; k<par->Mp1; k++) {
        gsl_vector_set(vdalpha, k, 0.0);
        gsl_vector_set(vdpsi4,  k, 0.0);
        for(l=0; l<par->Mp1; l++) {
            gsl_vector_set(vdalpha, k, gsl_vector_get(vdalpha, k) + gsl_matrix_get(Der, l,k) * gsl_vector_get(q[3][i][t], l));
            gsl_vector_set(vdpsi4,  k, gsl_vector_get(vdpsi4,  k) + gsl_matrix_get(Der, l,k) * gsl_vector_get(q[4][i][t], l));
        }
    }

    int j;
    for ( j = 0; j < par->Mp1; ++j ) {

        x = par->xMin + (i + LAM[j])*par->dx;

        D     = gsl_vector_get (q[0][i][t], j);
        S     = gsl_vector_get (q[1][i][t], j);
        tau   = gsl_vector_get (q[2][i][t], j);
        alpha = gsl_vector_get (q[3][i][t], j);
        psi4  = gsl_vector_get (q[4][i][t], j);

        p=0;
        grhdc_c2p(&D, &S, &tau, &psi4, &alpha, &rho, &v, &eps, &p, par->e_gamma);

        dalpha= gsl_vector_get (vdalpha, j) / par->dx;
        dpsi4 = gsl_vector_get (vdpsi4 , j) / par->dx;

        sqrtdetgamma = pow(psi4,1.5);
        sqrtmdetg    = sqrtdetgamma * alpha;
        W            = 1.0 / ( sqrt(1 - (v*v*psi4)) );
        W2hrho       = W*W*(p + rho + eps*rho);

        T00 = (W2hrho - p)/(alpha*alpha);
        T0i = (v*W2hrho)/(alpha);
        Tii = (v*v*W2hrho + p*3./(psi4));

        gsl_vector_set (Src[0], j, 0.0 );
        gsl_vector_set (Src[1], j, sqrtmdetg*((-T00*alpha*dalpha) + 0.5*Tii*dpsi4) );
        gsl_vector_set (Src[2], j, -sqrtmdetg*(T0i*dalpha) );
    }

    gsl_vector_set_zero( Src[3] ); // ALPHA     // DONT CHANGE GEOMETRY HERE
    gsl_vector_set_zero( Src[4] ); // PSI4

    par->evalTempS += par->Mp1;

    gsl_vector_free(vdalpha);
    gsl_vector_free(vdpsi4);

}

void grhdc_S_RK_3D ( gsl_vector **Src, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *Der, double *LAM ) {

    double x;
    double rho, v, p, eps, D, S, tau, alpha, psi4, dalpha, dpsi4, sqrtdetgamma, sqrtmdetg, W, W2hrho;
    double T00, T0i, Tii;

    gsl_vector *vdalpha = gsl_vector_calloc(par->Mp1);
    gsl_vector *vdpsi4  = gsl_vector_calloc(par->Mp1);

    int k,l;
    for(k=0; k<par->Mp1; k++) {
        gsl_vector_set(vdalpha, k, 0.0);
        gsl_vector_set(vdpsi4,  k, 0.0);
        for(l=0; l<par->Mp1; l++) {
            gsl_vector_set(vdalpha, k, gsl_vector_get(vdalpha, k) + gsl_matrix_get(Der, l,k) * gsl_vector_get(q[3][i][t], l));
            gsl_vector_set(vdpsi4,  k, gsl_vector_get(vdpsi4,  k) + gsl_matrix_get(Der, l,k) * gsl_vector_get(q[4][i][t], l));
        }
    }

    int j;
    for ( j = 0; j < par->Mp1; ++j ) {

        x = par->xMin + (i + LAM[j])*par->dx;

        D     = gsl_vector_get (q[0][i][t], j);
        S     = gsl_vector_get (q[1][i][t], j);
        tau   = gsl_vector_get (q[2][i][t], j);
        alpha = gsl_vector_get (q[3][i][t], j);
        psi4  = gsl_vector_get (q[4][i][t], j);

        p=0;
        grhdc_c2p(&D, &S, &tau, &psi4, &alpha, &rho, &v, &eps, &p, par->e_gamma);

        dalpha= gsl_vector_get (vdalpha, j) / par->dx;
        dpsi4 = gsl_vector_get (vdpsi4 , j) / par->dx;

        sqrtdetgamma = pow(psi4,1.5);
        sqrtmdetg    = sqrtdetgamma * alpha;
        W            = 1.0 / ( sqrt(1 - (v*v*psi4)) );
        W2hrho       = W*W*(p + rho + eps*rho);

        T00 = (W2hrho - p)/(alpha*alpha);
        T0i = (v*W2hrho)/(alpha);
        Tii = (v*v*W2hrho + p*3./(psi4));

        gsl_vector_set (Src[0], j, -2.0/(x)*alpha* D *v );
        gsl_vector_set (Src[1], j, -2.0/(x)*alpha* S *v + sqrtmdetg*((-T00*alpha*dalpha) + 0.5*(Tii)*dpsi4) );
        gsl_vector_set (Src[2], j, -2.0/(x)*alpha*tau*v - sqrtmdetg*(T0i*dalpha + 2.0/(x)*p*v));

    }

    gsl_vector_set_zero( Src[3] ); // ALPHA     // DONT CHANGE GEOMETRY HERE
    gsl_vector_set_zero( Src[4] ); // PSI4

    par->evalTempS += par->Mp1;

    gsl_vector_free(vdalpha);
    gsl_vector_free(vdpsi4);

}


// #### #### ####

void grhdc_Phy_Admiss ( DG_DATA *D, PARA *par ) {
/*
 *  for TOV star, check whether velocity exceeds 10-2 and other physics motivated criteria
 */

    double rho, p;

    if ( strcmp ( par->shape, "TOV" ) == 0 ) {
        int i, j;

        for ( i = 0; i < par->N; ++i )
            for ( j = 0; j < par->Mp1; ++j ) {
                if ( fabs(gsl_vector_get(D->pU[6][i], j)) > 0.003 ) {
                    D->beta[6][i] = NEED_AW;
//                    D->beta[6][par->N-i] = NEED_AW;   // symmetry!
                }

                rho = gsl_vector_get(D->pU[5][i], j);
                p   = gsl_vector_get(D->pU[8][i], j);

                if ( (rho < 0) || (p < 0) ) {
                    D->beta[8][i] = NEED_AW;
//                    D->beta[8][par->N-i] = NEED_AW;   // symmetry!
                }

                if ( rho > 1e-8 && rho < 3e-6 ) {
                    D->beta[5][i] = NEED_AW;
//                    D->beta[5][par->N-i] = NEED_AW;   // symmetry!
                }
            }

    }

}








































