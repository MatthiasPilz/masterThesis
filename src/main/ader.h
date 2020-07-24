/*
 * @file   ader.h
 * @author mapi
 * @date   03/2015
 *
 * Main header file of ader program. Includes all *.h files and defines the
 * master structures.
 */

#ifndef __ADER_H__
#define __ADER_H__


/*****************************************************************************/
/* Inclusion of standard c and system header files                           */
/*****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*****************************************************************************/
/* Inclusion of necessary GSL header files     				                 */
/*****************************************************************************/
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_vector.h"


/*****************************************************************************/
/* DEBUG flag - don't forget to turn off for serious runs!                   */
/*****************************************************************************/
//#define DEBUG

#ifndef DEBUG
#define ASSERT(n)
#else
#define ASSERT(n) \
if(!(n)) { \
printf("%s - Failed\n",#n); \
printf("  On %s\n",__DATE__); \
printf("  At %s\n",__TIME__); \
printf("  In File %s\n",__FILE__); \
printf("  At Line %d\n",__LINE__); \
exit(1);}
#endif


/*****************************************************************************/
/* CONVERGENCE test flag                                                     */
/*****************************************************************************/
#define CONVERGENCE


/*****************************************************************************/
/* OUTPUT flags                                                              */
/*****************************************************************************/
#define NEWLINE printf ( "\n" );

//#define __OUTPUT__
#ifndef __OUTPUT__
#define OUTPUT_EVO(n)
#else
#define OUTPUT_EVO(n) printf ( "# SUCCESS # - %s\n", n );
#endif


/*****************************************************************************/
/* Different DATA structs 									                 */
/* parameters with three '///' are given/precomputed                         */
/* parameters with two '//' will be calculated during the run                */
/*****************************************************************************/
typedef struct {

    double **Avg;           ///     cell averages
    double **Avg_adm;       // cell averages for admissibility tests
    double **Flux;          // fluxes
    double **Source;        // source
    double **hack_Flux;     // flux on boundary inside the star (only one value per variable)

    gsl_vector *F1_0;       ///     linear flux coefficients at 'xi'=0
    gsl_vector *F1_1;       ///     linear flux coefficients at 'xi'=1
    gsl_vector *S1;         ///     linear source coefficients

} FV_DATA;


typedef struct {

    gsl_vector ****gslData;     // stencil data...
    gsl_vector ****gslWs;       // stencil coefficients...
    gsl_matrix **Rec;			/// reconstruction matrices

    double *l_s;                /// lambda_s --> precomputable
    double **OIM;               /// oscillation indicator matrix --> precomputable

    double **w_s;               // temporary variable w_s as in eq.(17)
    double **w_tilde;           // temporary variable w_tilde as in eq.(17)
    double **sigma;             // oszillation indicator sigma as in eq.(18)

    int *shift;                 /// stencil shift for different orders

    gsl_vector ***WT;           // WENO-coefficients

} WENO_DATA;


typedef struct {

    gsl_vector ***Vold;     // data for iteration scheme
    gsl_vector ***Vnew;     // data for iteration scheme

    gsl_vector ***F;        // flux coefficients --> different for each project!
    gsl_vector ***S;        // source coefficients
    gsl_vector ***hack_F;   // flux coefficients on boundary inside the star
	
    gsl_matrix *I;          ///     precomputable
    gsl_matrix *B;          ///     precomputable
    gsl_matrix *Ms;         ///     precomputable

    double *LAM;            ///     Gauss-Legendre Nodes for output

// RK stuff
    gsl_matrix *D;          ///     precomputable derivative matrix (lagrange)
    double **Aij;           ///     Butcher Aij Matrix (7x7)
    gsl_vector ****k_RK;     // k values in RK scheme --> only one variable!
    gsl_vector ****q_RK;     // variable values in RK scheme...
    gsl_vector **F_RK;       // flux coefficients in RK scheme
    gsl_vector **S_RK;       // source coefficients in RK scheme
    double **Dense;         ///     precomputable dense output polynomials
	
} ADER_DATA;


typedef struct {

    gsl_vector ***U;        // polynomial data coefficients
    gsl_vector ***pU;       // predictor for U
    double **V;             // subcell averages
    double **V_grid;        // grid cell averages
    double **V_grid_pre;    // grid cell averages from predictor

    int **beta;             // troubled cell indicator
    int **beta_sub;         // troubled cell indicator for subcells
    int **beta_prev;        // previous troubled cells

    gsl_matrix *Int;        ///     interpolation matrix
    gsl_matrix *Rec;        ///     reconstruction matrix
    double *AVG;            ///     vector for calculation of direct average of cells
    double *LAM;            ///     Gauss-Legendre Nodes for output
    double *BLow;           ///     vector to compute lower boundary polynomial value
    double *BUp;            ///     vector to compute upper boundary polynomial value

    gsl_matrix *F;          ///     precomputable
    gsl_matrix *G00;        ///     precomputable
    gsl_matrix *G01;        ///     precomputable
    gsl_matrix *G10;        ///     precomputable
    gsl_matrix *G11;        ///     precomputable
    gsl_matrix *H;          ///     precomputable

} DG_DATA;


/*****************************************************************************/
/* Parameter struct contains ALL problem relevant variables                  */
/*****************************************************************************/
typedef struct {
	
    int M;				///     'order' of scheme
    int M_sub;          ///     'order' of ader weno subcell scheme
	int size;			// pow( M+1, 2 )
    int Mp1;            // M+1
    int n_stencil;      // number of stencils per reconstruction
    char scheme[255];   ///     choose between RK or DUMBSER

// grid setup for (so far only equidistant discretisation)
	int N;				/// 	number of cells on interval
	double xRange;		/// 	xRange of interval
	double xMin;		/// 	smallest x value
    double xMax;		// biggest x value
	double dx;			// cell size
    int NSC;            // number of subcells per cell (2M+1)
	
    double CFL;			// CFL factor
    double CFL_sub;     ///     CFL factor for subcell ADER-WENO
	
	int Nt;				// number of time steps
	double tEND;		/// 	final time
	double dt;			// time step size
	double time;		// current time
    int count;          // current time step

    char boundary[255]; ///     type of boundary condition

// WENO reconstruction variables
    int WENO_R;         ///     exponent in weight calculation  ( [1]: r = 8 )
    double WENO_eps;    ///     eps in weight calculation       ( [1]: eps = 1e-14 )
    double WENO_lam;    ///     lambda for centered stencils    ( [1]: lam = 1e5 )
	
// iteration scheme variables - for all problems identical	
    int iterMAX;		//	maximum number of iterations --> M+1
    int iterEXTRA;      ///     additional number of iterations
	double eps;			/// 	exit condition for iteration scheme
    double chanAvg;		// average value of accuracy in iteration scheme
	double changeMAX;	// worst value of accuracy in iteration scheme
    double iterAvg;     // average number of iterations
    double iterAvgTemp; // average number of iterations hack...
    int iterFail;		// flag if iteration didn't finish in M+1 steps
    int iterFailCount;  // last time step which set iterFail flag

// ader_dg variables
    int admiss_num;     // flag wheter numerical admissibility test failed;
    int admiss_phy;     // flag wheter physical admissibility test failed;
    double N_trouble;   // average number of troubled cells
    double N_neigh;     // average number of neighbours to troubled cells
    double kappa;       ///     M-factor in Qiu Shu Troubled Cell Indicator
    double MH2;         // kappa*dx*dx

	int output;			/// 	output every ..th timestep
    int Npoly;          ///     points per cell for polynomial output

// function evaluation statistics
    long evalTempF;     // helping variable for f evals
    long evalCountF;    // number of evaluations of f(u) total
    long evalPrediF;    // number of f evals for predictor
    long evalTempS;     // helping variable for S evals
    long evalCountS;    // number of evaluations of S(u) total
    long evalPrediS;    // number of S evals for predictor
    /// hack
        long evalCF;    // temp variable
        long evalPF;    // temp variable
        long evalCS;    // temp variable
        long evalPS;    // temp variable
	
// Runge-Kutta parameter
    int stages;         ///     stages of RK scheme
    int stage_flag;     // flag for const stages or according to M

// project specific variables	
    char project[255];  /// 	Problem definition
    int var;			// number of degrees of freedom
    char shape[255];    /// 	initial values
    int source_flag;    // sources?
    double speed;       ///     TODO ### ### ### ###
    double e_gamma;     ///     for euler equations
    int dimension;      ///     pseudo dimension (only for TOV)

    double hack_bound_R;///     evolve TOV star only in domain [0,hack_bound_R]

    int *sym_flags;     // Symmetry flags

// convergence test
    int conv_N_start;   ///     lowest resolution
    int conv_N_end;     ///     highest resolution
    int conv_M_start;   ///     lowest order
    int conv_M_end;     ///     highest order
    int change_N;       ///     change factor
    int change_M;       ///     change summand

    char file[255];     ///     parameter file name

} PARA;


typedef struct {
	
    void    (*f_var_num)           ( PARA *par );
    void    (*f_init_data)         ( FV_DATA *V, PARA *par );
    void    (*f_flux_coeff_i)      ( ADER_DATA *A, PARA *par, int i );
    void    (*f_flux_coeff_all)    ( ADER_DATA *A, PARA *par );
    void    (*f_source_coeff_i)    ( ADER_DATA *A, PARA *par, int i );
    void    (*f_source_coeff_all)  ( ADER_DATA *A, PARA *par );
    double  (*f_calc_speed)        ( double **V, PARA *par, int i );
    void    (*f_post_evolution)    ( double **V, PARA *par, char *name );
    void    (*f_convergence_output)( FV_DATA *V, PARA *par, char *name );
    void    (*f_phy_admiss)        ( DG_DATA *D, PARA *par );
    void    (*f_calc_F_RK)         ( gsl_vector **F, gsl_vector ****q, PARA *par, int i, int t );
    void    (*f_calc_S_RK)         ( gsl_vector **S, gsl_vector ****q, PARA *par, int i, int t, gsl_matrix *D, double *LAM);
    void    (*f_output_averages)   ( double **V, PARA *par, char *name );

} PROJ;


/*****************************************************************************/
/* Enum variables (as extern!?)                                              */
/*****************************************************************************/
enum BetaStatus { OK, AW_NEIGH, NEED_AW };


/*****************************************************************************/
/* Macros and Global Function Pointer                                        */
/*****************************************************************************/
#define SIGNOF(a) ( ( (a) < 0 )  ?  -1   : ( (a) > 0 ) )

void (*f_predictor_step)
        ( ADER_DATA *A, gsl_vector ***U, int **beta, PARA *par, PROJ *pro );


/*****************************************************************************/
/* Inclusion of all 'sub'-headers                                            */
/*****************************************************************************/
#include "../main/main.h"

#include "../ader_weno/ader_weno.h"

#include "../data/data.h"

#include "../grid/grid.h"

#include "../output/output.h"

#include "../projects/projects.h"
#include "../projects/advection/advection.h"
#include "../projects/euler/euler.h"
#include "../projects/wave/wave.h"
#include "../projects/burger/burger.h"
#include "../projects/grhdc/grhdc.h"

#endif
