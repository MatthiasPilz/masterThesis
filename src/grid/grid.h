// prototype declaration for grid subdir


// ### interpolation.c ###
void Interpolation                  ( gsl_vector ***U, double **V, gsl_matrix *Int, PARA *parU, PARA *parV );
void Interpolation_Prev_Troubled    ( DG_DATA *DG, double **V, gsl_matrix *Int, PARA *parU, PARA *parV );
void Interpolation_Good_Cells       ( DG_DATA *DG, double **V, gsl_matrix *Int, PARA *parU, PARA *parV );
void Reconstruction                 ( double **V, gsl_vector ***U, int **beta, gsl_matrix *Rec, PARA *parV, PARA *parU );
void Cell_Averages                  ( gsl_vector ***U, double **V, double *AVG, PARA *par );

// ### admissibility.c ###
void Numerical_Admissibility        ( DG_DATA *DG, double **V, PARA *par );
void Troubled_Cell_Indicator        ( DG_DATA *DG, PARA *par );
