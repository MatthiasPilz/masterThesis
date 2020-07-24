#include "ader.h"
#include "grhdc.h"

// ### ### ###

void Hack_Set_DG_Parameter ( PARA *par ) {
    // set parameters for small grid inside TOV star

    par->xMin = 0.0;
    par->N = (int)(par->hack_bound_R/par->dx) + 1;  // one ghost cell on the outside!?
    par->xRange = par->N * par->dx;

    ASSERT ( par->xRange < 8.0 );   // TOV star ~ r=8.1

}

// ### ### ###

void Hack_Save_TOV_Boundary_Flux ( FV_DATA *V, ADER_DATA *A, PARA *par ) {

}
