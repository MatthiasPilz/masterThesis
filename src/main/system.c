#include "ader.h"
#include "main.h"


void Delete_Data ( PARA *par ) {
    printf ( "already existing files in ../data/ :\n" );
    if ( system ( "ls -la ../data/" ) );
    printf ( "\n" );

    int flag = 0;
    char temp[255];
    char command[511];
    int r;

    do {
        printf ( "Do you want to delete any old data? (y/n) " );
        r = scanf ( "%s", temp );

        if ( !strcmp ( temp, "y" ) ) {
            flag = 1;

            sprintf ( command, "rm -r ../data/%s", par->project );
            if ( system ( command ) );
                printf ( "../data/%s deleted\n", par->project );
        }
        else if ( !strcmp ( temp, "n" ) ) {
            flag = 1;
            printf ( "No files deleted!\n" );
        }
        else
            printf ( "y or n?\n" );

    } while ( flag != 1 );

    if ( r != 1 )
        printf ( "r = %d\n", r );

    printf ( "\n" );

}

// ### ### ### ### ### ### ###

void Copy_Parameter_Data ( char *name, char *path ) {

    char command[255];
    sprintf ( command, "cp ../par/%s %s/%s", name, path, name );

    printf ( "%s\n", command );

    if ( !system ( command ) )
        printf ( "# SUCCESS # - copied parameter data.\n" );

}

void Setup_Convergence_Data_Structure ( PARA *par ) {

    ASSERT ( par );

    printf ( "\nSetting up data structure for convergence:\n" );

    // 'project'_convergence_M'M'
    char command[255];
    sprintf ( command, "mkdir ../data/%s_convergence_M%d", par->project, par->M_sub );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d\n", par->project, par->M_sub );

    // 'project'_convergence_M'M'/'shape'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s", par->project, par->M_sub, par->shape );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s\n", par->project, par->M_sub, par->shape );

    // 'project'_convergence_M'M'/'shape'/'scheme'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s", par->project, par->M_sub, par->shape, par->scheme );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s\n", par->project, par->M_sub, par->shape, par->scheme );

    // 'project'_convergence_M'M'/'shape'/'scheme'/CFL_'CFL'
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s/CFL_%.8f", par->project, par->M_sub, par->shape, par->scheme, par->CFL );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s/CFL_%.8f\n", par->project, par->M_sub, par->shape, par->scheme, par->CFL );

    // 'project'_convergence_M'M'/'shape'/'scheme'/CFL_'CFL'/N_'N'

    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d", par->project, par->M, par->shape, par->scheme, par->CFL, par->N );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/N_%d\n", par->project, par->M, par->shape, par->scheme, par->CFL, par->N );

    // make directory for iteration statistics
    sprintf ( command, "mkdir ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic", par->project, par->M_sub, par->shape, par->scheme, par->CFL );
    if ( !system ( command ) )
        printf ( "\t created path: ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic\n", par->project, par->M_sub, par->shape, par->scheme, par->CFL );

    // write header if at N_start
    if ( par->conv_N_start == par->N ) {
        sprintf ( command, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic/grid\n", par->project, par->M, par->shape, par->scheme, par->CFL_sub );
        FILE *fp = fopen ( command, "a" );
        fprintf ( fp, "M=%d\n", par->M );
        fprintf ( fp, "N\titer\tf_pre\tf_cor\tf_sum\n" );
        fclose ( fp );

        sprintf ( command, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic/subcell\n", par->project, par->M, par->shape, par->scheme, par->CFL_sub );
        FILE *fp1 = fopen ( command, "a" );
        fprintf ( fp1, "M=%d\n", par->M );
        fprintf ( fp1, "N\titer\tf_pre\tf_cor\tf_sum\n" );
        fclose ( fp1 );
    }

    // copy parameter file to according folder
    sprintf ( command, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f", par->project, par->M_sub, par->shape, par->scheme, par->CFL );
    Copy_Parameter_Data ( par->file, command );


}

// ### ### ### ### ### ### ###

void Delete_old_Iteration_File ( PARA *par ) {

    char command[255];
    int k;

    for ( k = par->conv_M_start; k <= par->conv_M_end; ++k ) {

        sprintf ( command, "rm -r ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic/*\n", par->project, k, par->shape, par->scheme, par->CFL_sub );
        if ( !system ( command ) )
            printf ( "\t deleted all files in: ../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic/\n", par->project, k, par->shape, par->scheme, par->CFL_sub );

    }

}

// #### #### ####

void Save_Iteration_Statistic ( PARA *par, PARA *subcell ) {

    char filename[255];

    // output for grid
    sprintf ( filename, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic/grid\n", par->project, par->M, par->shape, par->scheme, par->CFL_sub );

    FILE *fp = fopen ( filename, "a" );
    fprintf ( fp, "%d\t%.5f\t%.5f\t%.5f\t%.5f\n", par->N, (double)par->iterAvg/par->count, (double)par->evalPrediF/(par->N*par->count), (double)(par->evalCountF-par->evalPrediF)/(par->N*par->count), (double)(par->evalCountF)/(par->N*par->count) );
    fclose ( fp );

    // output for subcells
    sprintf ( filename, "../data/%s_convergence_M%d/%s/%s/CFL_%.8f/Iteration_Statistic/subcell\n", par->project, par->M, par->shape, par->scheme, par->CFL_sub );

    FILE *fp2 = fopen ( filename, "a" );
    fprintf ( fp2, "%d\t%.5f\t%.5f\t%.5f\t%.5f\n", par->N, subcell->iterAvg/(subcell->N_trouble + subcell->N_neigh), (double)subcell->evalPrediF/(subcell->N*(subcell->N_trouble + subcell->N_neigh)), (double)(subcell->evalCountF-subcell->evalPrediF)/(subcell->N*(subcell->N_trouble + subcell->N_neigh)), (double)subcell->evalCountF/(subcell->N*(subcell->N_trouble + subcell->N_neigh)) );
    fclose ( fp2 );

    printf ( "# SUCCESS # - saved iteration statistics\n" );
}
