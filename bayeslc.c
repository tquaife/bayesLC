/*
Given n land cover proportions, and a corresponding 
confusion matrix of size nxn calculate the posterior
mean and std deviations of the joint probability
distribution.

Based on Ed Cripps' MATLAB code and the statistical 
model described in:

Cripps et al. (2007) Modelling uncertainty in
satellite-derived land cover maps. Bayesian Analysis.

Requires my matrix i/o and strparse libraries and the GSL.

Tristan Quaife. 2007.
tquaife@geog.ucl.ac.uk
*/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

#include<matrixio.h>
#include<strparse.h>

#include<gsl_randist.h>
#include<gsl_rng.h>

/*
#define MAX_LC 50
#define MAX_FLX 7
*/
#define MAX_LC 15
#define MAX_FLX 25

#define MAX_LINE_LEN 1000

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif



/*
===================
Function prototypes
===================
*/

int conditional_ar( double*, double*, long*, int, int, int );
int count_sites( long *, int, int );
int extract_masked(  double *, double *, long *, int, int );
int extract_masked_long( long *, long *, long *, int, int );
void output_to_ascii_matrix( double *, long *, int, int, FILE * ) ;
void cl_parser( int, char **, double *, int *, int * );
void usage( char * );

/*
============================
Some basic utility functions
============================
*/


void check_dims( long nx, long ny, long cx, long cy )
/*Check the dimensions of two matrices/images.
Exit with failure if they are not equal.*/
{
	if( (nx != cx) || (ny != cy) ){
		fprintf( stderr, "%s: matrix dimensions do not match\n", __FILE__ );
		exit( EXIT_FAILURE );
	}

	return ;
}




long *load_ldMatrix( long *xsize, long *ysize, char *filename )
/*Wrapper function for fread_ascii_ldMatrix that handles
the file opening and exits with failure if open fails*/
{

	FILE * fp ;
	long *ptr ;
	
	if( ( fp = fopen( filename, "r" ) ) == NULL ){
		fprintf( stderr, "%s: unable to open file %s\n", __FILE__, filename );
		exit( EXIT_FAILURE );
	}
	
	ptr = fread_ascii_ldMatrix( xsize, ysize, fp ) ;
	
	fclose( fp );
	
	return( ptr );

}

double *load_ffMatrix( long *xsize, long *ysize, char *filename )
/*Wrapper function for fread_ascii_ffMatrix that handles
the file opening and exits with failure if open fails*/
{

	FILE * fp ;
	double *ptr ;
	
	if( ( fp = fopen( filename, "r" ) ) == NULL ){
		fprintf( stderr, "%s: unable to open file %s\n", __FILE__, filename );
		exit( EXIT_FAILURE );
	}
	
	ptr = fread_ascii_ffMatrix( xsize, ysize, fp ) ;
	
	fclose( fp );
	
	return( ptr );

}



/*
===============================
========= main ================
===============================
*/

int main( int argc, char **argv )
{


	double *lc[ MAX_LC ], *lc_vec[ MAX_LC ];
	double *priors[ MAX_LC ], *priors_vec[ MAX_LC ];
	double *joint_prob, *post_prob ;
	double *cfm, *fwd_cond_prob ;
	double *work1, *work2 ;
	double sum_jprob ;
	double *sum_post[MAX_LC+MAX_FLX], *sum_post_sqr[MAX_LC+MAX_FLX] ;
	double dcorr = 1.0 ;
	long nx, ny, cx, cy;
	long *msk, *counts, *counts_vec, *tmp_vec ;
	int length=5, n_sim=300, n_props=0, i, j, k, m, n;
	int nline = 0, n_sites = 0, n_fluxes = 0 ;
	
	double *flux[MAX_FLX][MAX_LC], *flux_vec[MAX_FLX][MAX_LC], flx ;
	double *out_flux_means_vec, *out_flux_stdvs_vec ;
	
	double spatial_integral[MAX_LC+MAX_FLX] ;	
	
	double spatial_integral_pftflx[MAX_LC*MAX_FLX] ;
	
	double *cell_area, *cell_area_vec ;
	double *land_area, *land_area_vec ;
	int weight_intgrl_by_cell_area = FALSE ;
	int weight_intgrl_by_land_area = FALSE ;
	
	unsigned int *work_int ;
	
	char line[ MAX_LINE_LEN ] ;	
	char tmpc[ MAX_LINE_LEN ] ;	
	
	FILE *ifp, *out_means[MAX_LC+MAX_FLX], *out_stdvs[MAX_LC+MAX_FLX] ;
	
	gsl_rng *rnd_num_gen ;
		
	
	/* 
	---------------------------
	Setup the GSL random number 
	generator, using 0 as seed. 
	---------------------------
	*/
	gsl_rng_env_setup(  );
	rnd_num_gen = gsl_rng_alloc( gsl_rng_default ) ;
	gsl_rng_set( rnd_num_gen, 0 );

	
	/*
	-----------------------------
	Get options from command line
	-----------------------------
	*/

	cl_parser( argc, argv, &dcorr, &length, &n_sim );


	
	/*
	------------------------------------
	Open input file and load input data.
	Anything following a '#' is treated
	as a comment and blank lines are 
	skipped.
	------------------------------------
	*/
	


	if( argc < 2 ){
		fprintf( stderr, "%s: requires at least one file name as an argument\n", *argv );
		exit( EXIT_FAILURE );
	}
	if( ( ifp = fopen( argv[argc-1], "r" ) ) == NULL ){
		fprintf( stderr, "%s: unable to open file %s\n", *argv, argv[argc-1] );
		exit( EXIT_FAILURE );
	}


	while( fgets( line, MAX_LINE_LEN, ifp ) != NULL ){
					
		strip_comments( line, '#' );	
		if( is_string_blank( line ) ) continue ;
		
		nline++ ;
		
		/*
		First line - 
		  names of files containing land cover proportions
		*/
		
		if( nline == 1 ){
			while( get_first_string_element( line, tmpc ) ){
				lc[ n_props ] = load_ffMatrix( &nx, &ny, tmpc );
				if( n_props > 0 ) check_dims( nx, ny, cx, cy );
				n_props++ ;
				cx = nx ;
				cy = ny ;
			}
			
		/*
		Second line - 
		  single file name: mask file
		*/
			
		}else if( nline == 2 ){
			get_first_string_element( line, tmpc ) ;
			msk = load_ldMatrix( &cx, &cy, tmpc );			
			check_dims( nx, ny, cx, cy );	
		
		/*
		Third line - 
		  single file name: confusion matrix file
		*/
		
		}else if( nline == 3 ){
			get_first_string_element( line, tmpc ) ;
			cfm = load_ffMatrix( &cx, &cy, tmpc );			
			check_dims( n_props, n_props, cx, cy );
		
		/*
		Fourth line - 
		  single file name: counts file
		
		if a second file name exists, this should
		contain the area (in arbitrary units) of
		each grid cell.
		
		if a third filename exists this should 
		contain the PROPORTION of the grid cell
		covered by land.
		
		*/
		
		}else if( nline == 4 ){
			get_first_string_element( line, tmpc ) ;
			counts = load_ldMatrix( &cx, &cy, tmpc );			
			check_dims( nx, ny, cx, cy );
		
			if( get_first_string_element( line, tmpc ) != 0 ){
				weight_intgrl_by_cell_area = TRUE ;
				cell_area = load_ffMatrix( &cx, &cy, tmpc );
				check_dims( nx, ny, cx, cy );
			}
	
			if( get_first_string_element( line, tmpc ) != 0 ){
				weight_intgrl_by_land_area = TRUE ;
				land_area = load_ffMatrix( &cx, &cy, tmpc );
				check_dims( nx, ny, cx, cy );
			}
		
	
		/*
		Fifth line - 
		  nprops worth of output filenames
		  for the posterior means
		*/
		
		}else if( nline == 5 ){
			n=0;
			while( get_first_string_element( line, tmpc ) ){
				if( ( out_means[ n ] = fopen( tmpc, "w" ) ) == NULL ){
					fprintf( stderr, "%s: cannot open file %s for output\n", *argv, tmpc );
					exit( EXIT_FAILURE );
				}
				n++; 
			}
			if( n != n_props ){
				fprintf( stderr, "%s: must specify as many ouput mean files as there are input proportions\n", *argv );
				exit( EXIT_FAILURE );
			}
		
	
		/*
		Sixth line - 
		  nprops worth of output filenames
		  for the posterior means
		*/
		
		}else if( nline == 6 ){
			n=0;
			while( get_first_string_element( line, tmpc ) ){
				if( ( out_stdvs[ n ] = fopen( tmpc, "w" ) ) == NULL ){
					fprintf( stderr, "%s: cannot open file %s for output\n", *argv, tmpc );
					exit( EXIT_FAILURE );
				}
				n++; 
			}
			if( n != n_props ){
				fprintf( stderr, "%s: must specify as many ouput stdv files as there are input proportions\n", *argv );
				exit( EXIT_FAILURE );
			}
		
		/*
		Line seven and up -
		  Each line should contain one flux file per lc/pft proportion
		  file, followed by output mean and stdv file names.
		*/
		
		}else if( nline >= 7 ){
			
			/*
			Read in one flux file per lc/pft proportion map (i.e. single pft flux maps)
			*/
			
			for( i=0; i<n_props; i++ ){
				if( get_first_string_element( line, tmpc ) == 0 ){
					fprintf( stderr, "%s: expecting one flux file for each LC proportion (%d) on line %d\n", *argv, n_props, nline );
					exit( EXIT_FAILURE );
				}
				flux[ n_fluxes ][ i ] = load_ffMatrix( &nx, &ny, tmpc );
				check_dims( nx, ny, cx, cy );	
			}
			
			
			/*
			Open a file for output of flux means
			*/
			
			if( get_first_string_element( line, tmpc ) == 0 ){
				fprintf( stderr, "%s: expecting a file for output of mean flux on line %d\n", *argv, nline );
				exit( EXIT_FAILURE );
			}
			if( ( out_means[ n_props+n_fluxes ] = fopen( tmpc, "w" ) ) == NULL ){
				fprintf( stderr, "%s: cannot open file %s for output\n", *argv, tmpc );
				exit( EXIT_FAILURE );
			}
			
			/*
			Open a file for output of flux stdvs
			*/
			
			if( get_first_string_element( line, tmpc ) == 0 ){
				fprintf( stderr, "%s: expecting a file for output of flux stdv on line %d\n", *argv, nline );
				exit( EXIT_FAILURE );
			}
			if( ( out_stdvs[ n_props+n_fluxes ] = fopen( tmpc, "w" ) ) == NULL ){
				fprintf( stderr, "%s: cannot open file %s for output\n", *argv, tmpc );
				exit( EXIT_FAILURE );
			}

			n_fluxes++ ;
			
		}		

	}
	
	fclose( ifp );

	if( nline < 6 ){
		fprintf( stderr, "%s: did not find the correct number of input lines.\nExpecting 6+, found %d\n", *argv, nline );
		exit( EXIT_FAILURE );
	
	}



	/*
	------------------------
	Allocate memory for, and
	calculate the priors
	------------------------
	*/
	
	for( i=0; i<n_props; i++ ){
		priors[i] = (double *)malloc( sizeof( double ) * nx * ny );
		conditional_ar( priors[i], lc[i], msk, length, nx, ny );
	}

	
		
	/*
	-----------------------------------
	Allocate space to hold the forward
	conditional, joint and posterior 
	probabilites and work space vectors
	-----------------------------------
	*/
	
	
	fwd_cond_prob = (double *)malloc( sizeof( double ) * n_props * n_props );
	joint_prob = (double *)malloc( sizeof( double ) * n_props * n_props );
	work_int = (unsigned int *)malloc( sizeof( unsigned int ) * n_props );
	post_prob = (double *)malloc( sizeof( double ) * n_props * n_props );
	tmp_vec = (long *)malloc( sizeof( long ) * n_props );
	work1 = (double *)malloc( sizeof( double ) * n_props );
	work2 = (double *)malloc( sizeof( double ) * n_props );


	/*
	----------------------------------
	Count the number of sites and make 
	masked vectors of the porportions
	and the priors
	
	Also read in grid cell areas if 
	present - otherwise set all areas
	to unity.
	----------------------------------
	*/

	n_sites = count_sites( msk, nx, ny ) ;
	
	counts_vec = (long *)malloc( sizeof( long ) * n_sites );
	extract_masked_long( counts_vec, counts, msk, nx, ny );


	/*site areas*/
	cell_area_vec = (double *)malloc( sizeof( double ) * n_sites );
	if( weight_intgrl_by_cell_area == TRUE ){
		extract_masked( cell_area_vec, cell_area, msk, nx, ny );
		free( cell_area );
	}else{
		for( i=0; i<n_sites; i++ )
			cell_area_vec[ i ] = 1.0 ;
	}
	
	land_area_vec = (double *)malloc( sizeof( double ) * n_sites );
	if( weight_intgrl_by_land_area == TRUE ){
		extract_masked( land_area_vec, land_area, msk, nx, ny );
		free( land_area );
	}else{
		for( i=0; i<n_sites; i++ )
			land_area_vec[ i ] = 1.0 ;
	}
	
	
	for( i=0; i<n_props; i++ ){
		lc_vec[i] = (double *)malloc( sizeof( double ) * n_sites );
		extract_masked( lc_vec[i], lc[i], msk, nx, ny );
		
		priors_vec[i] = (double *)malloc( sizeof( double ) * n_sites );
		extract_masked( priors_vec[i], priors[i], msk, nx, ny );
		
		sum_post[i] = (double *)malloc( sizeof( double ) * n_sites );
		sum_post_sqr[i] = (double *)malloc( sizeof( double ) * n_sites );
		
		/*zero output arrays*/
		for( j=0; j<n_sites; j++ ){
			*( sum_post[i]+j ) = 0. ;			
			*( sum_post_sqr[i]+j ) = 0. ;			
		}
	}

	/*
	----------------------------
	Allocations for flux vectors
	----------------------------
	*/

	for( i=0; i<n_fluxes; i++ ){
	
		for( j=0; j<n_props; j++ ){
			flux_vec[i][j] = (double *)malloc( sizeof( double ) * n_sites );
			extract_masked( flux_vec[i][j], flux[i][j], msk, nx, ny );
		}
		
		sum_post[i+n_props] = (double *)malloc( sizeof( double ) * n_sites );
		sum_post_sqr[i+n_props] = (double *)malloc( sizeof( double ) * n_sites );		

		/*zero output arrays*/
		for( j=0; j<n_sites; j++ ){
			*( sum_post[i+n_props]+j ) = 0. ;			
			*( sum_post_sqr[i+n_props]+j ) = 0. ;			
		}
		
		if( i==0 ){
			out_flux_means_vec = (double *)malloc( sizeof( double ) * n_sites );
			out_flux_stdvs_vec = (double *)malloc( sizeof( double ) * n_sites );
		
		}


	}


	/*
	------------------------
	Monte Carlo simulations:
	------------------------	
	*/


	for( i=0; i<n_sim; i++ ){
	
		/*
		--------------------------
		Zero the spatial integrals
		--------------------------
		*/
	
		for( j=0; j<n_fluxes+n_props; j++ )
			spatial_integral[ j ] = 0.0 ;
		for( j=0; j<n_fluxes*n_props; j++ )
			spatial_integral_pftflx[ j ] = 0.0 ;	
	
		/*
		-----------------------------------
		Generate the forward conditional
		probabilites based on the confusion
		matrix & the Dirichlet distribution
		------------------------------------
		*/
	
		for( j=0; j<n_props; j++ ){
			for( k=0; k<n_props; k++ ){
			 	*( work1 + k ) = *( cfm + k * n_props + j );
				*( work2 + k ) = 0. ;
			}
			gsl_ran_dirichlet( rnd_num_gen, n_props, work1, work2 );
			for( k=0; k<n_props; k++ ) *( fwd_cond_prob + k * n_props + j ) = *( work2 + k ) ;
		}		
	
	
		/*
		---------------------------------
		Calculate the joint probabilities
		for all sites, and from this the
		posterior distributions.
		---------------------------------
		*/
	
		
		for( m=0; m<n_sites; m++ ){
			
			/*joint probs, note multiplication
			through by dcorr and then priors*/
			
			for( j=0; j<n_props; j++ ){
				
				for( k=0; k<n_props; k++ ) 
					*( work1 + k ) = *( fwd_cond_prob + k * n_props + j ) * dcorr;
				
				gsl_ran_dirichlet( rnd_num_gen, n_props, work1, work2 );
				
				for( k=0; k<n_props; k++ ) 
					*( joint_prob + k * n_props + j ) = *( work2 + k ) * *( priors_vec[j]+m );
			
			
				/*also zero the tmp vector
				here for efficiency*/
				*(tmp_vec+j) = 0. ;
			
			}
			
			/*posterior probs*/
			
			for( j=0; j<n_props; j++ ){
			
				sum_jprob = 0. ;
				for( k=0; k<n_props; k++ ) sum_jprob += *( joint_prob+j*n_props+k ) ;
				for( k=0; k<n_props; k++ ){
					
					*( post_prob+j*n_props+k ) = *( joint_prob+j*n_props+k ) / sum_jprob ;

					
					/*work vectors initialised here
					to avoid using another loop*/

					*( work1+k ) = *( post_prob+j*n_props+k ) ;
					*( work_int+k ) = 0. ;
				}
				
				
				gsl_ran_multinomial( rnd_num_gen, n_props, round( *( counts_vec+m )* *( lc_vec[j]+m ) ), work1, work_int );
				
				
				for( k=0; k<n_props; k++ ) *( tmp_vec+k ) += *( work_int+k ) ;
				
			}
			
			
			
			/*calculate outputs*/
			for( j=0; j<n_props; j++ ){
			
				spatial_integral[ j ] += *( tmp_vec+j ) / (float) *( counts_vec+m ) * cell_area_vec[ m ] * land_area_vec[ m ] ;
			
				*( sum_post[j]+m ) += *( tmp_vec+j ) / (float) *( counts_vec+m ) ;
				*( sum_post_sqr[j]+m ) += ( *( tmp_vec+j ) / (float) *( counts_vec+m ) ) * ( *( tmp_vec+j ) / (float) *( counts_vec+m ) ) ;			
	
			}
			for( j=0; j<n_fluxes; j++ ){
				flx=0.;
				for( k=0; k<n_props; k++ ){
					flx += *( flux_vec[j][k]+m ) * *( tmp_vec+k ) / (float) *( counts_vec+m ) ;
					spatial_integral_pftflx[ n_props*j+k ] += *( flux_vec[j][k]+m ) * *( tmp_vec+k ) / (float) *( counts_vec+m ) * cell_area_vec[ m ] * land_area_vec[ m ] ;
				}
							
				spatial_integral[ n_props+j ] += flx * cell_area_vec[ m ] * land_area_vec[ m ] ;
				
				*( sum_post[n_props+j]+m ) += flx ;
				*( sum_post_sqr[n_props+j]+m ) += flx * flx ;			
				
			}
			
			
		}/*m, sites*/	
		
		
		
		/*
		------------------------
		Output spatial integrals
		------------------------
		*/
		
		for( j=0; j<n_fluxes+n_props; j++ )
			fprintf( stdout, "%f ", spatial_integral[ j ] ) ;
		
		for( j=0; j<n_fluxes*n_props; j++ )
			fprintf( stdout, "%f ", spatial_integral_pftflx[ j ] ) ;
			
		fprintf( stdout, "\n" );
		
		
		
		
	}/*i, simulations*/



	/*
	------------------------------------
	calculate means and stdv and output.
	************************************
	N.B. uses the lc_vec and prior_vec
	array as they're no longer needed.
	************************************
	***  lc_vec = means              ***
	***  priors_vec = stdvs          ***
	************************************
	------------------------------------
	*/
	
	for( m=0; m<n_sites; m++ ){
		for( j=0; j<n_props; j++ ){	
	
			*( lc_vec[j]+m ) = *( sum_post[j]+m ) / (float) n_sim ;
			*( priors_vec[j]+m ) = ( *( sum_post_sqr[j]+m ) - ( *( sum_post[j]+m ) * *( sum_post[j]+m ) ) / (float) n_sim ) / (float) n_sim ;
		}
	}

	
	for( j=0; j<n_props; j++ ){
		output_to_ascii_matrix( lc_vec[j], msk, nx, ny, out_means[j] ) ;
		output_to_ascii_matrix( priors_vec[j], msk, nx, ny, out_stdvs[j] ) ;
	}

	/*
	re-use the same vectors for the output
	of the flux statistics (if required)
	*/

	for( j=0; j<n_fluxes; j++ ){	
		for( m=0; m<n_sites; m++ ){
	
			*( out_flux_means_vec+m ) = *( sum_post[n_props+j]+m ) / (float) n_sim ;
			*( out_flux_stdvs_vec+m ) = sqrt( ( *( sum_post_sqr[n_props+j]+m ) - ( *( sum_post[n_props+j]+m ) * *( sum_post[n_props+j]+m ) ) / (float) n_sim ) / (float) n_sim );
		}

		output_to_ascii_matrix( out_flux_means_vec, msk, nx, ny, out_means[n_props+j] ) ;
		output_to_ascii_matrix( out_flux_stdvs_vec, msk, nx, ny, out_stdvs[n_props+j] ) ;
	}


	/*
	--------------
	Free up memory
	--------------
	*/
	
	free( msk );
	free( cfm );
	free( work1 );
	free( work2 );
	free( work_int );
	free( tmp_vec );
	free( post_prob );
	free( joint_prob );
	free( fwd_cond_prob );
	free( cell_area_vec );
	free( land_area_vec );

	
	for( i=0; i<n_props; i++ ){
		free( lc[i] );
		free( lc_vec[i] );
		free( priors[i] );
		free( priors_vec[i] );
		
		fclose( out_means[i] );
		fclose( out_stdvs[i] );
	}
	
	
	if( n_fluxes > 0 ){
		free( out_flux_means_vec );
		free( out_flux_stdvs_vec );
	}
	
	for( i=0; i<n_fluxes; i++ ){
	
		free( sum_post[n_props+i] );
		free( sum_post_sqr[n_props+i] );
	
		fclose( out_means[n_props+i] );
		fclose( out_stdvs[n_props+i] );
	
		for( k=0; k<n_props; k++ ){
			free( flux[i][k] );
			free( flux_vec[i][k] );
		}
	
	}
	
	/*
	--------------
	return success
	--------------
	*/
	
	return( EXIT_SUCCESS );

}



int conditional_ar( double *outpt, double *inpt, long *msk, int length, int nx, int ny )
/*======================================================
conditional autoregressive model.

This is used to set up the priors and is basically
a low pass filter with size 1+2*length
======================================================*/
{

	
	int i, ii, j, jj, n;
	double sum ;
		
	/*zero the ouput*/

	for( i=0; i<ny; i++ )
		for( j=0; j<nx; j++ )
			*outpt = 0.0 ;

	
	/*low pass the inpt matrix*/

	for( i=0; i<ny; i++ ){
	 for( j=0; j<nx; j++ ){
	  
	  if( *( msk+i*nx+j ) == 0 ) continue;
	  
	  n = 0 ;
	  sum = 0.0 ;
	  
	  for( ii=-length; ii<=length; ii++ ){
	    
	    if( (i+ii<0)||(i+ii>=ny) ) continue;
	    
	    for( jj=-length; jj<=length; jj++ ){
	    
	      if( (j+jj<0)||(j+jj>=nx) ) continue;
	      if( *( msk+(i+ii)*nx+(j+jj) ) == 0 ) continue;
	      
	      sum+=*( inpt+(i+ii)*nx+(j+jj) ) ;
	      n++ ;
	   
	    }
	  }
	 
	  *( outpt +i*nx+j ) = sum / (float) n ;
	 
	 }
	}


	return( 1 );
}



int count_sites( long *msk, int nx, int ny )
/*return the number of non-zero elements
in the matrix/image pointed to by msk*/
{
	int n=0, i ;
	
	for( i=0;i<nx*ny; i++ ) 
		if ( *( msk+i ) != 0 ) n++ ;
	
	return( n );

}



int extract_masked( double *vec, double *mat, long *msk, int nx, int ny )
/*extract the none masked elements of mat and place them in vec
(there must be enough space alloacted in vec to do this). The total
number of unmasked items is returned.*/
{


	int n=0, i ;
	
	for( i=0;i<nx*ny; i++ ) 
		if ( *( msk+i ) != 0 ) *( vec+n++ ) = *( mat+i ) ;
	
	return( n );

}


int extract_masked_long( long *vec, long *mat, long *msk, int nx, int ny )
/*A version of the above for processing long int data*/
{


	int n=0, i ;
	
	for( i=0;i<nx*ny; i++ ) 
		if ( *( msk+i ) != 0 ) *( vec+n++ ) = *( mat+i ) ;
	
	return( n );

}



void output_to_ascii_matrix( double *vec, long *msk, int nx, int ny, FILE *stream )
/*output a vector to a specified stream, placing 
values where there is no mask indicated*/
{

	int i, j, n=0 ;
	
	for( i=0;i<ny; i++ ){
		for( j=0;j<nx; j++ ){
			if( *(msk+i*nx+j) != 0 ){
				fprintf( stream, "%f ", *(vec+n++) );
			}else{
				fprintf( stream, "%f ", 0. );
			}

 
		}
		fprintf( stream, "\n" ) ;
	}
	
		
	return ;

}



void cl_parser( argc, argv, d, l, n )
/* Get options from the command line, these are used to control the
correlation parameter (dcorr) the car lenght (length) and the number
of monte carlo simulation to perfrom (n_sim)*/
int argc;
char **argv;
double *d ;
int *l ;
int *n ;
{

	int i;

	/*
	Sweep the line for a usage request first. This is required because the
	next sweep only goes up to argc-1 (to allow for the control file to be 
	placed at the end of the command line). However we also want the user
	to be able to just use the -u flag so we need to check for that first.
	*/

	for( i=1; i<argc; i++ ) if( !strncmp( argv[ i ], "-u", 2 ) ) usage( *argv ); 

	for( i=1; i<(argc-1); i++ ){

		if( *argv[ i ] == '-' ){

			/**/ if( !strncmp( argv[ i ], "-d", 2 ) ) *d = atof( argv[ ++i ] ); 
			else if( !strncmp( argv[ i ], "-l", 2 ) ) *l = atoi( argv[ ++i ] ); 
			else if( !strncmp( argv[ i ], "-n", 2 ) ) *n = atoi( argv[ ++i ] ); 
			else{
			
				fprintf( stderr, "%s: unknown option on command line: %s\n", argv[ 0 ], argv[ i ] );
				fprintf( stderr, "(use the option -u to see brief usage instructions)\n" );
				exit( EXIT_FAILURE );

			}
		
		}
	
	}

}



void usage( char *bin_name )
{

	fprintf( stderr, "\nusage: %s [options] control_file.txt\n\n", bin_name );
	
	fprintf( stderr, "The control file must be the last word on the command line.\n" );
	
	fprintf( stderr, "Options are:\n" );
	fprintf( stderr, "-d %%f        set the correlation factor to %%f\n" );
	fprintf( stderr, "-l %%d        set the length for the conditional auto regressive model to %%d\n" );
	fprintf( stderr, "-n %%d        set the number of simulations to %%d (the higher the better)\n" );
	fprintf( stderr, "-u           display this message\n" );

	fprintf( stderr, "\nThe format of the control file is a series of ascii text lines each specifying user\n" );
	fprintf( stderr, "defined input and output filenames. Anything after a # character is treated as a comment.\n" );
	fprintf( stderr, "Blank lines, or lines consisting only of a comment are ignored. Excluding these\n" );
	fprintf( stderr, "the following information is expected, in order:\n\n" );
	
	fprintf( stderr, "Line 1: a list of N files containing the input LC porportions\n" );
	fprintf( stderr, "Line 2: a mask file containing 0 where a mask is to be applied and 1 elsewhere\n" );
	fprintf( stderr, "Line 3: a file containing an NxN confusion matrix\n" );
	fprintf( stderr, "Line 4: a file containing the number of \"counts\"\n" );
	fprintf( stderr, "Line 5: a list of N files to contain the output means\n" );
	fprintf( stderr, "Line 6: a list of N files to contain the output stdvs\n" );

	fprintf( stderr, "\nAside from the confusion matrix (which is an NxN ascii matrix) all the input files are plain, headerless\n" );
	fprintf( stderr, "ascii text files containg a single matrix, each of which must have the same dimensions as the rest.\n" );


	exit( EXIT_FAILURE );
	return ;
}
