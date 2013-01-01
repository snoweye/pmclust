#include <R.h>
#include <Rinternals.h>

/*
   Example:
   W = [ 1.0 2.0 3.0 4.0
         5.0 6.0 7.0 8.0 ]
   y = [ -1.0 -2.0 -3.0 -4.0 ]
   W_new = [ 0.0 0.0 0.0 0.0
             4.0 4.0 4.0 4.0 ]
*/

/*
   This function adds y to W by column and returns Z.
   All W, and y should be in double.
*/
SEXP W_plus_y(SEXP R_W, SEXP R_y, SEXP R_nrow, SEXP R_ncol){
	SEXP R_Z;
	double *C_Z, *C_W, *C_y;
	int *C_nrow, *C_ncol;
	int i, j;

	C_W = REAL(R_W);
	C_y = REAL(R_y);
	C_nrow = INTEGER(R_nrow);
	C_ncol = INTEGER(R_ncol);

	PROTECT(R_Z = allocVector(REALSXP, *C_nrow * *C_ncol));
	C_Z = REAL(R_Z);

	for(j = 0; j < *C_ncol; j++){
		for(i = 0; i < *C_nrow; i++){
			*C_Z = *C_W + *C_y;
			C_Z++;
			C_W++;
		}
		C_y++;
	}

	UNPROTECT(1);
	return(R_Z);
} /* End of W_plus_y(). */

