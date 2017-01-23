#pragma once

//(Johannes Gerer) ON 20.Jan.2011 FROM  http://hesperia.gsfc.nasa.gov/~schmahl/nnls/
//or from http://www.netlib.org/toms/587?


//Computes non-negative least squares (Check nnls.c for arguments)
// minimize: | D.x - b | with respect to x >= 0
// a contains D(m,n) in column-major ordering, b(m), x(n)
// rnorm will contain the residuum norm
// double working w(n + m) 
// integer working index(n)
/*  RETURNS:   1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
/*             2     THE DIMENSIONS OF THE PROBLEM ARE BAD. */
/*                   EITHER M .LE. 0 OR N .LE. 0. */
/*             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. */
int nnls_c(double* a, const int* m, const int* n, double* b, 
	double* x, double* rnorm, double* w, int* index);

