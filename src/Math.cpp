#include "precompiled.h"
#include "Math.h"
#include "utility.h"
#include "integer_types.h"
#include "exceptions.h"

namespace fipster {

	

	void interpolate::test(){
		vector<double> x,y;
		vector<vector<double>::iterator> c;

		//prepare x, y and c
		x.resize(100,10000000);
		int offset=10;
		auto start=x.begin()+10;
		auto stop=x.end()-10;
		for(auto i=start;i!=stop;i++)
			*i=(i-start)*(i-start);
		std::generate_n(std::back_inserter(y),offset,[](){return -100;});
		std::generate_n(std::back_inserter(y),stop-start,rand);
		
		double a=y[offset],b=oneD(start,stop,y,offset,*start);
		FIPSTER_ASSERT(a==b);
		try{
			//check if the following throws and thus never reach the assertion
			oneD(start,stop,y,offset,*start-0.1);
			FIPSTER_ASSERT(0);
		}catch(...){}

		auto i=start+(stop-start)/2;
		a=y[offset+i-start],b=oneD(start,stop,y,offset,*i);
		FIPSTER_ASSERT(a==b);
		if(i+1!=stop && i!=start)
			for(int j=0;j<2;j++){
				double n=2*j-1+*i;
				a=y[offset+i-start-1+j]+(y[offset+i-start+j]-y[offset+i-start-1+j])/(*(i+j)-*(i-1+j))*(n-*(i-1+j));
				b=oneD(start,stop,y,offset,n);		
				FIPSTER_ASSERT(a==b);
			}

		a=y[offset+stop-start-1],b=oneD(start,stop,y,offset,*(stop-1));
		FIPSTER_ASSERT(a==b);
		try{
			//check if the following throws and thus never reach the assertion
			oneD(start,stop,y,offset,*(stop-1)+0.1);
			FIPSTER_ASSERT(0);
		}catch(...){}
	}


//FROM http://people.sc.fsu.edu/~jburkardt/cpp_src/subset/subset.html
//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
}
//****************************************************************************80


//Returns Ln[Gamma[xx]]
double gammln(const double xx) {
	int j;
	double x,tmp,y,ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}

//Returns n! as double
double factorial(const uint n) { 
	static double a[171];
	static bool init=true;
	if (init) {
		init = false;
		a[0] = 1.;
		for (int i=1;i<171;i++) a[i] = i*a[i-1];
	}
	if (n > 170) throw("factrl out of range");
	return a[n];
}

//Returns Ln[n!]
double factln(const int n) {
	static const int NTOP=2000;
	static double a[NTOP];
	static bool init=true;
	if (init) {
		init = false;
		for (int i=0;i<NTOP;i++) a[i] = gammln(i+1.);
	}
	if (n < 0) throw("negative arg in factln");
	if (n < NTOP) return a[n];
	return gammln(n+1.);
}

//Returns the binomial coefficient N over K as double
double binom(const int n, const int k) {
	if (n<0 || k<0 || k>n) throw("bad args in binom");
	if (n<171) return floor(0.5+factorial(n)/(factorial(k)*factorial(n-k)));
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}


//round to nearest integer represented as double
double round(double r) {
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

//calculates the relative difference of two floating point numbers
	double reldiff(double a, double b, bool abs_if_zero=false //use abs difference, if b is zero
								 ){
		if(abs_if_zero && b==0)
			return a;
	return 2.0*abs(b-a)/(abs(a)+abs(b));
}

//
//double maxreldiff( double2& a, double2 &b )
//{
//	if(a.length2d!=b.length2d)
//		return -1.0;
//
//	double r=0,t;
//	for(int i = 0; i<a.length ; ++i ) 
//		for(int j=0; j < a.cols; ++j) 
//	{
//		t = reldiff(a(i)(j),b(i)(j));
//		if(t>0.00001){
//			a.trace(toS(t));
//			a.trace(toS(i)+", "+toS(j)+": "+toS(a(i)(j))+" "+toS(b(i)(j)));
//
//		}
//		r=max(r,t);
//
//	}
//	return r;
//}

}
