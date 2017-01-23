
#include <tuple>
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <ctime>

using namespace std;

typedef vector<double> vd;

//calculates z = alpha*x + beta*y
void multiply_add(vd& z, const vd& x, const vd& y,double alpha, double beta){

	bool	alpha_zero	= alpha==0, 
			alpha_one	= alpha==1,
			beta_zero	= beta==0,
			beta_one	= beta==1;

	for(vd::size_type i=0;i<z.size();i++)
		z[i] = (alpha_zero?0:(alpha_one?1:alpha))*x[i] + 
				(beta_zero?0:(beta_one?1:beta))*y[i];
}

//calculates z = alpha*x + beta*y
void multiply_add2(vd& z, const vd& x, const vd& y,double alpha, double beta){

	bool	alpha_zero	= alpha==0, 
			alpha_one	= alpha==1,
			beta_zero	= beta==0,
			beta_one	= beta==1;

	for(vd::size_type i=0;i<z.size();i++)
		z[i] = alpha*x[i] + beta*y[i];
}


template<char flags>
void abc(vd& z, const vd& x, const vd& y,double alpha, double beta){

	bool	alpha_zero	= !!(flags & 1), 
			alpha_one	= !!(flags & 1<<1),
			beta_zero	= !!(flags & 1<<2),
			beta_one	= !!(flags & 1<<3);

	for(vd::size_type i=0;i<z.size();i++)
		z[i] = (alpha_zero?0:(alpha_one?1:alpha))*x[i] + 
				(beta_zero?0:(beta_one?1:beta))*y[i];
}

//calculates z = alpha*x + beta*y
void multiply_add3(vd& z, const vd& x, const vd& y,double alpha, double beta){

	bool	alpha_zero	= alpha==0, 
			alpha_one	= alpha==1,
			beta_zero	= beta==0,
			beta_one	= beta==1;

	char flags;
	flags=alpha_zero+(alpha_one<<1)+(beta_zero<<2)+(beta_one<<3);
	switch(flags){
	case 0:
		abc<0>(z,x,y,alpha,beta);break;
	case 1:
		abc<1>(z,x,y,alpha,beta);break;
	case 2:
		abc<2>(z,x,y,alpha,beta);break;
	case 3:
		abc<3>(z,x,y,alpha,beta);break;
	case 4:
		abc<4>(z,x,y,alpha,beta);break;
	case 5:
		abc<5>(z,x,y,alpha,beta);break;
	case 6:
		abc<6>(z,x,y,alpha,beta);break;
	case 7:
		abc<7>(z,x,y,alpha,beta);break;
	case 8:
		abc<8>(z,x,y,alpha,beta);break;
	}
}

/* timings
3: 0,216	7,48	2,005
2: 0,305	14,6	2,644	
1: 0,853
*/

int testmain()
{
	int N=1000000;
	vd z(N),y(N),x(N);
	generate_n(x.begin(),N,rand);
		generate_n(y.begin(),N,rand);
		bool a;
		cin>>a;
	clock_t t=0;
		t-=clock();
	for(int i=0;i<3000;i++)
	{
		double beta = (rand()>RAND_MAX/2) ^ a ? 1 : 0;
		double alpha = (rand()>RAND_MAX/2) ^ a ? 1 : 0;
		multiply_add3(z,x,y,alpha,beta);
	}
		t+=clock();
	cout<<setprecision(10)<<(double(t)/CLOCKS_PER_SEC)<<endl;
	exit(0);
	return 0;
}

