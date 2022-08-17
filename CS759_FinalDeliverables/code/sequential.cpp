/* Chad Bustard, Sequential finite difference method */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <sys/time.h>
//#include "omp.h"
using namespace std;

template<typename period = std::micro, typename rep = size_t>
class stopwatch {
public:
    using count_type = std::chrono::duration<rep,period>;

    stopwatch() {
        begin = end = std::chrono::high_resolution_clock::time_point::min();
    }
    inline void start(void) {
        begin = std::chrono::high_resolution_clock::now();
    }
    inline void stop(void) {
        end = std::chrono::high_resolution_clock::now();
    }
    inline rep count(void) {
        using namespace std::chrono;
        return duration_cast<count_type>(end - begin).count();
    }
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> begin;
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
};


int main() {
int N = 1000; //number of grid points in each direction
//set up grid of points in the domain 0<=x,y<=1
double h = 1./((double) N);
double u[N+1][N+1]; //u(x,y)
double f[N+1][N+1]; //f(x,y)
double g[N+1][N+1]; //g(x,y)


float previous; //holds previous iteration of u(i,j)
float maximum; //maximum absolute difference between iterations of u(i,j). 
float tolerance = 0.1; //tolerance for Gauss-Seidel method. When maximum (see above) between iterations of u is below this tolerance, stop the loop
float difference;
int count = 0; 

/*Boundary conditions:
	f(x,y) = 0	
	u(x,y=0) = 100-200x
	u(x=0,y) = 100-200y
	u(x,y=1) = -100+200x
	u(x=1,y) = -100+200y
*/
int i,j;
double x,y;
for(i = 0; i < N+1; i++) {
	for(j = 0; j < N+1; j++) {
		x = (double) i * h;
		y = (double) j * h;
		f[i][j] = 0;

		u[i][j] = 1.0; //set initial u = 1

		if (j == 0) 
		{ 
			u[i][j] = 100.0 - 200.0*x;
		}
		else if (j == N) 
		{
			u[i][j] = -100.0 + 200.0*x;
		}
		if (i==0) 
		{
			u[i][j] = 100.0 - 200.0*y;
		}
		else if (i == N) 
		{ 
			u[i][j] = -100.0 + 200.0*y;
		}	
	}
}

//write initial guess for u(x,y) to file
ofstream fout1("initial.txt");

if(fout1.is_open()) {
	cout << "Writing initial u(x,y) to file 'initial.txt'" << endl;
	for (i = 0; i < N+1; i++) {
		for (j = 0; j< N+1; j++) {
			fout1 << u[i][j] << " ";
		}
		fout1 << endl;
	}
	fout1.close();
}

std::cout << "Starting sequential Gauss-Seidel method \n";

stopwatch<std::milli,float> sw;

sw.start();

//Loop:
do
{
	maximum = 0;
	for(i = 1; i < N; i++) {
		for(j = 1; j < N; j++) {
			previous = u[i][j];
			u[i][j] = 0.25*(u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] - pow(h,2)*f[i][j]);
			difference = fabs(previous-u[i][j]); 
			if (difference > maximum) 
				maximum = difference;
		}
	}
	count ++; //count number of iterations until convergence
} while (maximum > tolerance && count < 10000);

sw.stop();

if (count < 10000) {
    std::cout << "Converged \n";
    std::cout << "Number of iterations: " << count << "\n";
    std::cout << "Elapsed Time: " << sw.count() << " ms\n";



//write final u(x,y) to file
ofstream fout("uFinal.txt");

if(fout) {
	cout << "Writing final u(x,y) to file 'uFinal.txt'" << endl;
	for (i = 0; i < N+1; i++) {
		for (j = 0; j< N+1; j++) {
			fout << u[i][j] << " ";
		}
		fout<<endl;
	}
	fout.close();
}
}
else {
	std::cout << "Didn't converge within 10000 iterations \n" << endl;
}
}

