#include <math.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "omp.h"
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
    const int N = 500; //number of grid points in each direction

//set up grid of points in the domain 0<=x,y<=1
    double h = 1. / ((double) N);
//    printf("%f\n", h);
    double u[N + 1][N + 1]; //u(x,y)
    double f[N + 1][N + 1]; //f(x,y)
//    double g[N + 1][N + 1]; //g(x,y)

    float maximum; //maximum absolute difference between iterations of u(i,j).
    const float tolerance = 0.001; //tolerance for Gauss-Seidel method. When maximum (see above) between iterations of u is below this tolerance, stop the loop
    int count = 0;

    /*Boundary conditions:
     f(x,y) = 0
     u(x,y=0) = 100-200x
     u(x=0,y) = 100-200y
     u(x,y=1) = -100+200x
     u(x=1,y) = -100+200y
     */
    for (int i = 0; i < N + 1; i++) {
        for (int j = 0; j < N + 1; j++) {
            double x = i * h;
            double y = j * h;
            f[i][j] = 0;
            u[i][j] = 1.0; //set initial u = 1

            if (j == 0) {
                u[i][j] = 100.0 - 200.0 * x;
            } else if (j == N) {
                u[i][j] = -100.0 + 200.0 * x;
            }
            if (i == 0) {
                u[i][j] = 100.0 - 200.0 * y;
            } else if (i == N) {
                u[i][j] = -100.0 + 200.0 * y;
            }
        }
    }

//write initial guess for u(x,y) to file
    ofstream fout1("initParallel2.txt");

    if (fout1) {
        cout << "Writing initial u(x,y) to file 'initParallel2.txt'" << endl;
        for (int i = 0; i < N + 1; i++) {
            for (int j = 0; j < N + 1; j++) {
                fout1 << u[i][j] << " ";
            }
            fout1 << endl;
        }
    }

    int num = 8; //number of threads
    omp_set_num_threads(num);

    std::cout << "Starting parallel Gauss-Seidel method \n";

    stopwatch<std::milli,float> sw;

    sw.start();

//Loop:
    do {
        maximum = 0;
#pragma omp parallel for shared(u) reduction(max:maximum)
        for (int i = 1; i < N; i++) {
            for (int j = 1; j < N; j++) {
                const double previous = u[i][j];
                u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - pow(h, 2) * f[i][j]);
                const double difference = fabs(previous - u[i][j]);
                if (difference > maximum) maximum = difference;
            }
        }
        count++; //count number of iterations until convergence
    } while (maximum > tolerance && count < 10000);

    sw.stop();

    if (count < 10000) {
	std::cout << "Converged \n";
    	std::cout << "Number of iterations: " << count << "\n";
    	std::cout << "Elapsed Time: " << sw.count() << " ms\n";

	//write final u(x,y) to file
    	ofstream fout("uFinal_Parallel2.txt");

    	if (fout) {
       	 cout << "Writing final u(x,y) to file 'uFinal_Parallel2.txt'" << endl;
       	 for (int i = 0; i < N + 1; i++) {
       	     for (int j = 0; j < N + 1; j++) {
       	         fout << u[i][j] << " ";
       	     }
       	     fout << endl;
       	 }
    	}
    }
    else {
	std::cout << "Didn't converge within 1000 iterations \n" << endl;
    }
}
