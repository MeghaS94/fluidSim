#ifndef SOLVERUTILS_H    // To make sure you don't declare the function more than once by including the header 
#define SOLVERUTILS_H

#include <iostream>

using namespace std;

#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}
//array value accessor
#define IX(i,j) ((i)+(N+2)*(j))

void add_source ( int N, float * x, float * s, float dt );

// diffusion calculation at each grid cell
// this leads to unstable states
void diffuse_bad ( int N, int b, float * x, float * x0, float diff, float dt );

// stable method : find the densities which when diffused backward in time, yield the densities we started with.
// x0[IX(i,j)] = x[IX(i,j)] 
// solve the resulting system of linear equation using the gauss-seidel relaxation technique
void diffuse ( int N, int b, float * x, float * x0, float diff, float dt);

// I do not understand this fully!! check the implementation work it out on paper
// find where the cell center would be at the previoius time step
// linear backtrace
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt );

void project ( int N, float * u, float * v, float * p, float * div );

void setBoundaryCondition( int N, int b, float * x );
 
void densityStep ( int N, float * x, float * x0, float * u, float * v, float diff, float dt, float* gridPoints, float* gridDensities );

void velocityStep ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt, float* vel_pos );

#endif