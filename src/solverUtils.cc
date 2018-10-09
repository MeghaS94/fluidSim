/*
    Class/ file that contains all the util functions and the main solver code
*/

// Density increases due to sources
// This function adds sources to the current state of the density grid. The soures are input as an array 's'ss

#include "../includes/solverUtils.h"
#include <math.h>

void add_source ( int N, float * x, float * s, float dt )
{
    int i, size=(N+2)*(N+2);
    for ( i=0 ; i<size ; i++ )
    {
        // you can remove dt!
        x[i] += dt*s[i];
    }
}

// diffusion calculation at each grid cell
// this leads to unstable states
void diffuse_bad ( int N, int b, float * x, float * x0, float diff, float dt )
{
    int i, j;
    float a=dt*diff*N*N;
    for ( i=1 ; i<=N ; i++ ) 
    {
        for ( j=1 ; j<=N ; j++ ) 
        {
            x[IX(i,j)] = x0[IX(i,j)] + a*(x0[IX(i-1,j)] + x0[IX(i+1,j)] + x0[IX(i,j-1)] + x0[IX(i,j+1)] - 4*x0[IX(i,j)]);
        }
    }
    // boundary conditions
    // set_bnd ( N, b, x );
}

// stable method : find the densities which when diffused backward in time, yield the densities we started with.
// x0[IX(i,j)] = x[IX(i,j)] 
// solve the resulting system of linear equation using the gauss-seidel relaxation technique
// better method to solve the system og lin eqs : conjugate gradiet method
void diffuse ( int N, int b, float * x, float * x0, float diff, float dt)
{
    int i, j, k;
    float a=dt*diff*N*N; // rate * time step
    a= 0.0001;
    // cout << "a is :" << a << "dt, diff,N are : " << dt << ", " << diff << ", " << N << endl;
    for ( k=0 ; k<20 ; k++ ) 
    {
        for ( i=1 ; i<=N ; i++ ) 
        {
            for ( j=1 ; j<=N ; j++ ) 
            {
                if (!isnan(x[IX(i,j)]) and !isnan( x0[IX(i,j)]) and  !isnan( x[IX(i-1,j)]) and !isnan( x[IX(i+1,j)]) and !isnan(x[IX(i,j-1)] ) and !isnan( x[IX(i,j+1)]) )
                {
                    x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)] +x[IX(i+1,j)] +x[IX(i,j-1)] +x[IX(i,j+1)]))/(1+4*a);
                    // std::cout << x[IX(i,j)] << ", "<< x0[IX(i,j)] << "," <<  x[IX(i-1,j)] << ", "<< x[IX(i+1,j)] << ", " << x[IX(i,j-1)] << ","<< x[IX(i,j+1)] << std::endl;
                }
               
            }
        }
    // boundary condition handling 
    set_bnd ( N, b, x );
    }

    // int index = 0;
	// 	for (int i=0;i<N+2;i++)
	// 	{
	// 		for (int j=0;j<N+2;j++)
	// 		{
	// 			gridPoints[index] = float(i);///float(N);
	// 			gridPoints[index+1] = float(j);///float(N);
	// 			gridPoints[index+2] = 0.0;  
	// 			gridPoints[index+3] = x[IX(i,j)];
	// 			// std::cout << gridPoints[index]<< std::endl;
	// 			index = index + 4;
	// 		}
	// 	}
}

// find where the cell center would be at the previoius time step
// linear backtrace
void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*N;
    for ( i=1 ; i<=N ; i++ ) 
    {
        for ( j=1 ; j<=N ; j++ ) 
        {
            // distance travelled from i,j, backwards along the velocity field
            x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)]; // distance = u * dt
            // boundary conditions
            if (x<0.5)
            {
                x=0.5;
            }
            if (x>N+0.5)
            {
                x=N+ 0.5;
            }  
            if (y<0.5)
            {
                y=0.5;
            }  
            if (y>N+0.5)
            {
                y=N+ 0.5;
            }
            
            i0=(int)x; i1=i0+1; 
            j0=(int)y; j1=j0+1;
            
            s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
            // bilinear interpolation
            // can be changed to hermite interpolation
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
                            s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
        }
    }
    // boundary condition
    set_bnd ( N, b, d );
}

void set_bnd( int N, int b, float * x )
{
    int i;
    for ( i=1 ; i<=N ; i++ ) 
    {
        x[IX(0,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
        x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
        x[IX(i,0)] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
    }
    x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
    x[IX(0 ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]);
    x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);
}
 
void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt, float * vel_pos, float *  dens_grid )
{
    // std::cout << "stepping density" << std::endl;
    add_source ( N, x, x0, dt );        
    SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt);
    SWAP ( x0, x );  
    advect ( N, 0, x, x0, u, v, dt );
    
    // update densities in the grid
    int index = 0;
    // int index1 = 0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			// vel_pos[index+6] = 
            dens_grid[index+3] = x[IX(i,j)];
			index = index + 4;
		}
	}    

}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt )
{
    add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
    SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
    SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
    // project ( N, u, v, u0, v0 );
    SWAP ( u0, u ); SWAP ( v0, v );
    advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
    // project ( N, u, v, u0, v0 );

    // vorticity confinement?
    // update velocities in the grid
}

