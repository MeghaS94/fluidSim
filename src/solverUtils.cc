/*
    Class/ file that contains all the util functions and the main solver code
*/

// Density increases due to sources
// This function adds sources to the current state of the density grid. The soures are input as an array 's'ss

#include "../includes/solverUtils.h"
#include <math.h>
#include <algorithm>

void add_source ( int N, float * x, float * x_prev, float dt )
{
    int i, size=(N+2)*(N+2);
    for ( i=0 ; i<size ; i++ )
    {
        // you can remove dt!
        x[i] += dt*x_prev[i];
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
    // setBoundaryCondition ( N, b, x );
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
    for ( k=0 ; k<20 ; k++ ) 
    {
        for ( i=1 ; i<=N ; i++ ) 
        {
            for ( j=1 ; j<=N ; j++ ) 
            {
                if (!isnan(x[IX(i,j)]) and !isnan( x0[IX(i,j)]) and  !isnan( x[IX(i-1,j)]) and !isnan( x[IX(i+1,j)]) and !isnan(x[IX(i,j-1)] ) and !isnan( x[IX(i,j+1)]) )
                {
                    x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)] +x[IX(i+1,j)] +x[IX(i,j-1)] +x[IX(i,j+1)]))/(1+4*a);
                }
            }
        }
    setBoundaryCondition ( N, b, x );
    }
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
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)]) + s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
        }
    }
    setBoundaryCondition ( N, b, d );
}

// velocity should conserve mass -> flow in to the cell = flow out of the cell
// Hodge decomposition -> velocity field = mass conserving field + gradient field(gradient field is scalar.
// mass conserving field = current velocity field - gradient field.
// gradient field is calculated by solving the poissons equation

void project ( int N, float * u, float * v, float * p, float * div )
{
    int i, j, k;
    float h;
    h = 1.0/N;
    for ( i=1 ; i<=N ; i++ ) 
    {
        for ( j=1 ; j<=N ; j++ ) 
        {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)] + v[IX(i,j+1)]-v[IX(i,j-1)]);
            p[IX(i,j)] = 0;
        }
    }
    setBoundaryCondition ( N, 0, div ); setBoundaryCondition ( N, 0, p );
    for ( k=0 ; k<20 ; k++ ) 
    {
        for ( i=1 ; i<=N ; i++ ) 
        {
        for ( j=1 ; j<=N ; j++ ) 
            {
                p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)] + p[IX(i,j-1)] + p[IX(i,j+1)])/4;
            }
        }
    setBoundaryCondition ( N, 0, p );
    }
    for ( i=1 ; i<=N ; i++ ) 
    {
        for ( j=1 ; j<=N ; j++ ) 
        {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    setBoundaryCondition ( N, 1, u ); setBoundaryCondition ( N, 2, v );
}

void setBoundaryCondition( int N, int b, float * x )
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
 
void densityStep ( int N, float * x, float * x0, float * u, float * v, float diff, float dt, float * vel_pos, float *  dens_grid )
{
    add_source ( N, x, x0, dt );        
    SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt);
    SWAP ( x0, x );  
    advect ( N, 0, x, x0, u, v, dt );
    
    // update densities in the grid
    int index = 0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
            dens_grid[index+3] = x[IX(i,j)];
			index = index + 4;
		}
	}    

}

void velocityStep ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt, float * vel_pos )
{
    // source is added by the mouse interactions
    diffuse ( N, 1, u, u0, visc, dt );
    diffuse ( N, 2, v, v0, visc, dt );
    project ( N, u, v, u0, v0 );
    SWAP ( u0, u ); SWAP ( v0, v );
    advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
    project ( N, u, v, u0, v0 );

    // TO DO : add vorticity confinement
    // update velocities in the grid
    float maxu = -1000.0f, maxv = -1000.0f;
    float minu = 1000.0f, minv = 1000.0f;
    float modified_u, modified_v;
    int index = 0;
    for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
            // TO DO : normalise values before adding them
            vel_pos[index+3] = float(i) + u[IX(i,j)];
            vel_pos[index+4] = float(j) + v[IX(i,j)];
            vel_pos[index+5] = 0.0f;
			index = index + 6;
		}
	}   
}

