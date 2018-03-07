#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grid.hh"
#include "physical_constants.hh"

GRID::GRID()
{
}


//------------------------------------------------------------   
// Read in the model file in ascii format
//------------------------------------------------------------ 
void GRID::Read(const char *infile)
{

    int my_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );

    verbose = 0;

    // Only the master MPI task reads the model from disk.

    if( my_rank == 0 )
    {

        verbose = 1;

        FILE *in = fopen(infile,"r");
        if( in == NULL )
        {
            printf("ERROR:  Can't open setup file %s\n",infile);
            exit(1); 
        }
    
        // basic grid parameters

        fscanf(in,"%d %lf %lf\n",&n_x,&dv,&t_begin);
        n_zones = n_x*n_x*n_x;
        dx = dv*t_begin*DAY_TO_SEC;
        x_cen   = n_x/2.0*dx;
        t_begin = t_begin;
    
        // estimate size of grid
    
        // allocate memory

        rho = new double[n_zones];     // mass density (g cm^(-3))
        ni_frac = new double[n_zones]; // 56ni mass fraction
        mu_e = new double[n_zones];    // number of electrons per particle
        edep = new double[n_zones];    // energy deposited due to Compton scatters
        vel = new double[n_zones];     // velocity of zone
    
        // calculate sums

        double tmass   = 0;
        double ke      = 0;
        ni_mass        = 0;
        double vol = dx*dx*dx;
    
        // read size of file

        double x1,x2,x3;
        int ind = 0;
        for (int i=0;i<n_x;i++)
        {
            for (int j=0;j<n_x;j++)
            {
                for (int k=0;k<n_x;k++)
                {
                    fscanf(in,"%lf %lf %lf",&x1,&x2,&x3);
                    rho[ind]     = x1;
                    ni_frac[ind] = x2;
                    mu_e[ind]    = x3;
          
                    tmass    += rho[ind]*vol;
                    ni_mass  += rho[ind]*ni_frac[ind]*vol;
                    double vx = (i*dx - x_cen)/(t_begin*DAY_TO_SEC);
                    double vy = (j*dx - x_cen)/(t_begin*DAY_TO_SEC);
                    double vz = (k*dx - x_cen)/(t_begin*DAY_TO_SEC);
                    double vv = vx*vx + vy*vy + vz*vz;
                    vel[ind] = sqrt(vv);
                    ke    += 0.5*rho[ind]*vol*vv;
                    ind++;
                }
            }
        }

        fclose(in);
    
        printf("# Model read\n");
        printf("# Total mass = %.3e (%.3e Msun)\n",tmass,tmass/M_SUN);
        printf("# 56ni  mass = %.3e (%.3e Msun)\n",ni_mass,ni_mass/M_SUN);
        printf("# Kinetic E  = %.3e ergs\n",ke);

    }

    // Now we exchange the model to all the MPI tasks.

    MPI_Bcast( &n_zones, 1, MPI_INT   , 0, MPI_COMM_WORLD );
    MPI_Bcast( &n_x    , 1, MPI_INT   , 0, MPI_COMM_WORLD );
    MPI_Bcast( &dx     , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &dv     , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &x_cen  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &t_begin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &ni_mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    double* rho_buffer     = new double [ n_zones ];
    double* ni_frac_buffer = new double [ n_zones ];
    double* mu_e_buffer    = new double [ n_zones ];

    if( my_rank == 0 )
    {
        for( int i = 0; i < n_zones; ++ i )
        {
            rho_buffer    [ i ] = rho[i];
            ni_frac_buffer[ i ] = ni_frac[i];
            mu_e_buffer   [ i ] = mu_e[i];
        }
    }
    else
    {
        for( int i = 0; i < n_zones; ++ i )
        {
            rho_buffer    [ i ] = 0.0;
            ni_frac_buffer[ i ] = 0.0;
            mu_e_buffer   [ i ] = 0.0;
        }
    }

    MPI_Bcast( rho_buffer    , n_zones, MPI_FLOAT, 0, MPI_COMM_WORLD );
    MPI_Bcast( ni_frac_buffer, n_zones, MPI_FLOAT, 0, MPI_COMM_WORLD );
    MPI_Bcast( mu_e_buffer   , n_zones, MPI_FLOAT, 0, MPI_COMM_WORLD );

    if( my_rank != 0 )
    {
        rho = new double[n_zones];     // mass density (g cm^(-3))
        ni_frac = new double[n_zones]; // 56ni mass fraction
        mu_e = new double[n_zones];    // number of electrons per particle
        edep = new double[n_zones];    // energy deposited due to Compton scatters
        vel = new double[n_zones];     // velocity of zone
        for( int i = 0; i < n_zones; ++ i )
        {
            rho[i]     = rho_buffer    [ i ];
            ni_frac[i] = ni_frac_buffer[ i ];
            mu_e[i]    = mu_e_buffer   [ i ];
        }
    }

    delete [] rho_buffer;
    delete [] ni_frac_buffer;
    delete [] mu_e_buffer;
    
} 


//------------------------------------------------------------   
// Get the index of the i,j,k zone
//------------------------------------------------------------ 
int GRID::Get_Index(int i, int j, int k)
{
  return i*n_x*n_x + j*n_x + k;
}

//------------------------------------------------------------   
// Locate the zone index given real x,y,z coordinates
//------------------------------------------------------------ 
int GRID::Get_Zone_Index(double *x)
{
  int ix = round((x[0]+x_cen)/dx); 
  int iy = round((x[1]+x_cen)/dx);
  int iz = round((x[2]+x_cen)/dx);
  
  // check for off grid
  if ((ix < 0)||(ix >= n_x)) return -1;
  if ((iy < 0)||(iy >= n_x)) return -1;
  if ((iz < 0)||(iz >= n_x)) return -1;

  return Get_Index(ix,iy,iz);
}

//------------------------------------------------------------   
// Homologously expand grid by a factor e
//------------------------------------------------------------   
void GRID::Expand(double e)
{
  const double e3inv = 1.0/(e*e*e);
  dx = dx*e;
  x_cen = x_cen*e;
  for (int i=0;i<n_zones;i++)
    rho[i]= rho[i]*e3inv;
}
