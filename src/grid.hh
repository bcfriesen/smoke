#ifndef _GRID_H
#define _GRID_H

//-------------------------------------------------
// Class to hold all the zones
//-------------------------------------------------
class GRID
{

public:

  double *rho;     // mass density (g cm^(-3))
  double *ni_frac; // 56ni mass fraction
  double *mu_e;    // number of electrons per particle
  double *edep;    // energy deposited due to Compton scatters
  double *vel;     // velocity of zone

private:

  int verbose;
  int n_zones;      // total number of zones
  int n_x;          // number of zones in x-dimension
  double dx;        // linear size of zone in cm
  double dv;        // size of zone in velocity coords (cm/s)
  double x_cen;     // location of grid center
  double t_begin;   // start time of the model
  double ni_mass;   // total 56ni mass


public:

  GRID();


  void Read(const char *infile);
  int  Get_Zone_Index(double *x);
  void Expand(double);
  
  int    Get_Index(int i, int j, int k);
  double Get_Nickel_Mass(int i) { return dx*dx*dx*rho[i]*ni_frac[i]; }
  double Get_Nickel_Mass()      { return ni_mass; }
  int    Get_n_zones() {return n_zones;}
  int    Get_n_x()     {return n_x;}
  double Get_x_cen()   {return x_cen;}
  double Get_t_begin() {return t_begin; }
  double Get_dx() {return dx; }
  double Get_Density(int i) { return rho[i]; }
  double Get_mu_e(int i) {return mu_e[i]; }

  double Get_vel(int i) {return vel[i]; }

  // add energy to the zone during a Compton scatter
  void add_edep(int i, double energy) {
    edep[i]= edep[i] + energy;
    return;
  }

  double Get_edep(int i) {
    return edep[i];
  }
};

#endif
