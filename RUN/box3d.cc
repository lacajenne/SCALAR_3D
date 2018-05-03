
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;

const int DIM   = 3;  // number of lattice dimensions
const int DIR_T = 0;  // indices for directions
const int DIR_X = 1;
const int DIR_Y = 2;

const int NEW_RUN = 0;
const int REJECT  = 0;  // Metropolis
const int ACCEPT  = 1;

// tables of nearest neighbour sites
// for all directions(sites have a single index: 0..VOL-1),
// pos.(p) and neg.(n)
int **nnp;
int **nnn;

// parameters of the model
// B2: coefficient of the standard kinetic term, usually 1
// B4: coefficient of the box^2 term:  1/2 B4 (box phi)^2
// G, H: 1/2 G phi^2 + 1/4 H phi^4
// J: external current
// f0: initial value of the field 4
double B2, B4, G, H, J, f0;

// lattice extent variables
int T, L, VOL;

// momentum units 
// 2 pi/ L, 2 pi / T 
double ps0, pt0;

// random number stuff
long   int idum;
void   init_rg(long int); // random number initialization
double uniform_rg();      // uniform random number in [0,1]
double gaussian_rg();     // gaussian random numbers with zero ave. and unit variance

// HMC functions
double action( double* field );
void   eval_force( double* f_array, double* field );

// lattice geometry
int  coord_to_index( int t, int x, int y);  // find site index from coordinates
void init_lattice(); // assign the table of nearest neighbours

// calculation of some observables
double ave_field(double *field);
double ave_field_square(double *field);
double kinetic_condensate(double *field);

// IO utility functions
void SaveLattice(double *field, char *filename);
void LoadLattice(double *field, char *filename);
void Log(char *msg, char *filename);


int main( int argc, char *argv[])
{
  int run_status; // restart? 0 for new run, 1 for restart
  int evals;      // 0: do not evaluate the S(t) field
  bool eval_sfield; 

  int nconf;  //  number of field configurations to be generated
  int nequil; //  number of thermalization field configurations (not used for meas.)

  char logfile[128]; // name of file with results for observables
  char latfile[128]; // name of file with lattice field configuration
  char logbuf[512];  // buffer for results to be piped into logfile

  // correlator data files ( C(t) and momentum two pt)
  FILE *corr_file;
  char file_name_corr[128];
  FILE *corr_file_t;
  char file_name_corr_t[128];

  // HMC variables
  double E, En, dE, S, Sn;
  double dt;
  int md_steps;
  int nsave = 15;
  int status;

  // take values of parameters from the command line
  if ( argc != 15 ){
    cerr << "usage : <L> <T> <B2> <B4> <G> <H> <J> <f0> <nconf> <nequil> <dt> <md_steps> <run_status(new:0, restart:1)> <eval S(t)(0: NO)>\n";
    exit(-1);
  }

  L          = atoi(argv[1]);
  T          = atoi(argv[2]);
  B2         = atof(argv[3]);
  B4         = atof(argv[4]);
  G          = atof(argv[5]);
  H          = atof(argv[6]);
  J          = atof(argv[7]);
  f0         = atof(argv[8]);
  nconf      = atoi(argv[9]);
  nequil     = atoi(argv[10]);
  dt         = atof(argv[11]);
  md_steps   = atoi(argv[12]);
  run_status = atoi(argv[13]);
  evals      = atoi(argv[14]);

  // evaluate the S(t) field?
  eval_sfield = (evals != 0);

  // evaluate lattice volume
  VOL = L * L * T;

  // units of lattice momentum
  ps0 = 2.0 * M_PI / ((double)L);
  pt0 = 2.0 * M_PI / ((double)T);

  // display simulation parameters on screen
  cout << "\nRUN PARAMETERS:\n\n";
  cout << "VOL    : " << L <<  " x " << L << " x " << T << "\n";
  cout << "B2     : " << B2 << "\n";  
  cout << "B4     : " << B4 << "\n";
  cout << "G      : " << G << "\n";
  cout << "H      : " << H << "\n";
  cout << "J      : " << J << "\n";
  cout << "f0     : " << f0 << "\n";
  cout << "nconf  : " << nconf << "\n";
  cout << "nequil : " << nequil << "\n";
  cout << "md dt  : " << dt << "\n";
  cout << "md NS  : " << md_steps << "\n";
  cout << "eval S : " << evals << "\n";

  // set name of file with S(t)
  if(eval_sfield) {
    sprintf(file_name_corr, "B3D_C_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
    corr_file = fopen(file_name_corr, "w");
    sprintf(file_name_corr_t, "B3D_S_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
    corr_file_t = fopen(file_name_corr_t, "w");
  }

  // new run?
  if(run_status==NEW_RUN){
    cout << "\nNEW RUN ...\n\n";
  } else {
    cout << "\nRESTART ...\n\n";
  }

  // correct for the decorrelation (see later)
  nconf *= nsave;

  // initialize the lattice geometry and the random numbers
  init_lattice();
  init_rg(-time(NULL));

  // the fields
  double *field  = new double[VOL];
  double *nfield = new double[VOL];
  double *mom    = new double[VOL];
  double *momt   = new double[VOL];
  double *force  = new double[VOL];
  double *corr   = new double[L];

  // read field configuration from file or start anew
  if(run_status==NEW_RUN){
    // initialize field
    for(int site=0; site<VOL; site++){
      double r0 = 0.05;
      double r = 1.0 + (2.0*uniform_rg()-1.0)*r0;
      field[site] = f0 * r;
    }
  } else {
    // load last saved field configuration
    sprintf(latfile, "B3D_L_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
    LoadLattice(field, latfile);
  }

  // initialize trial field and action
  for(int site=0; site<VOL; site++) nfield[site] = field[site];
  S = action(field);

  // loop over configurations
  int acc_count = 0;
  while( acc_count < nequil + nconf ){

    acc_count += 1;

    for(int a=0; a<VOL; a++){
      nfield[a] = field[a];
      mom[a]    = gaussian_rg();
    }

    E = S;
    for(int a=0;a<VOL;a++) E += 0.5 * mom[a] * mom[a];
    eval_force(force, nfield);
    for(int a=0;a<VOL;a++) momt[a] = mom[a] + 0.5 * dt * force[a];

    for(int step=0;step<md_steps;step++){
      for(int a=0;a<VOL;a++) nfield[a] += dt * momt[a];
      eval_force(force, nfield);
      for(int a=0;a<VOL;a++) momt[a]   += dt * force[a];
    }

    for(int a=0;a<VOL;a++) mom[a] = momt[a] - 0.5 * dt * force[a];

    Sn = action( nfield );
    En = Sn;
    for(int a=0;a<VOL;a++) En += 0.5 * mom[a] * mom[a];

    dE = En - E;
    status = REJECT;
    if( dE < 0.0 || uniform_rg() <= exp(-dE) ) status = ACCEPT;

    if(acc_count < nequil){
      sprintf(logbuf, "%.16f", exp(-dE));
      sprintf(logfile, "B3D_EH_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
      Log(logbuf, logfile);
      
      sprintf(logbuf, "%.16f", 1.0*status);
      sprintf(logfile, "B3D_A_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
      Log(logbuf, logfile);
    }

    // it the configuration is accepted, update the field
    if(status == ACCEPT){
      S = Sn;
      for(int a=0;a<VOL;a++) field[a] = nfield[a];
    }

    // save observables only after thermalization
    //  (and only every nsave configurations, to reduce correlation)
    if(acc_count > nequil){
      if(acc_count%nsave == 0){

    	// average field
    	sprintf(logbuf, "%.16f", ave_field(field));
	sprintf(logfile, "B3D_F_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
    	Log(logbuf, logfile);
	
	// kinetic condensate
    	sprintf(logbuf, "%.16f", kinetic_condensate(field));
	sprintf(logfile, "B3D_K_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
    	Log(logbuf, logfile);

	if(eval_sfield) {
	  
	  // S(t) field, averaged over 2 dimensions
	  // used to estimate C(t)
	  double sdata[T];
          for(int t=0; t<T; t++) {
            sdata[t] = 0.0;
            for(int x=0; x<L; x++) {
              for(int y=0; y<L; y++) {
                int ind = coord_to_index(t, x, y);
                sdata[t] += field[ind];
              }
            }
            sdata[t] /= (double)L;
            fprintf(corr_file_t, "%f\t", sdata[t]);
          }
          fprintf(corr_file_t, "\n");
	  
	  // fourier transform of the field, for several values
	  //  of the lattice momentum

	  double sp, phC, phS;
	  double ftR_0 = 0.0;
	  double ftR_1 = 0.0;
	  double ftI_1 = 0.0;
	  double ftR_2 = 0.0;
	  double ftI_2 = 0.0;
	  double ftR_3 = 0.0;
	  double ftI_3 = 0.0;

	  double KC = sqrt((double)VOL);
	  
	  for(int a=0; a<VOL; a++) {
	    
	    int t, x, y, b;
	    t = a % T;
	    b = a / T;
	    x = b % L;
	    b = b / L;
	    y = b % L;
	    
	    ftR_0 += field[a] / KC;
	    
	    sp = x * ps0;
	    phC = cos(sp);
	    phS = sin(sp);
	    ftR_1 += phC * field[a] / KC;
	    ftI_1 += phS * field[a] / KC;

	    sp = (x + y) * ps0;
	    phC = cos(sp);
	    phS = sin(sp);
	    ftR_2 += phC * field[a] / KC;
	    ftI_2 += phS * field[a] / KC;

	    sp = (x + y) * ps0 + t * pt0;
	    phC = cos(sp);
	    phS = sin(sp);
	    ftR_3 += phC * field[a] / KC;
	    ftI_3 += phS * field[a] / KC;
	    
	  }
	  
	  fprintf(corr_file, "%f\t", ftR_0);
	  fprintf(corr_file, "%f\t", ftR_1);
	  fprintf(corr_file, "%f\t", ftI_1);
	  fprintf(corr_file, "%f\t", ftR_2);
	  fprintf(corr_file, "%f\t", ftI_2);
	  fprintf(corr_file, "%f\t", ftR_3);
	  fprintf(corr_file, "%f\t", ftI_3);
	  
	  fprintf(corr_file, "\n");

	}
      }
    }
  }
 
  
  // close the correlator files properly
  if(eval_sfield){
    fflush(corr_file);
    fclose(corr_file);
    fflush(corr_file_t);
    fclose(corr_file_t);
  }

  // save the last field configuration on file
  sprintf(latfile, "B3D_L_L%d_T%d_B2_%.3f_B4_%.3f_G%.3f_H%.3f_J%.3f.dat ", L, T, B2, B4, G, H, J);
  SaveLattice(field, latfile);

  delete[] field;
  delete[] nfield;
  delete[] mom;
  delete[] momt;
  delete[] force;
  delete[] corr;

} // end of main


void init_rg(long int id)
{
  idum = id;
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


double uniform_rg()
{
#ifdef GCC_RAND
  return rand()/((double)RAND_MAX);
#else

  //tratto da http://www.library.cornell.edu/nr/cbookcpdf.html

  int j;
  long int k;
  static long int idum2 = 123456789;
  static long int iy = 0;
  static long iv[NTAB];
  double temp;

  if(idum <= 0){
    if(-idum < 1)
      idum = 1;
    else
      idum = -idum;
    idum2 = idum;
    for(j = NTAB+7; j >= 0; j--){
      k = idum/IQ1;
      idum = IA1*(idum-k*IQ1)-k*IR1;
      if(idum < 0) idum += IM1;
      if(j < NTAB) iv[j] = idum;
    }
    iy = iv[0];
  }
  k = idum/IQ1;
  idum = IA1*(idum-k*IQ1)-k*IR1;
  if(idum < 0) idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  if(idum2 < 0) idum2 += IM2;
  j = iy/NDIV;
  iy = iv[j] - idum2;
  iv[j] = idum;
  if(iy < 1) iy += IMM1;
  if((temp=AM*iy) > RNMX) return RNMX;
  return temp;
#endif
}

double gaussian_rg()
{
  // Genera x con distribuzione exp(-.5*x^2)dx in (-\infty,+\infty).

  double v1, v2, zsq, R;
  static int flag = 0;
  static double gset;
  if(flag == 0){
    flag = 1;
    do{
      v1 = 2.0*uniform_rg()-1.0;
      v2 = 2.0*uniform_rg()-1.0;
      zsq = v1*v1 + v2*v2;
    } while(zsq > 1 || zsq == 0);
    R = sqrt(-2.*log(zsq)/zsq);
    gset = v1 * R;
    return v2 * R;
  }

  flag = 0;
  return gset;
}


inline int coord_to_index( int t, int x, int y)
{
  // return the site index from the coordinates (t, x, y)
  return t + T * x + T*L * y;
}


void init_lattice() {

  nnp = new int*[VOL];
  nnn = new int*[VOL];

  for(int a=0;a<VOL;a++){

    nnp[a] = new int[DIM]; // one n.n for each pos. dir. for each lattice point (a)
    nnn[a] = new int[DIM]; // one n.n for each neg. dir. for each lattice point (a)

    // extract the coordinates (t, x, y, z) from the site index
    int t, x, y, b;
    t = a % T;
    b = a / T;
    x = b % L;
    b = b / L;
    y = b % L;

    // identify nearest neighbours implementing periodic b.c.

    nnp[a][DIR_T] = t == (T-1) ? coord_to_index(0, x, y)   : coord_to_index(t+1, x, y);
    nnn[a][DIR_T] = t == 0     ? coord_to_index(T-1, x, y) : coord_to_index(t-1, x, y);

    nnp[a][DIR_X] = x == (L-1) ? coord_to_index(t, 0, y)   : coord_to_index(t, x+1, y);
    nnn[a][DIR_X] = x == 0     ? coord_to_index(t, L-1, y) : coord_to_index(t, x-1, y);

    nnp[a][DIR_Y] = y == (L-1) ? coord_to_index(t, x, 0)   : coord_to_index(t, x, y+1);
    nnn[a][DIR_Y] = y == 0     ? coord_to_index(t, x, L-1) : coord_to_index(t, x, y-1);

  }
}

inline double action( double* field )
{
  // this method evaluates the action on a field configuration

  double S = 0.0;
  double f, tb;

  for(int a=0; a<VOL; a++){
    f = field[a];
    S += B2 * DIM * f * f;
    for(int mu=0; mu<DIM; mu++) S -= B2 * f * field[nnp[a][mu]];
    S += 0.5*G*f*f + 0.25*H*f*f*f*f + J*f;

    if(B4 != 0.0) {
      tb = 0.0;
      for(int mu=0; mu<DIM; mu++) tb += field[nnp[a][mu]] + field[nnn[a][mu]] - 2.0*f;
      S += B4/2.0 * tb * tb;
    }

  }
  return S;
}


inline void eval_force( double* f_array, double* field )
{
  // this method evaluates the HMC force, i.e.
  //  - dS/dphi(x)

  double f;
  double tb;

  for(int a=0; a<VOL; a++){

    f = field[a];
    f_array[a] = -2*B2*DIM*f - G*f - H*f*f*f - J;

    for(int mu=0;mu<DIM;mu++){
      f_array[a] += B2 * (field[nnp[a][mu]] + field[nnn[a][mu]]);
    }
    
    if(B4 != 0.0) {
      tb = 4.0 * DIM * DIM * f;
      for(int mu=0;mu<DIM;mu++){
	int p_mu = nnp[a][mu];
	int n_mu = nnn[a][mu];
	tb += -4.0 * DIM * (field[p_mu] + field[n_mu]);
	for(int nu=0;nu<DIM;nu++){
	  int np_mu_nu = nnp[n_mu][nu];
	  int pp_mu_nu = nnp[p_mu][nu];
	  int pn_mu_nu = nnn[p_mu][nu];
	  int nn_mu_nu = nnn[n_mu][nu];
	  tb += field[pp_mu_nu] + field[pn_mu_nu] + field[np_mu_nu] + field[nn_mu_nu];
	}
      }
      f_array[a] += -B4 * tb;
    }

  }
}

inline double ave_field(double *field)
{
  double ret = 0.0;
  for(int site=0; site<VOL; site++) ret += field[site];
  ret /= VOL;
  return ret;
}

inline double ave_field_square(double *field)
{
  double ret = 0.0;
  for(int site=0; site<VOL; site++) ret += field[site] * field[site];
  ret /= VOL;
  return ret;
}

inline double kinetic_condensate(double *field)
{
  double ret = 0.0;
  double fs = ave_field_square(field);
  for(int site=0; site<VOL; site++) {
    for(int mu=0; mu<DIM; mu++) {
      int np = nnp[site][mu];
      double d = field[np] - field[site];
      ret += d*d;
    }
  }
  ret /= VOL;
  ret += G * fs;
  return ret;
}

inline void SaveLattice(double *field, char *filename)
{
  FILE *file;
  file = fopen(filename, "w");
  for(int i=0;i<VOL;i++){
    fprintf(file, "%.40f\n", field[i]);
  }
  fflush(file);
  fclose(file);
}


inline void LoadLattice(double *field, char *filename)
{
  FILE *file;
  file = fopen(filename, "r");
  for(int i=0;i<VOL;i++){
    int p = 0;
    p = fscanf(file, "%lf\n", &field[i]);
  }
  fflush(file);
  fclose(file);
}


inline void Log(char *msg, char *filename){
  FILE *file;
  file = fopen(filename, "a");
  fprintf(file, "%s\n", msg);
  fflush(file);
  fclose(file);
}
