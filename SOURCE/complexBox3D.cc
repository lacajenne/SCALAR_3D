
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
double B2, B4, H, f0;

// lattice extent variables
int L, VOL;

// momentum units 
double ps0;

// random number stuff
long   int idum;
void   init_rg(long int); // random number initialization
double uniform_rg();      // uniform random number in [0,1]
double gaussian_rg();     // gaussian random numbers with zero ave. and unit variance

// HMC functions
double action( double** field );
void   eval_force( double** f_array, double** field );

// lattice geometry
int  coord_to_index( int t, int x, int y);  // find site index from coordinates
void init_lattice(); // assign the table of nearest neighbours

// calculation of some observables
double ave_field(double **field, int comp);
double kinetic_condensate(double **field);

// IO utility functions
void SaveLattice(double **field, char *filename);
void LoadLattice(double **field, char *filename);
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
  if ( argc != 12 ){
    cerr << "usage : <L> <B4> <B2> <H> <f0> <nconf> <nequil> <dt> <md_steps> <run_status(new:0, restart:1)> <eval S(t)(0: NO)>\n";
    exit(-1);
  }

  L          = atoi(argv[1]);
  B4         = atof(argv[2]);
  B2         = atof(argv[3]);
  H          = atof(argv[4]);
  f0         = atof(argv[5]);
  nconf      = atoi(argv[6]);
  nequil     = atoi(argv[7]);
  dt         = atof(argv[8]);
  md_steps   = atoi(argv[9]);
  run_status = atoi(argv[10]);
  evals      = atoi(argv[11]);

  // evaluate the S(t) field?
  eval_sfield = (evals != 0);

  // evaluate lattice volume
  VOL = L * L * L;

  // unit of lattice momentum
  ps0 = 2.0 * M_PI / ((double)L);

  // display simulation parameters on screen
  cout << "\nRUN PARAMETERS:\n\n";
  cout << "VOL    : " << L <<  " x " << L << " x " << L << "\n";
  cout << "B4     : " << B4 << "\n";
  cout << "B2     : " << B2 << "\n";
  cout << "H      : " << H << "\n";
  cout << "f0     : " << f0 << "\n";
  cout << "nconf  : " << nconf << "\n";
  cout << "nequil : " << nequil << "\n";
  cout << "md dt  : " << dt << "\n";
  cout << "md NS  : " << md_steps << "\n";
  cout << "eval S : " << evals << "\n";

  // set name of file with S(t)
  if(eval_sfield) {
    sprintf(file_name_corr, "CB3D_C_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
    corr_file = fopen(file_name_corr, "w");
    sprintf(file_name_corr_t, "CB3D_S_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
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
  double **field;
  field  = new double*[VOL];
  
  double **nfield;
  nfield = new double*[VOL];
  
  double **mom;
  mom    = new double*[VOL];
  
  double **momt;
  momt   = new double*[VOL];
  
  double **force;
  force  = new double*[VOL];
  
  for(int a=0; a<VOL; a++) {
    field[a]  = new double[2];
    nfield[a] = new double[2];
    mom[a]    = new double[2];
    momt[a]   = new double[2];
    force[a]  = new double[2];
  }
  
  // read field configuration from file or start anew
  if(run_status==NEW_RUN){
    // initialize field
    for(int site=0; site<VOL; site++){
      double r0 = 0.05;
      double r = 1.0 + (2.0*uniform_rg() - 1.0) * r0;
      double s = 1.0 + (2.0*uniform_rg() - 1.0) * r0;
      field[site][0] = f0 * r;
      field[site][1] = f0 * s;
    }
  } else {
    // load last saved field configuration
    sprintf(latfile, "CB3D_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
    LoadLattice(field, latfile);
  }
  
  // initialize trial field and action
  for(int site=0; site<VOL; site++) {
    nfield[site][0] = field[site][0];
    nfield[site][1] = field[site][1];
  }
  S = action(field);
  
  // loop over configurations
  int acc_count = 0;
  while( acc_count < nequil + nconf ){
    
    acc_count += 1;
    
    for(int a=0; a<VOL; a++){
      nfield[a][0] = field[a][0];
      nfield[a][1] = field[a][1];
      mom[a][0] = gaussian_rg();
      mom[a][1] = gaussian_rg();
    }
    
    E = S;
    for(int a=0; a<VOL; a++) {
      E += mom[a][0] * mom[a][0];
      E += mom[a][1] * mom[a][1];
    }
    
    eval_force(force, nfield);
    
    for(int a=0;a<VOL;a++) {
      momt[a][0] = mom[a][0] + 0.5 * dt * force[a][0];
      momt[a][1] = mom[a][1] + 0.5 * dt * force[a][1];
    }
    
    for(int step=0; step<md_steps; step++){
      for(int a=0; a<VOL; a++) {
	nfield[a][0] += dt * momt[a][0];
	nfield[a][1] += dt * momt[a][1];
      }
      
      eval_force(force, nfield);
      
      for(int a=0;a<VOL;a++) {
	momt[a][0] += dt * force[a][0];
	momt[a][1] += dt * force[a][1];
      }
    }
    
    for(int a=0; a<VOL; a++) {
      mom[a][0] = momt[a][0] - 0.5 * dt * force[a][0];
      mom[a][1] = momt[a][1] - 0.5 * dt * force[a][1];
    }
    
    Sn = action(nfield);
    En = Sn;
    for(int a=0; a<VOL; a++) {
      En += 0.5 * mom[a][0] * mom[a][0];
      En += 0.5 * mom[a][1] * mom[a][1];
    }
    
    dE = En - E;
    status = REJECT;
    if( dE < 0.0 || uniform_rg() <= exp(-dE) ) status = ACCEPT;
    
    if(acc_count < nequil){
      sprintf(logbuf, "%.16f", 1.0*status);
      sprintf(logfile, "CB3D_A_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
      Log(logbuf, logfile);
    }
    
    // it the configuration is accepted, update the field
    if(status == ACCEPT){
      S = Sn;
      for(int a=0;a<VOL;a++) {
	field[a][0] = nfield[a][0];
	field[a][1] = nfield[a][1];
      }
    }
    
    // save observables only after thermalization
    //  (and only every nsave configurations, to reduce correlation)
    if(acc_count > nequil){
      if(acc_count%nsave == 0){

    	// average field 
    	sprintf(logbuf, "%.16f", ave_field(field, 0));
	sprintf(logfile, "CB3D_F0_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
    	Log(logbuf, logfile);
    	sprintf(logbuf, "%.16f", ave_field(field, 1));
	sprintf(logfile, "CB3D_F1_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
    	Log(logbuf, logfile);
	
	// kinetic condensate
    	sprintf(logbuf, "%.16f", kinetic_condensate(field));
	sprintf(logfile, "CB3D_K_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
    	Log(logbuf, logfile);
	
	if(eval_sfield) {
	  // S(t) field, averaged over 2 dimensions
	  // used to estimate C(t)
	  double sdata[2][L];
          for(int t=0; t<L; t++) {
            sdata[0][t] = 0.0;
	    sdata[1][t] = 0.0;
            for(int x=0; x<L; x++) {
              for(int y=0; y<L; y++) {
                int ind = coord_to_index(t, x, y);
                sdata[0][t] += field[ind][0];
		sdata[1][t] += field[ind][1];
              }
            }
            sdata[0][t] /= ((double)L);
	    sdata[1][t] /= ((double)L);
            fprintf(corr_file_t, "%f\t", sdata[0][t]);
	    fprintf(corr_file_t, "%f\t", sdata[1][t]);
          }
          fprintf(corr_file_t, "\n");
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
  sprintf(latfile, "CB3D_L%d_B4_%.3fB2_%.3f_H%.3f.dat ", L, B4, B2, H);
  SaveLattice(field, latfile);

  delete[] field;
  delete[] nfield;
  delete[] mom;
  delete[] momt;
  delete[] force;

}


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
  return t + L * x + L*L * y;
}


void init_lattice() {

  nnp = new int*[VOL];
  nnn = new int*[VOL];

  for(int a=0;a<VOL;a++){

    nnp[a] = new int[DIM]; // one n.n for each pos. dir. for each lattice point (a)
    nnn[a] = new int[DIM]; // one n.n for each neg. dir. for each lattice point (a)

    // extract the coordinates (t, x, y, z) from the site index
    int t, x, y, b;
    t = a % L;
    b = a / L;
    x = b % L;
    b = b / L;
    y = b % L;

    // identify nearest neighbours implementing periodic b.c.

    nnp[a][DIR_T] = t == (L-1) ? coord_to_index(0, x, y)   : coord_to_index(t+1, x, y);
    nnn[a][DIR_T] = t == 0     ? coord_to_index(L-1, x, y) : coord_to_index(t-1, x, y);

    nnp[a][DIR_X] = x == (L-1) ? coord_to_index(t, 0, y)   : coord_to_index(t, x+1, y);
    nnn[a][DIR_X] = x == 0     ? coord_to_index(t, L-1, y) : coord_to_index(t, x-1, y);

    nnp[a][DIR_Y] = y == (L-1) ? coord_to_index(t, x, 0)   : coord_to_index(t, x, y+1);
    nnn[a][DIR_Y] = y == 0     ? coord_to_index(t, x, L-1) : coord_to_index(t, x, y-1);

  }
}

inline double action( double** field )
{
  // this method evaluates the action on a field configuration
  
  double S = 0.0;
  double f0, f1, chim, tb0, tb1;
  
  for(int a=0; a<VOL; a++){
    
    f0 = field[a][0];
    f1 = field[a][1];
    
    chim = f0*f0 + f1*f1;
    
    S += 2.0 * B2 * DIM * chim;
    for(int mu=0; mu<DIM; mu++) {
      S -= 2.0 * B2 * f0 * field[nnp[a][mu]][0];
      S -= 2.0 * B2 * f1 * field[nnp[a][mu]][1];
    }
    
    S += 0.5 * H * chim * chim;
    
    if(B4 != 0.0) {
      tb0 = 0.0;
      tb1 = 0.0;
      for(int mu=0; mu<DIM; mu++) {
	tb0 += field[nnp[a][mu]][0] + field[nnn[a][mu]][0] - 2.0*f0;
	tb1 += field[nnp[a][mu]][1] + field[nnn[a][mu]][1] - 2.0*f1;
      }
      S += B4/2.0 * (tb0*tb0 + tb1*tb1);
    }
    
  }
  
  return S;
}


inline void eval_force( double** f_array, double** field )
{
  // this method evaluates the HMC force, i.e.
  //  - dS/dphi(x)

  double f0, f1;
  double tb0, tb1;

  for(int a=0; a<VOL; a++){

    f0 = field[a][0];
    f1 = field[a][1];
    
    f_array[a][0] = -4.0*B2*DIM*f0 - 2.0 * H * (f0*f0 + f1*f1) * f0; 
    f_array[a][1] = -4.0*B2*DIM*f1 - 2.0 * H * (f0*f0 + f1*f1) * f1; 
    
    for(int mu=0;mu<DIM;mu++){
      f_array[a][0] += 2.0 * B2 * (field[nnp[a][mu]][0] + field[nnn[a][mu]][0]);
      f_array[a][1] += 2.0 * B2 * (field[nnp[a][mu]][1] + field[nnn[a][mu]][1]);
    }
    
    if(B4 != 0.0) {
      tb0 = 4.0 * DIM * DIM * f0;
      tb1 = 4.0 * DIM * DIM * f1;
      for(int mu=0;mu<DIM;mu++){
	int p_mu = nnp[a][mu];
	int n_mu = nnn[a][mu];
	tb0 -= 4.0 * DIM * (field[p_mu][0] + field[n_mu][0]);
	tb1 -= 4.0 * DIM * (field[p_mu][1] + field[n_mu][1]);
	for(int nu=0;nu<DIM;nu++){
	  int np_mu_nu = nnp[n_mu][nu];
	  int pp_mu_nu = nnp[p_mu][nu];
	  int pn_mu_nu = nnn[p_mu][nu];
	  int nn_mu_nu = nnn[n_mu][nu];
	  tb0 += field[pp_mu_nu][0] + field[pn_mu_nu][0] + field[np_mu_nu][0] + field[nn_mu_nu][0];
	  tb1 += field[pp_mu_nu][1] + field[pn_mu_nu][1] + field[np_mu_nu][1] + field[nn_mu_nu][1];
	}
      }
      f_array[a][0] -= B4 * tb0;
      f_array[a][1] -= B4 * tb1;
    }
    
  }
}

inline double ave_field(double **field, int comp)
{
  double ret = 0.0;
  for(int site=0; site<VOL; site++) ret += field[site][comp];
  ret /= VOL;
  return ret;
}

inline double kinetic_condensate(double **field)
{
  double ret = 0.0;
  for(int site=0; site<VOL; site++) {
    for(int mu=0; mu<DIM; mu++) {
      int np = nnp[site][mu];
      double d0 = field[np][0] - field[site][0];
      double d1 = field[np][1] - field[site][1];
      ret += d0*d0 + d1*d1;
    }
  }
  ret /= VOL;
  return ret;
}

inline void SaveLattice(double **field, char *filename)
{
  FILE *file;
  file = fopen(filename, "w");
  for(int i=0;i<VOL;i++){
    fprintf(file, "%.20f\t%.20f\n", field[i][0], field[i][1]);
  }
  fflush(file);
  fclose(file);
}


inline void LoadLattice(double **field, char *filename)
{
  FILE *file;
  file = fopen(filename, "r");
  for(int i=0;i<VOL;i++){
    int p = 0;
    p = fscanf(file, "%lf\t%lf\n", &field[i][0], &field[i][1]);
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
