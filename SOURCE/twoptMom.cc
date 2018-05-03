#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

// (t, x, y)
#define CHAN0R 0 // (0, 0, 0) R
#define CHAN1R 1 // (0, k, 0) R
#define CHAN1I 2 // (0, k, 0) I
#define CHAN2R 3 // (0, k, k) R
#define CHAN2I 4 // (0, k, k) I
#define CHAN3R 5 // (p, k, k) R
#define CHAN3I 6 // (p, k, k) I

using namespace std;

const int NCHAN = 7;
const int NMOM = 4;

double ps0;
double pt0;

void stat( const vector<double>& data, double& ave, double& error );
double latmom(int nt, int nx, int ny);

int main(int argc, char *argv[]){

  if(argc != 4){
    cerr << "Usage: L T file\n";
    return 1;
  }
  
  int L = atoi(argv[1]);
  int T = atoi(argv[2]);

  ifstream data_file(argv[3]);
  if(!data_file){
    cerr << "Wrong parameters.\n";
    return 1;
  }
  
  ps0 = 2.0 * M_PI / ((double)L);
  pt0 = 2.0 * M_PI / ((double)T);

  vector<double>* data = new vector<double>[NCHAN];
  vector<double>* pdata = new vector<double>[NMOM];
  double avef[NCHAN];

  double tmp;
  int flag = 1;
  while(flag){

    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN0R].push_back(tmp);
    
    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN1R].push_back(tmp);
    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN1I].push_back(tmp);
    
    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN2R].push_back(tmp);
    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN2I].push_back(tmp);
    
    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN3R].push_back(tmp);
    data_file >> tmp;
    if(data_file.eof()){ flag = 0; break; }
    data[CHAN3I].push_back(tmp);
  }
  
  int nconf = data[0].size();
  
  for(int nc=0; nc<nconf; nc++) {
    for(int k=0; k<NMOM; k++) {
      pdata[k].push_back(0.0);
    }
  }
  
  for(int k=0; k<NCHAN; k++) {
    avef[k] = 0.0;
    for(int nc=0; nc<nconf; nc++) {
      avef[k] += data[k].at(nc);
    }
    avef[k] /= ((double)nconf);
  }


  for(int nc=0; nc<nconf; nc++) {
    
    double tmpR, tmpI;
    
    pdata[0].at(nc) = data[0].at(nc)*data[0].at(nc) - avef[0]*avef[0];

    tmpR = data[1].at(nc) * data[1].at(nc) - avef[1]*avef[1];
    tmpI = data[2].at(nc) * data[2].at(nc) - avef[2]*avef[2];
    pdata[1].at(nc) = tmpR + tmpI;

    tmpR = data[3].at(nc) * data[3].at(nc) - avef[3]*avef[3];
    tmpI = data[4].at(nc) * data[4].at(nc) - avef[4]*avef[4];
    pdata[2].at(nc) = tmpR + tmpI;

    tmpR = data[5].at(nc) * data[5].at(nc) - avef[5]*avef[5];
    tmpI = data[6].at(nc) * data[6].at(nc) - avef[6]*avef[6];
    pdata[3].at(nc) = tmpR + tmpI;

  }


  double mom;
  double ave, err;
  
  mom = latmom(0, 0, 0);
  stat(pdata[0], ave, err);
  cout << mom << " " << ave << " " << err << endl;

  mom = latmom(0, 1, 0);
  stat(pdata[1], ave, err);
  cout << mom << " " << ave << " " << err << endl;

  mom = latmom(0, 1, 1);
  stat(pdata[2], ave, err);
  cout << mom << " " << ave << " " << err << endl;

  mom = latmom(1, 1, 1);
  stat(pdata[3], ave, err);
  cout << mom << " " << ave << " " << err << endl;

  data_file.close();
  
  cerr << "ave on " << nconf << "\n";
  
  return 0;
}


void stat( const vector<double>& data, double& ave, double& error )
{
  double* bin_ave;
  unsigned long int ndata = data.size();
  bin_ave = new double[ndata];
  
  ave = 0.0;
  for(unsigned long int n = 0; n < ndata; n++) ave += data[n];
  ave /= ndata;
	
  error = 0.0;
  for(unsigned long int bin_size = 1; bin_size <= ndata; bin_size ++){
    int nbins = ndata/bin_size;
    for(unsigned long int i = 0; i < nbins; i++){
      bin_ave[i] = 0.0;
      for(unsigned long int k = 0; k < bin_size; k++)
	bin_ave[i] += data[k+i*bin_size];
      bin_ave[i] /= bin_size;
    }
    double tmp = 0.0;
    for(unsigned long int i = 0; i < nbins; i++)
      tmp += ( bin_ave[i] - ave ) * ( bin_ave[i] - ave );
    tmp = sqrt(tmp)/nbins;
    if(tmp > error) error = tmp;
  }
  
  delete[] bin_ave;
}

double latmom(int nt, int nx, int ny)
{
  double res = 0.0;
  double pt = pt0 * nt;
  double px = ps0 * nx;
  double py = ps0 * ny;
  res += sin(pt/2.0) * sin(pt/2.0);
  res += sin(px/2.0) * sin(px/2.0);
  res += sin(py/2.0) * sin(py/2.0);
  res = 2.0 * sqrt(res);
  return res;
}
