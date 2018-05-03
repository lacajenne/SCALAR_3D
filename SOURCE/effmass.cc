#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

void stat( const vector<double>& data, double& ave, double& error );

int main(int argc, char *argv[]){

  if(argc != 3){
    cerr << "Usage: NT file\n";
    return 1;
  }
	
  int NT = atoi(argv[1]);
  
  ifstream data_file(argv[2]);
  if(!data_file || NT <= 0){
    cerr << "Wrong parameters.\n";
    return 1;
  }
  
  vector<double>* data = new vector<double>[NT];
  vector<double>* cdata = new vector<double>[NT];
  double tmp;
	
  int flag = 1;
  while(flag){
    for(int t = 0; t < NT; t++){
      data_file >> tmp;
      if(data_file.eof()){ flag = 0; break; }
      data[t].push_back(tmp);
      cdata[t].push_back(0.0);
    }
  }
  
  int nconf = data[0].size();
  
  double s[NT];
  for(int t=0; t<NT; t++) {
    s[t] = 0.0;
    for(int nc=0; nc<nconf; nc++) {
      s[t] += 1.0/nconf * data[t].at(nc);
    }
  }
  
  for(int t=0; t<NT; t++) {
    for(int nc=0; nc<nconf; nc++) {
      for(int t0=0; t0<NT/2; t0++) {
	cdata[t].at(nc) += 2.0/NT * (data[t0].at(nc) * data[(t0+t)%NT].at(nc) - s[t0] * s[(t0+t)%NT]);
      }
    }
  }
  
  for(int nc=0; nc<nconf; nc++) {
    for(int t=0; t<NT/2; t++) {
      cdata[t].at(nc) = 0.5 * (cdata[t].at(nc) + cdata[NT-1-t].at(nc)); 
    }
  }
  
  double C[NT/2], E[NT/2];
  for(int t=0; t<NT/2; t++) stat(cdata[t], C[t], E[t]);
  
  for(int t=0; t<NT/2-1; t++) {
    double m = log(C[t]/C[t+1]);
    double err = (E[t]/C[t])*(E[t]/C[t]) + (E[t+1]/C[t+1])*(E[t+1]/C[t+1]);
    err = sqrt(err);
    cout << t << "    " << m << "    " << err << endl;
  }

  
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

