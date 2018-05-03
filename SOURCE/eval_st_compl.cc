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

  if(argc != 4){
    cerr << "Usage: NT comp file\n";
    return 1;
  }
  
  int NT = atoi(argv[1]);
  int comp = atoi(argv[2]);
  
  ifstream data_file(argv[3]);
  if(!data_file || NT <= 0 || (comp != 0 && comp != 1)) {
    cerr << "Wrong parameters.\n";
    return 1;
  }
  
  vector<double>* data = new vector<double>[NT];
  double tmp_1, tmp_2;
	
  int flag = 1;
  while(flag){
    for(int t = 0; t < NT; t++){
      data_file >> tmp_1;
      data_file >> tmp_2;
      if(data_file.eof()){ flag = 0; break; }
      if(comp == 0) data[t].push_back(tmp_1);
      if(comp == 1) data[t].push_back(tmp_2);
    }
  }
  
  int nconf = data[0].size();
  
  double ave, err;
  for(int t = 0; t < NT; t++){
    stat(data[t], ave, err);
    cout << t << " " <<  ave << " " << err << "\n";
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

