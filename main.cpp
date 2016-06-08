// Please see license.txt for licensing and copyright information //
// // Author: Paul Zimmerman, University of Michigan //
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "gstring.h"


using namespace std;


int main(int argc, char* argv[]){
  string inpfile;
  string xyzfile;
  string nprocs;
  switch (argc){
  case 1:
    inpfile="inpfileq";
    xyzfile="initial.xyz";
    nprocs="1";
    break;
  case 2:
    inpfile="inpfileq";
    xyzfile=argv[1];
    nprocs="1";
    break;
  case 3:
    inpfile="inpfileq";
    xyzfile=argv[1];
    nprocs=argv[2];
    break;
  default:
    cout << "Invalid command line options." << endl;
    return -1;
  }

#if 0
  double* A = new double[4]; 
  A[0] = A[3] = 2.;
  A[1] = A[2] = 0.;
  double* B = new double[4];
  B[0] = B[3] = 0.;
  B[1] = B[2] = 0.;

  mat_times_mat(B,A,A,2);
 
  printf(" done testing dgemm \n");
  exit(1);
#endif

  clock_t startTime = clock();

  cout << " MAIN: nprocs: " << nprocs << endl;
  int nnprocs = atoi(nprocs.c_str());
  int name = atoi(xyzfile.c_str());
//  cout << " MAIN: inpfile: " << inpfile << " xyzfile: " << xyzfile << endl;
  GString gstr;
//  gstr.init(inpfile, xyzfile);
  gstr.init(inpfile, name, nnprocs);
  gstr.String_Method_Optimization();


  clock_t endTime = clock();
  clock_t clockTicksTaken = endTime - startTime;
  double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;
  printf("\n\n\n total time: %8.5f \n",timeInSeconds);

  return 0;
}
