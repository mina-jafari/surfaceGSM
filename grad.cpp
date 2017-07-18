// Please see license.txt for licensing and copyright information //
// // Author: Paul Zimmerman, University of Michigan //
#include "grad.h"
//#include "mopac.h"

using namespace std;

void Gradient::update_knnr()
{
#if !USE_KNNR
  return;
#endif
  if (knnr1.npts>0 && knnr_inited)
    int newpts = knnr1.add_extra_points();
  else
  {
    printf("  loading KNNR files \n");
    int kfound = knnr1.begin(runNum,natomsg);
    if (kfound>0)
      res_t = 1000;
    knnr_inited = 1;
  }
  return;
}


//type 0: QM gradient only
//type 1: knnr if accurate
//type 2: knnr only
//type 3: knnr if exact

double Gradient::grads(double* coords, double* grad, double* Ut, int type)
{
  //printf(" grads(%i) k: %i \n",type,knn_k);

#if USE_KNNR
  if (type==2 && knnr_active==0)
  {
    printf("\n ERROR: grad type 2 (knnr only) while knnr_active==0 \n");
    exit(1);
  }
#endif

  if (knnr_active && knnr_inited==0 && type>0) update_knnr();

  if (knnr_active==3) type = 3; //force type 3

  wrote_grad = 0;

  double errknn = 10.;
#if USE_KNNR
  if (knnr_active && type>0)
  {
    //printf(" about to call grad_knnr() \n");
    errknn = knnr1.grad_knnr(coords,energy,grad,Ut,knn_k);
    energy *= 627.5;
    //printf(" E(kNN): %4.3f",energy);
    xyz_grad = 0;
    if (energy - V0 < -500. || energy - V0 > 500.)
    {
      printf(" kf");
      errknn = 99.;
    }
  }
#else
  energy = -99.;
  type = 0;
#endif

  //if (type==3) printf(" t3errknn: %5.4f",errknn);

  int do_exact = 0;
  if (type==0) do_exact = 1;
  else if (type==1 && errknn > KNNR_MAX_DIST) do_exact = 1;
  else if (type==2) do_exact = 0;
  else if (type==3 && errknn > 0.01) do_exact = 1;

  if (always_do_exact) do_exact = 1;
  //do_exact = 1; //debug
  if (do_exact)
  {
    printf(" eg"); fflush(stdout);
    int success = external_grad(coords,grad);
    xyz_grad = 1;
  }
  else
    printf(" kg"); 

#if 0
  if (xyz_grad)
  {
    printf(" Grad: \n");
    for (int i=0;i<natoms;i++)
      printf(" %s %12.10f %12.10f %12.10f \n",anames[i].c_str(),grad[3*i+0],grad[3*i+1],grad[3*i+2]);
  }
#endif

#if 0
  printf(" XYZ: \n");
  for (int i=0;i<natoms;i++)
    printf(" %s %4.3f %4.3f %4.3f \n",anames[i].c_str(),coords[3*i+0],coords[3*i+1],coords[3*i+2]);
#endif


  return energy;
}

int Gradient::external_grad(double* coords, double* grad)
{

#if QCHEM
//  printf(" gqc"); fflush(stdout);
  energy = qchem1.grads(coords,grad);
//  printf(" gqce"); fflush(stdout);
#elif USE_ASE
  //printf(" grad ase \n"); fflush(stdout);
  energy = ase1.grads(coords,grad);
  //printf(" done grad ase \n"); fflush(stdout);
#else
  char* pbsPath;
  pbsPath = getenv ("PBSTMPDIR");
  string pdir = "";
  if (pbsPath!=NULL)
  {
    string pstr(pbsPath);
    pdir = pstr + "/";
  }
  Mopac mp1; 
  mp1.alloc(natomsg);
  mp1.reset(natomsg,anumbers,anames,coords);
  energy = mp1.grads(pdir+"mxyzfile"+runends);
  for (int i=0;i<3*natomsg;i++)
    grad[i] = mp1.grad[i];
  mp1.freemem();
#endif
  gradcalls++;

  int success = 1;
  if (V0==0.0)  V0 = energy;
  if (energy-V0>1000. || energy-V0 < -1000.)  success = 0;

//  printf(" write_on: %i success: %i energy: %6.5f \n",write_on,success,energy);

#if WRITE_FILES
  if (write_on && success)
  {
    string nstr = StringTools::int2str(gradcalls+res_t,4,"0");
    string filename = "scratch/qcsave"+runName0+"."+nstr;
    write_xyz_grad(coords,grad,filename);
  }
#endif

//printf(" grad ending early \n");
//exit(1);

  return success;
}


void Gradient::write_xyz_grad(double* coords, double* grad, string filename)
{
  wrote_grad = 1;

  ofstream xyzfile;
  string xyzfile_string = filename+".xyz";
  xyzfile.open(xyzfile_string.c_str());
  xyzfile.setf(ios::fixed);
  xyzfile.setf(ios::left);
  xyzfile << setprecision(15);

  xyzfile << natomsg << endl << energy/627.5 << endl;
  for (int i=0;i<natomsg;i++)
    xyzfile << anames[i] << " " << coords[3*i+0] << " " << coords[3*i+1] << " " << coords[3*i+2] << endl;

  xyzfile.close();

  ofstream gradfile;
  string gradfile_string = filename+".grad";
  gradfile.open(gradfile_string.c_str());
  gradfile.setf(ios::fixed);
  gradfile.setf(ios::left);
  gradfile << setprecision(15);

  gradfile << natomsg << endl << energy/627.5 << endl;
  for (int i=0;i<natomsg;i++)
    gradfile << anames[i] << " " << grad[3*i+0] << " " << grad[3*i+1] << " " << grad[3*i+2] << endl;

  gradfile.close();

  return;
}

void Gradient::init(string infilename, int natoms0, int* anumbers0, string* anames0, int run, int rune, int ncpu, int knnr_level)
{

#if QCHEM && USE_ASE
  printf(" Cannot use both QCHEM and ASE \n");
  exit(-1);
#endif

  V0 = 0.;

  xyz_grad = 0;
  always_do_exact = 0;
  write_on = 1;
  wrote_grad = 0;
  res_t = 0;
  gradcalls = 0;
  nscffail = 0;
  fcounter = 0;

  natoms = natoms0;
  natomsg = natoms;
  anumbers = new int[natoms+1];
  anames = new string[natoms+1];

  for (int i=0;i<natoms0;i++)
    anumbers[i] = anumbers0[i];
  for (int i=0;i<natoms0;i++)
    anames[i] = anames0[i];

  for (int i=0;i<natoms0;i++)
  if (anumbers[i]<1)
    natomsg--;

  runNum = run;
  runend = rune;
  string nstr = StringTools::int2str(run,4,"0");
  runName0 = StringTools::int2str(runNum,4,"0")+"."+StringTools::int2str(runend,4,"0");
  runends = nstr;

#if QCHEM
  qchem1.init(infilename,natomsg,anumbers,anames,run,rune);
  qchem1.ncpu = ncpu;
#endif
#if USE_ASE
  ase1.init(infilename,natomsg,anumbers,anames,run,rune);
  ase1.ncpu = ncpu;
#endif

  knnr_active = 0;
  knn_k = KNN_K;
#if USE_KNNR
  if (knnr_level)
  {
    knnr1.printl = 0;
    int kfound = 0;
#if 1
    knnr_inited = 0;
#else
    printf("  loading KNNR files \n");
    kfound = knnr1.begin(runNum,natomsg);
#endif
    knnr_active = knnr_level;
   // knnr1.test_points();
    if (kfound>0)
      res_t = 1000;
  }
#endif


#if QCHEM
  printf("  grad initiated: Q-Chem mode \n");
#elif USE_ASE
  printf("  grad initiated: ASE mode \n");
#else
  printf("  grad initiated: Mopac mode \n");
#endif

  //printf(" grad init knnr_active: %i \n",knnr_active);

  return;
}



