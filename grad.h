// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //
#ifndef GRAD_H
#define GRAD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <vector>
#include <cstring>
#include <math.h>

#include "stringtools.h"
#include "qchem.h"
#include "ase.h"
#include "knnr.h"

class Gradient 
{

    private:

        int runNum;
        int runend;
        string runends;
        string runName;
        string runName0;

        int natoms;
        int natomsg;
        string* anames;
        int* anumbers;

        int fcounter;
        QChem qchem1;
        ASE ase1;

        int knn_k;
        KNNR knnr1;
        int knnr_inited;


    public:

        double grads(double* coords, double* grad, double* Ut, int type);
        void init(string infilename, int natoms, int* anumbers, string* anames,
                int run, int rune, int ncpu, int use_knnr);
        void update_knnr();
        void freemem();
        void write_xyz_grad(double* coords, double* grad, string filename);
        int external_grad(double* coords, double* grads);

        int knnr_active;
        int always_do_exact;
        int write_on;
        int wrote_grad;
        int xyz_grad;
        int gradcalls;
        int nscffail;
        double V0;

        double energy0;
        double energy;
        int res_t; //restart found files

};

#endif
