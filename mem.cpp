// Please see license.txt for licensing and copyright information //
// Author: Paul Zimmerman, University of Michigan //
#include "icoord.h"

void ICoord::freemem(){

    delete [] grad;

    for (int i=0;i<max_bonds;i++)
        delete [] bonds[i];
    // delete [] bonds;

    for (int i=0;i<max_angles;i++)
        delete [] angles[i];
    // delete [] angles; 
    // delete [] anglev;

    for (int i=0;i<max_nonbond;i++)
        delete [] nonbond[i];

    // delete [] coordn;


    delete [] anumbers;
    delete [] amasses;
    delete [] anames;
    delete [] coords0;
    delete [] coordsts;
    delete [] coords;

    delete [] xyzic;

    return;

}


void ICoord::alloc_mem(){

    farBond = 1.0;
    use_xyz = 0;

    frozen = NULL;
    xyzic = new int[natoms];

    // pgrad = new double[3*natoms];
    grad = new double[3*natoms];
    // for (int i=0;i<3*natoms;i++) pgrad[i]=0;
    for (int i=0;i<3*natoms;i++) grad[i]=0;

    nbonds = 0;
    max_bonds=natoms*natoms+8;
    bonds = new int*[max_bonds];
    for (int i=0;i<max_bonds;i++)
        bonds[i]=new int[2];
    bondd = new double[max_bonds];
    for (int i=0; i<max_bonds; i++)
    {
        bondd[i] = 0.0;
        for (int j=0; j<2; j++)
        {
            bonds[i][j] = 0;
        }
    }

    nangles = 0;
    max_angles=natoms*24;
    angles = new int*[max_angles];
    for (int i=0;i<max_angles;i++)
        angles[i]=new int[3];
    for (int i=0; i<max_angles; i++)
    {
        for (int j=0; j<3; j++)
        {
            angles[i][j] = 0;
        }
    }

    anglev = new double[max_angles];
    for (int i=0; i<max_angles; i++)
    {
        anglev[i] = 0.0;
    }

    ntor = 0;
    max_torsions=natoms*80;
    torsions = new int*[max_torsions];
    for (int i=0;i<max_torsions;i++)
        torsions[i] = new int[4];
    torv = new double[max_torsions];

    nimptor = 0;
    max_imptor=natoms*10;
    imptor = new int*[max_imptor];
    for (int i=0;i<max_imptor;i++)
        imptor[i] = new int[4];
    imptorv = new double[max_imptor];


    n_nonbond = 0;
    max_nonbond = natoms*natoms;
    nonbond = new int*[max_nonbond];
    for (int i=0;i<max_nonbond;i++)
        nonbond[i]=new int[2];
    nonbondd = new double[max_nonbond];

    coordn = new int[natoms];
    for (int i=0;i<natoms;i++)
        coordn[i]=0;

    ffR = new double[natoms];
    ffeps = new double[natoms];

    return;
}

