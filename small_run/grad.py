#!/usr/bin/env python
import sys 
#sys.path.append('/export/zimmerman/paulzim/ase')
# this is for parallel gsm, so VASP does not read numOfThreads from submit script
import os
os.environ['OMP_NUM_THREADS'] = '1'

from ase.lattice.surface import surface
from ase import Atoms,Atom
from ase.visualize import view
from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import write
from ase.io import read
from ase.lattice.surface import fcc111,add_adsorbate,fcc100,fcc110,hcp0001,bcc111
import os
import sys
import shutil

nargs = len(sys.argv)
argv1 = sys.argv[1]
argv2 = sys.argv[2]

os.system("cp scratch/structure"+argv1+" scratch/structure"+argv1+".xyz")
fname = 'scratch/structure'+argv1+".xyz"
tmpDir = "."#os.getenv('PBSTMPDIR')
folder = 'scratch/vasp'+argv1

#lattice and initial atom setup for unit cell
#slab = read("initial.xyz")
#slab.set_cell([12.720823288,12.720823288,21.0])
#slab2 = hcp0001('Ru', size=(3,3,2), vacuum=10)
#define adsorbate
#molecule = Atoms('CO', [(0., 0., 0.), (0., 0., 1.16)])
#add adsorbate to the surface
#add_adsorbate(slab2, molecule, 1.7, 'hcp')

slab = fcc111('Cu', size=(2,2,2), vacuum=5.0)
molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., 1.1)])
add_adsorbate(slab, molecule, 1.85, 'ontop')

#print("slab1", slab1)
#print("slab2", slab2)
#sys.exit()

if os.path.isfile(fname):
  #print("\nfname:", fname)
  #current position read in
  slabatoms = read(fname)
  #print("after read in")
  slab.set_positions(slabatoms.get_positions())
  #for atom in slab:
  #  print("atom:", atom)

#mask = [atom.tag > 2 for atom in slab] # MB: this needs to be adapted!!!
#slab.set_constraint(FixAtoms(mask=mask))
print "Now submit VASP",argv1,argv2 
calc = Vasp(xc='PBE',lreal='Auto',kpts=[1,1,1],ispin=2,ismear=0,sigma=0.01,algo='fast',istart=0,ediff=0.00001,lwave = False,lcharg = False) #npar=8
slab.set_calculator(calc)

cwd = os.getcwd()

os.chdir(tmpDir)
if not os.path.exists(folder):
    os.makedirs(folder)
os.chdir(folder)
energy = - slab.get_potential_energy()
grads = - slab.get_forces()

f = open('GRAD'+argv1, 'w')
f.write(str(energy))
f.write('\n')
f.write(str(grads))
f.write('\n')
f.close()

shutil.copy2('GRAD'+argv1, cwd+'/scratch')
os.chdir(cwd)
