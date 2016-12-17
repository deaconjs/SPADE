import sys
sys.path.append('../..')
from MolecularSystem import System
from BlissSuperimpose import superimpose
x = System(None)
y = System(None)
x.load_pdb('atoms1.pdb')
y.load_pdb('atoms2.pdb')
superimpose([x,y])
x.save_pdb('atoms3.pdb')
y.save_pdb('atoms4.pdb')




