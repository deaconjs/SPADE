import sys
import os
sys.path.append(os.getcwd())
from MolecularComponents.classMolecule import Molecule
import string
import parms

class Water(Molecule):
    def __init__(self, PDBlines):
        Molecule.__init__(self, PDBlines)
        self.res_number = string.atoi(PDBlines[0][23:26])
        self.key = '%s%s'%(self.res_number, self.chain_name)
