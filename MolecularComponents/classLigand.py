import sys
import os
sys.path.append(os.getcwd())
from MolecularComponents.classMolecule import Molecule
import parms

class Ligand(Molecule):
    def __init__(self, PDBlines, parent):
        Molecule.__init__(self, PDBlines)
        self.parent = parent
        for atom in self.atoms:
            atom.parent = self
        self.res_number = int(PDBlines[0][23:26])
        self.key = '%s%s'%(self.res_number, self.chain_name)

        
