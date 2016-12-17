import sys
import os
sys.path.append(os.getcwd())
from MolecularComponents.classMolecule import Molecule
from MolecularComponents.classPoint import Point
import string

class Nucleotide(Molecule):
    def __init__(self, PDBlines, parent, atoms):
        Molecule.__init__(self, PDBlines)
        self.parent = parent
        self.atoms = atoms
        for atom in atoms:
            atom.parent = self
        line = PDBlines[0]
        self.res_type1  = string.strip(line[17:20])
        self.res_number = string.atoi(line[23:26])

        self.is_Nterm = 0
        self.is_Cterm = 0
        self.has_central_pt = 0    
        self.atoms_dict = {}
        for atom in self.atoms:
            self.atoms_dict[atom.atom_type] = atom
            if atom.atom_type[:2] == 'C3':
                self.central_pt = Point(atom.x, atom.y, atom.z)
                self.central_atom = atom
                self.x=atom.x
                self.y=atom.y
                self.z=atom.z
                self.has_central_pt = 1
        self.is_Nterm = 0
        self.is_Cterm = 0
        self.vtk_arg_list = {}
        

        trace_args = {'color':[0.1,0.1,1.0],
                      'opacity':1.0}
        self.vtk_arg_list['trace'] = trace_args
        self.visible = 0
        self.transparent = 0

