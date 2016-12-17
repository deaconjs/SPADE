import sys
import os
sys.path.append(os.getcwd())

from MolecularComponents.classMolecule import Molecule
from MolecularComponents.classAtom import Atom
from MolecularComponents.classPoint import Point
import string
import copy

class AminoAcid(Molecule):
    def __init__(self, parent, PDBlines, atoms):
        Molecule.__init__(self, PDBlines)
        self.atoms = atoms
        for atom in self.atoms:
            atom.parent = self
        self.parent = parent
        # load amino acid properties
        types = {'ALA':'A',
                 'CYS':'C','CYD':'C','CYX':'C','CYZ':'C',
                 'ASP':'D','GLU':'E','PHE':'F','GLY':'G',
                 'HIS':'H','HID':'H','HIE':'H',
                 'ILE':'I','LYS':'K','LEU':'L',
                 'MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R',
                 'SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}
        self.res_type1  = types[self.res_type]
        line = PDBlines[0]
        self.res_number = string.atoi(line[23:26])
        # store alpha carbon coordinates
        self.has_central_pt = 0
        self.atoms_dict = {}
        for atom in self.atoms:
            self.atoms_dict[atom.atom_type] = atom
            if atom.atom_type == 'CA':
                self.central_pt = Point(atom.x, atom.y, atom.z)
                self.central_atom = atom
                self.x=atom.x
                self.y=atom.y
                self.z=atom.z
                self.has_central_pt = 1
        try:
            self.atoms_dict['CA']
        except KeyError:
            print "Warning: %s%d has no apparent alpha carbon"%(self.res_type, self.res_number)

        types = {'ALA':71.09,'CYS':103.15,'CYD':103.15,'CYX':103.15,'CYZ':103.15,
                 'ASP':115.09,'GLU':129.12,'PHE':147.18,
                 'GLY':57.05,'HIS':137.14,'HID':137.14,'HIE':137.14,'ILE':113.16,'LYS':128.17,'LEU':113.16,
                 'MET':131.19,'ASN':114.11,'PRO':97.12,'GLN':128.14,'ARG':156.19,
                 'SER':87.08,'THR':101.11,'VAL':99.14,'TRP':186.21,'TYR':163.18}
        self.mw = types[self.res_type]
        # these get fixed later on in initialization by class Protein
        self.is_Nterm = 0
        self.is_Cterm = 0
        # create the vtk_args_list
        # [onoff, width, color, splines, sides, specular, specularpower]
        self.vtk_arg_list = {}                  # this gets overridden because residue doesnt need everything in Molecule
        trace_args = {'color':[0.1,0.1,1.0],
                     'opacity':1.0}
        self.vtk_arg_list['trace'] = trace_args
        self.visible = 0
        self.transparent = 0

        
    def get_res_number(self):
        return self.res_number
        
    # beta carbon calculation/accessors
    # header data storage (if sequence)
    def get_average_sidechain_position(self):
        atom_count = 0.0
        avg_x = avg_y = avg_z = 0.0
        for atom in self.atoms:
            if atom.atom_type in ['N','O','CA','H']:
                pass
            else:
                avg_x = avg_x + atom.x
                avg_y = avg_y + atom.y
                avg_z = avg_z + atom.z
                atom_count = atom_count + 1
        if atom_count == 0:
            return [self.central_pt.x,self.central_pt.y, self.central_pt.z]
        else:
            return [avg_x/atom_count,avg_y/atom_count, avg_z/atom_count]
    
    def get_average_atom_feature(self, feature):
        """ returns the average value of the given feature (from the 'features' dict)
            first checks to see that the feature exists, and that it is a number for each atom
            if either check failes, returns None
        """
        average_val = 0.0
        for atom in self.atoms:
            try:
                average_val += string.atof(atom.features[feature])
            except ValueError:
                print 'Feature %s is not a number for atom %s'%(feature, atom.atom_num)
                return None
            except KeyError:
                print 'Feature %s does not exist for atom %s'%(feature, atom.atom_num)
                return None
        return average_val / (len(self.atoms)+0.0)



