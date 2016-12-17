import sys
import os
import re
sys.path.append(os.getcwd())

from MolecularComponents.classPoint import Point
import string

class Atom(Point):
    def __init__(self, parent, PDBline):
        self.pdb_line = PDBline
        self.parent      = parent                       # parent molecule
        self.build()
        
    def build(self, initial_build=1):
        # this could probably be rewritten for robustness to learn the
        # useage of column patterns in any particular pdb file
        PDBline = self.pdb_line
        record           = string.strip(PDBline[0:6])   # record name
        self.atom_number = int(PDBline[6:11])   # atom number
        self.atom_type   = string.strip(PDBline[13:16]) # atom type
        # if extra_type is present, test to see if its a number or a letter
        self.extra_type1 = string.strip(PDBline[11:12])
        if self.extra_type1 != "":
            if re.compile('[0-9]').search(self.extra_type1):
                # if its a number, append it to atom_type -- present in Amber protonated files
                self.atom_type = self.atom_type + self.extra_type1
            elif re.compile('[A-Z]').search(self.extra_type1):
                # if a capital letter, its probably the first letter of the type. Just add it there.
                self.atom_type = self.extra_type1 + self.atom_type
        self.extra_type2 = string.strip(PDBline[12:13]) # sometimes a number appears here
        if self.extra_type2 != "":
            if re.compile('[0-9]').search(self.extra_type2):
                self.atom_type = self.atom_type + self.extra_type2
            elif re.compile('[A-Z]').search(self.extra_type2):
                self.atom_type = self.extra_type2 + self.atom_type
        #if PDBline[15:16] == " ":
        #    self.branch_des = 0
        #else:
        #    self.branch_des  = string.atoi(PDBline[15:16])  # branch designator
        #self.alt_ind     = string.strip(PDBline[16:17]) # alternate location indicator
        self.res_type    = string.strip(PDBline[17:20]) # residue type
        self.chain_name  = string.strip(PDBline[21:22]) # chain name
        self.res_number  = int(PDBline[22:26])  # residue number
        self.x           = float(PDBline[30:38])  # x-coor
        y = string.split(PDBline[38:46])
        if len(y) == 2:
            self.y           = float(y[1])  # y-coor
        else:
            self.y           = float(y[0])
        z = string.split(PDBline[46:54])
        if len(z) == 2:
            self.z           = float(z[1])
        else:
            self.z           = float(z[0])
        # set occupancy value, if available
        self.occupancy = None
        if len(PDBline) >= 60:
            token = string.strip(PDBline[54:60])
            if len(token) > 0:
                self.occupancy   = float(token)  # occupancy
        # set the b_factor, if available
        self.b_factor = None
        if len(PDBline) >= 66:
            token = string.strip(PDBline[60:66])
            if len(token) > 0:
                self.b_factor = float(token)
        # set the segment term, if available
        self.segment = None
        if len(PDBline) >= 76:
            token = string.strip(PDBline[73:76])
            if len(token) > 0:
                self.segment     = string.strip(token) # segment identifier
        # set the charge term, if available
        #self.pdb_charge = None
        #if len(PDBline) >= 80:
        #    token = string.strip(PDBline[79:80])
        #    if len(token) > 0:
        #        self.pdb_charge  = string.atof(token)  # charge from the pdb file

        if initial_build:
            self.selected = 1
            self.visible = 0
            self.bonds_list = []        # the list of bonded atoms gets filled in class molecule, in define_connectivity
            self.vtk_arg_list = {}
            radii  = {'C':1.75,'O':1.4,'N':1.55,'S':1.8,'P':2.0,'H':1.17,'Z':3.0}
            self.radius = radii.get(self.atom_type[0], 1.7)
            atoms_args = {'color':[0.8,0.8,0.8],
                          'opacity':1.0}
            volume_args = {'color':atoms_args['color'],
                           'opacity':1.0}
            self.vtk_arg_list['atoms'] = atoms_args
            self.vtk_arg_list['volume'] = volume_args
        
            self.features = {}        
            self.data     = {}
        
    def reduce_for_storage(self):
        line = self.pdb_line
        keep_list = ['parent', 'selected', 'visible', 'vtk_arg_list', 'features', 'data']
        features_del_list = ['b_factor']
        #for key in self.__dict__.keys():
        #    if key not in keep_list:
        #        del self.__dict__[key]
        #del self.features['b_factor']
        self.pdb_line = line
        
    def rebuild_from_storage(self):
        self.build(0)
        
    def get_selected(self):
        return self.selected

    def pdb_print(self):
        if self.res_type in ['WAT', 'HOH', 'H2O', 'IP3']:
            print "HETATM%5d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n"%(self.atom_number,self.atom_type,self.res_type,self.chain_name,self.res_number,self.x,self.y,self.z, " ")
        else:
            print "ATOM %6d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n"%(self.atom_number,self.atom_type,self.res_type,self.chain_name,self.res_number,self.x,self.y,self.z, " ")

    def get_pdb_line(self):
        if self.res_type in ['WAT', 'HOH', 'H2O', 'IP3']:
            return "HETATM%5d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n"%(self.atom_number,self.atom_type,self.res_type,self.chain_name,self.res_number,self.x,self.y,self.z, " ")
        else:
            return "ATOM %6d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n"%(self.atom_number,self.atom_type,self.res_type,self.chain_name,self.res_number,self.x,self.y,self.z, " ")

    def select(self):
        self.selected = 1

    def deselect(self):
        self.selected = 0

    def point_line_distance(self, a1, a2):
        """
        calculate p4, the atom on the line between atoms p1 and p2 to which atom p3 is perpendicular
        then return the distance between p3 and p4.
           1_____3_____2         is relevant, but, 
        3  1___________2         and
           1___________2   3     are not.
        
        returns -1 if p3 is not perpendicular to p1-p2 between them. 
        
        however, if 3 is the alpha carbon of 1 or 2, shielding is surely not as bad when
        the alpha carbon is just to the inside of the line. so, a simple weighting function
        attempts to fix the problem by doubling the distance. Therefore lower shielding is
        reported when considering the alpha carbons of the amino acids being investigated.
        """
        score = Point.point_line_distance(self, p1.central_point, p2.central_point)
        if self.atom_type == 'CA':
            if self.res_number == a1.res_number and self.chain_name == a1.chain_name:
                return 2*score
            elif self.res_number == a2.res_number and self.chain_name == a2.chain_name:
                return 2*score





