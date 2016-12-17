import sys
import os
sys.path.append(os.getcwd())
from MolecularComponents.classMolecule import Molecule
from MolecularComponents.classFutamuraHash import FutamuraHash
import string
import parms
    
class Polymer(Molecule):
    def __init__(self, PDBlines, parent):
        Molecule.__init__(self, PDBlines)
        self.parent = parent
        self.residues = []
        self.residue_dict = {}
        self.selected = 1
        self.key = self.chain_name
        self.default_chain_color = self.get_chain_color_word()
        trace_args   = {'visualize':1,
                       'currently_on':0,
                       'representation':'line',
                       'width':0.4,
                       'splines':6,
                       'sides':5,
                       'specular':0.1,
                       'specular_power':10}
        self.vtk_arg_list['trace']  = trace_args

    def add_residue(self, residue):
        self.residues.append(residue)
        self.residue_dict[residue.res_number] = residue
        
    def len(self):
        if self.residues:
            return len(self.residues)
        else:
            return 0

    def pdb_print(self):
        for rez in self.residues:
            rez.pdb_print()

    def select(self):
        Molecule.select(self)
        self.selected = 1
        for rez in self.residues:
            rez.select()
            
    def deselect(self):
        Molecule.deselect(self)
        self.selected = 0
        for rez in self.residues:
            rez.deselect()

    def get_chain_color(self):                          # returns three doubles
        if parms.get('colorblind'):
            colormap = parms.get('colorblind_colormap')
            colormap_count = len(colormap)
        else:
            colormap = parms.get('colormap')
            colormap_count = len(colormap)
        index = 0
        for pchain in self.parent.ProteinList:
            if self.chain_name == pchain.chain_name:
                color_name = 'chain%s'%(index%12+1)
                return colormap[color_name]
            index = index + 1
        for nchain in self.parent.NucleotideChainList:
            if self.chain_name == nchain.chain_name:
                color_name = 'chain%s'%(index%12+1)
                return colormap[color_name]
            index = index + 1
            
    def get_chain_color_word(self):                          # returns the word color
        if parms.get('colorblind'):
            colormap = parms.get('colorblind_colormap')
        else:
            colormap = parms.get('colormap')
        index = 0
        for pchain in self.parent.ProteinList:
            if self.chain_name == pchain.chain_name:
                return 'chain%s'%(index%12+1)
            index = index + 1
        for nchain in self.parent.NucleotideChainList:
            if self.chain_name == nchain.chain_name:
                return 'chain%s'%(index%12+1)
            index = index + 1

    def locate_termini(self):
        self.residues[0].is_Nterm = 1
        self.residues[len(self.residues)-1].is_Cterm = 1
        if self.__module__ == 'MolecularComponents.classProtein':
            max_adjacent_distance = 4.0
        elif self.__module__ == 'MolecularComponents.classNucleotideChain':
            max_adjacent_distance = 8.0
        for rez_index in range(1, len(self.residues)-1):
            # if the number of the last residue is not one less than that of the current
            if self.residues[rez_index].res_number != self.residues[rez_index-1].res_number+1:
                # then its an N-terminus
                self.residues[rez_index].is_Nterm = 1
            # if the next residue is more than some # of angstroms away its probably a break in the chain
            # this distance has not been thoroughly examined
            # first make sure both residues have central atoms
            if self.residues[rez_index].has_central_pt and self.residues[rez_index-1].has_central_pt:
                if self.residues[rez_index].central_pt.get_distance_from(self.residues[rez_index-1]) > max_adjacent_distance:
                    self.residues[rez_index].is_Nterm = 1
            elif self.residues[rez_index-1].has_central_pt == 0:
                # if, however, the one N-terminal to this one doesn't even have an alpha carbon,
                # then this one is an N-terminus
                self.residues[rez_index].is_Nterm = 1
            if self.residues[rez_index].res_number != self.residues[rez_index+1].res_number-1:
                # then its an C-terminus
                self.residues[rez_index].is_Cterm = 1
            # now test the distance to the last residue
            if self.residues[rez_index].has_central_pt and self.residues[rez_index+1].has_central_pt:
                if self.residues[rez_index].central_pt.get_distance_from(self.residues[rez_index+1]) > max_adjacent_distance:
                    self.residues[rez_index].is_Cterm = 1
            elif self.residues[rez_index+1].has_central_pt == 0:
                # if, however, the one N-terminal to this one doesn't even have an alpha carbon,
                # then this one is an N-terminus
                self.residues[rez_index].is_Cterm = 1

    def res_number_to_residues_index(self, res_number):
        """ This function is made obsolete by the polymer's residue dictionary
        """
        index = 0
        for rez in self.residues:
            if rez.res_number == res_number:
                return index
            index = index + 1
        else:
            return -1

    def assign_atom_distances(self, mindist=0.0, maxdist=12.0):
        f_hash = FutamuraHash(self, 0.0, maxdist, 1)
        atom_list = self.get_central_atom_list()
        for cp in atom_list:
            cp.data['nearby_cps'] = f_hash.get_atoms_between(cp.atom_number, mindist, maxdist)

    def get_central_point_list(self):                              # returns a list of point objects
        central_point_list = []
        for rez in self.residues:
            if rez.has_central_pt==1:
                central_point_list.append(rez.central_pt)
            else:
                print 'missing central atom for rez %s'%(rez.res_number)
        return central_point_list

    def get_central_atom_list(self):
        central_atom_list = []
        for rez in self.residues:
            if rez.has_central_pt==1:
                central_atom_list.append(rez.central_atom)
        return central_atom_list
    




