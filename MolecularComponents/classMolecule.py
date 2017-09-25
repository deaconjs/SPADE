import sys
import os
sys.path.append(os.getcwd())
from MolecularComponents.classAtom import Atom
from MolecularComponents.classPoint import Point
from MolecularComponents.classFutamuraHash import FutamuraHash
import string
import math
import parms
verbose = 0
import copy

class Molecule:
    def __init__(self, PDBlines):
        self.atoms = []
        # create Atom objects from the lines
        # calculate the centroid at the same time
        centroid = [0.0,0.0,0.0]
        self.pdb_lines = PDBlines
        for line in PDBlines:
            new_atom = Atom(self, line)
            new_atom.data['parent_molecule'] = self
            self.add_atom(new_atom)
            centroid[0] = centroid[0] + new_atom.x
            centroid[1] = centroid[1] + new_atom.y
            centroid[2] = centroid[2] + new_atom.z
        self.centroid = Point(centroid[0]/len(self.atoms),centroid[1]/len(self.atoms),centroid[2]/len(self.atoms))

        line = PDBlines[0]
        self.chain_name = string.strip(line[21:22])
        self.res_type   = string.strip(line[17:20])
        self.selected   = 1             # start everything out as selected
        self.visible    = 1             # and visible
        self.atom_points = 'None'   # holds vtkPoint
        # stuff for distances for bonds
        self.C_types  = {'C':1.80,'O':1.65,'N':1.70,'S':2.05,'P':2.10,'H':1.35,'Z':2.5}
        self.O_types  = {'C':1.65,'O':1.45,'N':1.65,'S':1.70,'P':1.85,'H':1.25,'Z':2.5}
        self.N_types  = {'C':1.70,'O':1.65,'N':1.70,'S':2.05,'P':2.10,'H':1.30,'Z':2.5}
        self.S_types  = {'C':2.05,'O':1.70,'N':2.05,'S':1.75,'P':2.10,'H':1.65,'Z':2.5}
        self.P_types  = {'C':2.10,'O':1.85,'N':2.10,'S':2.10,'P':2.30,'H':1.75,'Z':2.5}
        self.H_types  = {'C':1.35,'O':1.25,'N':1.30,'S':1.65,'P':1.75,'H':1.05,'Z':2.5}
        self.Z_types  = {'C':2.5 ,'O':2.50,'N':2.50,'S':2.50,'P':2.50,'H':2.50,'Z':2.5}
        self.common_atoms_list = ['C','O','N','S','P','H','Z']
        self.vtk_arg_list = {}
        volume_args = {'visualize':0,
                       'currently_on':0,
                       'specular':0.3,
                       'specular_power':30,
                       'representation':'surface'}
        atoms_args = {'visualize':0,
                      'currently_on':0,
                      'width':0.1,
                      'sticks_width':0.1,
                      'sticks_sides':9,
                      'representation':'wireframe'}

        self.vtk_arg_list['volume'] = volume_args
        self.vtk_arg_list['atoms'] = atoms_args
        self.x_table = None

        self.features = {}
        self.data     = {}

    def add_atom(self, atom):
        if self.atoms:
            self.atoms.append(atom)
        else:
            self.atoms = []
            self.atoms.append(atom)
            
    def get_bounds(self):
        x = 0.0
        y = 0.0
        z = 0.0
        count = 0
        for atom in self.atoms:
            x = x + atom.x
            y = y + atom.y
            z = z + atom.z
            count = count + 1
                    
        center_x = x / float(count)
        center_y = y / float(count)
        center_z = z / float(count)
        
        min_x = max_x = center_x
        min_y = max_y = center_y
        min_z = max_z = center_z
        
        for atom in self.atoms:
            if atom.x >= max_x:
                max_x = atom.x
            if atom.x < min_x:
                min_x = atom.x
            if atom.y >= max_y:
                max_y = atom.y
            if atom.y < min_y:
                min_y = atom.y
            if atom.z >= max_z:
                max_z = atom.z
            if atom.z < min_z:
                min_z = atom.z
                
        return [min_x,max_x,min_y,max_y,min_z,max_z]

    def len(self):
        if self.atoms:
            return len(self.atoms)
        else:
            return 0

    def pdb_print(self):
        for atom in self.atoms:
            atom.pdb_print()

    def get_pdb_lines(self):
        lines = []
        for atom in self.atoms:
            lines.append(atom.get_pdb_line())
        return lines
    
    def select(self):
        self.selected = 1
        for atom in self.atoms:
            atom.select()

    def deselect(self):
        self.selected = 0
        for atom in self.atoms:
            atom.deselect()

    def get_bonds_list(self, target_atom, hydrogens_on=1):
        ai1 = 0
        # first find this atom, from the stored index
        for index in range(0, len(self.atoms)):
            if self.atoms[index].atom_number == target_atom.atom_number:
                ai1 = index
        return self.get_bonds_list_from_int(ai1, hydrogens_on)

    # ai1 should be the integer index of the atom in its parent monomer or small molecule        
    def get_bonds_list_from_int(self, ai1, hydrogens_on=1):
        x1 = self.atoms[ai1].x
        y1 = self.atoms[ai1].y
        z1 = self.atoms[ai1].z
        bonds_list = []
        # add a bond to the N-terminal amino acid, if present
        if self.__module__ == "MolecularComponents.classAminoAcid":
            if self.atoms[ai1].atom_type == 'N' and self.is_Nterm == 0:
                last_ind = self.parent.res_number_to_residues_index(self.res_number)-1
                for atom in self.parent.residues[last_ind].atoms:
                    if atom.atom_type == 'C':
                        x2 = atom.x
                        y2 = atom.y
                        z2 = atom.z
                        # make a coordinate halfway between them
                        dis = math.sqrt((x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0)
                        if dis < self.min_bonded_distance(atom.atom_type[0], self.atoms[ai1].atom_type[0]):
                            if x1 >= x2:
                                x = x1-(0.5*(x1-x2))
                            else:
                                x = x1+(0.5*(x2-x1))
                            if y1 >= y2:
                                y = y1-(0.5*(y1-y2))
                            else:
                                y = y1+(0.5*(y2-y1))
                            if z1 >= z2:
                                z = z1-(0.5*(z1-z2))
                            else:
                                z = z1+(0.5*(z2-z1))
                            # and store it in bonds_list
                            bonds_list.append([x,y,z])
            if self.atoms[ai1].atom_type == 'C' and self.is_Cterm == 0:
                next_ind = self.parent.res_number_to_residues_index(self.res_number)+1
                for atom in self.parent.residues[next_ind].atoms:
                    if atom.atom_type == 'N':
                        x2 = atom.x
                        y2 = atom.y
                        z2 = atom.z
                        dis = math.sqrt((x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0)
                        # make a coordinate halfway between them
                        if dis < self.min_bonded_distance(atom.atom_type[0], self.atoms[ai1].atom_type[0]):
                            if x1 >= x2:
                                x = x1-(0.5*(x1-x2))
                            else:
                                x = x1+(0.5*(x2-x1))
                            if y1 >= y2:
                                y = y1-(0.5*(y1-y2))
                            else:
                                y = y1+(0.5*(y2-y1))
                            if z1 >= z2:
                                z = z1-(0.5*(z1-z2))
                            else:
                                z = z1+(0.5*(z2-z1))
                            # and store it in bonds_list
                            bonds_list.append([x,y,z])
        # now do the internal bonds
        atom1 = self.atoms[ai1].atom_type[0]
        for ai2 in range(0,len(self.atoms)):
            if hydrogens_on == 0 and (self.atoms[ai2].atom_type[0] == 'H' or atom1[0] == 'H'):
                continue
            if ai1 != ai2:
                x2 = self.atoms[ai2].x
                y2 = self.atoms[ai2].y
                z2 = self.atoms[ai2].z
                dis = math.sqrt((x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0)
                # make a coordinate halfway between them
                atom2 = self.atoms[ai2].atom_type[0]
                if (atom1 not in self.common_atoms_list) or (atom2 not in self.common_atoms_list):
                  disthresh =  1.8
                else:
                    # dictionary style makes the following switch operation more concise
                    thresh_array = {'C':self.C_types[atom2],
                                    'H':self.H_types[atom2],
                                    'O':self.O_types[atom2],
                                    'N':self.N_types[atom2],
                                    'S':self.S_types[atom2],
                                    'P':self.P_types[atom2],
                                    'Z':self.Z_types[atom2]}
                    disthresh = thresh_array[atom1]

                if dis < disthresh:
                    if x1 >= x2:
                        x = x1-(0.5*(x1-x2))
                    else:
                        x = x1+(0.5*(x2-x1))
                    if y1 >= y2:
                        y = y1-(0.5*(y1-y2))
                    else:
                        y = y1+(0.5*(y2-y1))
                    if z1 >= z2:
                        z = z1-(0.5*(z1-z2))
                    else:
                        z = z1+(0.5*(z2-z1))
                    # and store it in bonds_list
                    bonds_list.append([x,y,z])
        return bonds_list                
                
    def min_bonded_distance(self, atom1, atom2):
        common_atoms_list = ['C','O','N','S','P','H','Z']
        if atom1 not in common_atoms_list or atom2 not in common_atoms_list:
            return 1.8
        
        C_types  = {'C':1.80,'O':1.65,'N':1.70,'S':2.05,'P':2.10,'H':1.35,'Z':2.5}
        O_types  = {'C':1.65,'O':1.45,'N':1.65,'S':1.70,'P':1.85,'H':1.25,'Z':2.5}
        N_types  = {'C':1.70,'O':1.65,'N':1.70,'S':2.05,'P':2.10,'H':1.30,'Z':2.5}
        S_types  = {'C':2.05,'O':1.70,'N':2.05,'S':1.75,'P':2.10,'H':1.65,'Z':2.5}
        P_types  = {'C':2.10,'O':1.85,'N':2.10,'S':2.10,'P':2.30,'H':1.75,'Z':2.5}
        H_types  = {'C':1.35,'O':1.25,'N':1.30,'S':1.65,'P':1.75,'H':1.05,'Z':2.5}
        Z_types  = {'C':2.5 ,'O':2.50,'N':2.50,'S':2.50,'P':2.50,'H':2.50,'Z':2.5}

        # concise switch statement style
        return_array = {'C':C_types[atom2],
                        'O':O_types[atom2],
                        'N':N_types[atom2],
                        'S':S_types[atom2],
                        'P':P_types[atom2],
                        'H':H_types[atom2],
                        'Z':Z_types[atom2]}
        
        return return_array[atom1]

    def get_central_pt(self):
        # for use with classes amino acid and nucleotide only
        new_pt = Point(self.central_pt.x, self.central_pt.y, self.central_pt.z)
        return new_pt
    
    #normally applied to monomers, ligands, and waters - coor_string is a list of triples in the
    # same order as the atoms are in the pdb.
    # format -
    # res_type res_number chain_name float float float ...
    def update_cfs(self, coor_string, distances_moved):
        # check to see that the type and number are the same
        coor_list = string.split(coor_string.strip())
        if coor_list[0] != self.res_type or int(coor_list[1]) != self.res_number:
            if len(string.strip(self.chain_name)) > 0:
                if coor_list[2] != self.chain_name:
                    print "problem loading: %s%d and %s%d different"%(coor_list[0], int(coor_list[1]), self.res_type, self.res_number)
                    return -1
            else:
                print "problem loading: %s%d and %s%d different"%(coor_list[0], int(coor_list[1]), self.res_type, self.res_number)
                return -1
        if len(string.strip(self.chain_name)) == 0:
            offset = 2
        else:
            offset = 3
        if len(coor_list) != (len(self.atoms)*3)+offset:
            print "problem loading: %d and %d different numbers of atoms?"%(len(coor_list), (len(self.atoms)*3)+offset)
            return -1
        distances = []
        for i in range(0,(len(self.atoms))):
            distances_moved.append(math.sqrt(((self.atoms[i].x-float(coor_list[i*3+offset]))**2)+
                                             ((self.atoms[i].y-float(coor_list[i*3+offset+1]))**2)+
                                             ((self.atoms[i].z-float(coor_list[i*3+offset+2]))**2)))
            self.atoms[i].x = float(coor_list[i*3+offset])
            self.atoms[i].y = float(coor_list[i*3+offset+1])
            self.atoms[i].z = float(coor_list[i*3+offset+2])

    def build_futamura_intersection_table(self, solvent_radius):
        outside_barrier = 4.0
        grid_spacing = 4*solvent_radius
        grid = FutamuraHash(self, outside_barrier, grid_spacing)
        T = grid.T
        block_assignments = grid.atom_block_assignments
        print 'locating intersections'
        # now locate the intersections
        x_table = {}
        radii  = {'C':solvent_radius + 1.75,
                  'O':solvent_radius + 1.4,
                  'N':solvent_radius + 1.55,
                  'S':solvent_radius + 1.8,
                  'P':solvent_radius + 2.0,
                  'H':solvent_radius + 1.17,
                  'Z':solvent_radius + 3.0}
        default_distance = solvent_radius + 1.8
            
        for atom in self.atoms:
            x_table['%s'%(atom.atom_number)] = []
            r1 = radii.get(atom.atom_type[0], default_distance)
            block = block_assignments['%s'%(atom.atom_number)]
            key_tokens = string.split(block)
            keys = [string.atoi(key_tokens[0]), string.atoi(key_tokens[1]), string.atoi(key_tokens[2])]
            # put 'this' block first, so that intersection table accesses search here first
            for second_atom in T[block]:
                if atom != second_atom:
                    r2 = radii.get(second_atom.atom_type[0], default_distance)
                    if atom.dist(second_atom) <= r1+r2:
                        x_table['%s'%(atom.atom_number)].append([second_atom.x,second_atom.y,second_atom.z, r2, second_atom.res_number, second_atom.chain_name, second_atom.atom_type])
            
            start_array = [0,0,0]
            end_array   = [0,0,0]
            # figure out starts and ends
            counts = [grid.volume_count_x, grid.volume_count_y, grid.volume_count_z]
            for ind in [0,1,2]:
                if keys[ind] == 0:
                    start_array[ind] = 0
                    end_array[ind]   = 2
                elif keys[ind] == counts[ind] - 1:
                    start_array[ind] = keys[ind]-1
                    end_array[ind]   = keys[ind]+1
                else:
                    start_array[ind] = keys[ind]-1
                    end_array[ind]   = keys[ind]+2

            for i in range(start_array[0], end_array[0]):
                for j in range(start_array[1], end_array[1]):
                    for k in range(start_array[2], end_array[2]):
                        key2 = '%s %s %s'%(i,j,k)
                        if key2 == block:
                            continue            # did this one earlier
                        if key2 in T.keys():
                            for second_atom in T[key2]:
                                if atom != second_atom:
                                    r2 = radii.get(second_atom.atom_type[0], default_distance)
                                    if atom.dist(second_atom) <= r1+r2:
                                        x_table['%s'%(atom.atom_number)].append([second_atom.x,second_atom.y,second_atom.z, r2, second_atom.res_number, second_atom.chain_name, second_atom.atom_type])
        self.x_table = x_table

    def fill_densities(self, distance_cutoff=9.0, force_rewrite=0):
        filename = self.parent.get_filename_by_extension('.ads', self.chain_name)
        create_new = 0
        if force_rewrite:
            create_new = 1
        else:
            try:
                densities_file = open(filename)
            except IOError:
                create_new = 1
        if create_new:
            if self.x_table == None:
                self.build_futamura_intersection_table(distance_cutoff)
            print 'filling densities'
            atom_count = 0.0
            densities = []
            max_count = 0
            for atom in self.atoms:
                data_list = self.x_table['%s'%(atom.atom_number)]
                atom_count = 0
                for data_ind in range(len(data_list)):
                    atom_data = self.x_table['%s'%(atom.atom_number)]
                    p = Point(data_list[data_ind][0], data_list[data_ind][1], data_list[data_ind][2])
                    d = atom.dist(p)
                    if d < distance_cutoff:
                        atom_count = atom_count + 1
                densities.append(atom_count)
                if atom_count > max_count:
                    max_count = atom_count
            for i in range(len(self.atoms)):
                self.atoms[i].features['atomic_density'] = (densities[i]+0.0)/max_count
                
            densities_file = open(filename, 'w')
            for atm in self.atoms:
                densities_file.write('%s\n'%(atm.features['atomic_density']))
            densities_file.close()
        else:           # else read the contacts_file to fill the contact_list
            for atom in self.atoms:
                buffer = densities_file.readline()
                if len(buffer) == 0:
                    break
                atom.features['atomic_density'] = string.atof(buffer)
            densities_file.close()


