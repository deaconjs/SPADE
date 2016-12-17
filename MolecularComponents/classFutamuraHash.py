import math
import string

class FutamuraHash:
    def __init__(self, mol, outside_barrier=4.0, grid_spacing=4*1.4, alpha_carbon_code=0, core_cutoff=5, neighbor_thresh=3):
        """ alpha_carbon_code is 0 if all atoms are added to the hash table,
            is 1 if all alpha carbons are added
            and is 2 if only the core alpha carbons are added
        """
        atom_list = []
        for atom in mol.atoms:
            atom_list.append(atom)
            
        # first figure out what the bounding box is for the molecule
        self.molecule = mol
        min_x = max_x = mol.centroid.x
        min_y = max_y = mol.centroid.y
        min_z = max_z = mol.centroid.z
        for atom in mol.atoms:
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
        # create an empty barrier around the bounding box, so that contours are meaningful
        # at the boundaries
        min_x = min_x - outside_barrier
        max_x = max_x + outside_barrier
        min_y = min_y - outside_barrier
        max_y = max_y + outside_barrier
        min_z = min_z - outside_barrier
        max_z = max_z + outside_barrier
        
        # create a grid, starting at the minimum corner, with spacing of distance grid_spacing
        # figure out how many will be required (round up)
        self.volume_count_x = int(math.ceil((max_x - min_x) / grid_spacing))
        self.volume_count_y = int(math.ceil((max_y - min_y) / grid_spacing))
        self.volume_count_z = int(math.ceil((max_z - min_z) / grid_spacing))
        
        # update the bounding box corners
        max_x = min_x + (self.volume_count_x * grid_spacing)
        max_y = min_y + (self.volume_count_y * grid_spacing)
        max_z = min_z + (self.volume_count_z * grid_spacing)
        # these indices hold the coordinates being used 
        self.volume_indices_x = []
        self.volume_indices_y = []
        self.volume_indices_z = []
        for i in range(self.volume_count_x):
            self.volume_indices_x.append(min_x + (i*grid_spacing))
        for i in range(self.volume_count_y):
            self.volume_indices_y.append(min_y + (i*grid_spacing))
        for i in range(self.volume_count_z):
            self.volume_indices_z.append(min_z + (i*grid_spacing))

        # T is a hash table that will store all of the atoms by their grid blocks
        # block_assignments stores keys to T by atom number (which should be unique)
        self.T = {}
        self.atom_block_assignments = {}
        if alpha_carbon_code == 0:
            atom_list = mol.atoms
        elif alpha_carbon_code == 1:
            atom_list = mol.get_central_atom_list()
        elif alpha_carbon_code == 2:
            atom_list = mol.get_core_alpha_carbons(core_cutoff, neighbor_thresh)
        for atom in atom_list:
            x_val = y_val = z_val = 0
            # figure out which bin it goes in
            for x_ind in range(0, self.volume_count_x):
                if atom.x < self.volume_indices_x[x_ind]:
                    break
                else:
                    x_val = x_ind
            for y_ind in range(self.volume_count_y):
                if atom.y < self.volume_indices_y[y_ind]:
                    break
                else:
                    y_val = y_ind
            for z_ind in range(self.volume_count_z):
                if atom.z < self.volume_indices_z[z_ind]:
                    break
                else:
                    z_val = z_ind
            
            # if the bin is empty, add a list and append this atom
            key = '%s %s %s'%(x_val, y_val, z_val)
            if key not in (self.T.keys()):
                self.T[key] = []
                self.T[key].append(atom)
            else:
                self.T[key].append(atom)
            self.atom_block_assignments['%s'%(atom.atom_number)] = key

    def get_grid_blocks_to_search(self, atom_number):
        block = self.atom_block_assignments['%s'%(atom_number)]
        key_tokens = string.split(block)
        keys = [string.atoi(key_tokens[0]), string.atoi(key_tokens[1]), string.atoi(key_tokens[2])]
        start_array = [0,0,0]
        end_array   = [0,0,0]
        # figure out starts and ends
        counts = [self.volume_count_x, self.volume_count_y, self.volume_count_z]
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
        search_keys = []
        for i in range(start_array[0], end_array[0]):
            for j in range(start_array[1], end_array[1]):
                for k in range(start_array[2], end_array[2]):
                    key2 = '%s %s %s'%(i,j,k)
                    search_keys.append(key2)
        return search_keys

    def get_nearby_atoms(self, atom_number, core_only=0):
        active_blocks = self.get_grid_blocks_to_search(atom_number)
        nearby_atoms = []
        for block_key in active_blocks:
            if block_key in self.T.keys():
                for second_atom in self.T[block_key]:
                    if second_atom.atom_number != atom_number:
                        nearby_atoms.append(second_atom.atom_number)
        return nearby_atoms

    def get_atoms_between(self, atom_number, min_range, max_range, core_only=0):
        """ returns a list of dictionaries with atom_number and distance keys
        """
        nearby_atoms = self.get_nearby_atoms(atom_number)
        return_atoms = []
        atom = self.molecule.atom_dict[atom_number]
        for other_atom_number in nearby_atoms:
            other_point = self.molecule.atom_dict[other_atom_number]
            dst = math.sqrt(((atom.x-other_point.x)**2) + ((atom.y-other_point.y)**2) + ((atom.z-other_point.z)**2))
            if dst >= min_range and dst <= max_range:
                the_atom = self.molecule.atom_dict[other_atom_number]
                return_atoms.append(the_atom)
        return return_atoms






