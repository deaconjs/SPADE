import sys
import os
import string
import parms
import copy
import numpy
import math
import random
import re
sys.path.append(os.getcwd())
from MolecularComponents.classAminoAcid import AminoAcid
from MolecularComponents.classPolymer import Polymer
from MolecularComponents.classPoint import Point
from MolecularComponents.classFutamuraHash import FutamuraHash
from MolecularComponents.classMolecule import Molecule

from MathFunc import *
from HBondFunc import *

verbose = 0

class Protein(Polymer):
    def __init__(self, PDBlines, parent):
        self.atomsHydList=[]
        Polymer.__init__(self,PDBlines, parent)
        # store a list element for each amino acid in PDBlines
        last_line = None
        current_residue = []
        atoms = []
        atom_cnt = 0
        offset = 0
        min_lengths = {'ALA':5,  'CYS':6,  'CYX':6,  'CYD':6, 'CYZ':6, 'ASP':8,  'GLU':9, 'PHE':11,
                       'GLY':4,  'HIS':10, 'HID':10, 'HIE':10,
                       'ILE':8,  'LYS':9,  'LEU':8,  'MET':8, 'ASN':8,  'PRO':7, 'GLN':9,
                       'ARG':11, 'SER':6,  'THR':7,  'VAL':7, 'TRP':14, 'TYR':12}
        for line in PDBlines:
            current_line = line[25:29]
            if current_line == last_line:
                # handle duplicate atom records
                atomname = line[12:16]
                splitname = string.split(string.strip(atomname))
                if len(splitname) == 2:          # have a 'CA A' or some such thing
                    if splitname[1] != 'A':      # dismiss B or C duplicates
                        last_line = current_line
                        atom_cnt += 1
                        continue
                    else:
                        line = line[:15] + ' ' + line[16:]    # just get rid of the 'A'
                        self.atoms[atom_cnt].atom_type = string.split(self.atoms[atom_cnt].atom_type)[0]
                else:
                    # if the line has a number, followed by a letter, it must be duplicate
                    if line[14] in ['1','2','3','4'] and line[15] in ['A','B','C','D']:
                        if line[15] in ['B','C','D']:             # get rid of the Bs Cs and Ds
                            last_line = current_line
                            atom_cnt += 1
                            continue
                        else:
                            line = line[:15] + ' ' + line[16:]
                            self.atoms[atom_cnt].atom_type = self.atoms[atom_cnt].atom_type[:-1]
                # just append to the previously created one
                current_residue.append(line)
                atoms.append(self.atoms[atom_cnt])
            else:
                # make a new amino acid and append it to residues list
                if len(current_residue)>0:
                    if len(current_residue) >= min_lengths[current_residue[0][17:20]]:
                        new_residue = AminoAcid(self, current_residue, atoms)
                        self.add_residue(new_residue)
                # reset the current list and append the new line to it
                current_residue = []
                atoms = []
                current_residue.append(line)
                atoms.append(self.atoms[atom_cnt])
            last_line = current_line
            atom_cnt += 1
        # append the last one created
        if len(current_residue) >= min_lengths[current_residue[0][17:20]]:
            new_residue = AminoAcid(self, current_residue, atoms)
            self.add_residue(new_residue)
        # this initializes proteolysis type, so that the chain knows if a calculation has been done
        self.proteolysis_type = 'none'
        self.locate_termini()
        self.residues_dict = {}
        for res in self.residues:
            self.residues_dict['%s'%(res.res_number)] = res
        #self.fill_pseudo_sidechains(0)
        #self.fill_neighbors_lists()
        #self.assign_ss_from_header()
        
    def get_sequence(self):                             # returns a list of characters
        my_sequence = ''
        for rez in self.residues:
            my_sequence = my_sequence + rez.res_type1
        return my_sequence

    def print_sequence(self):
        x = self.get_sequence()
        print x
                
    def fill_proteolysis_fragments(self, plys_type, misses_allowed=100000):
        """ fill self.fragment_list. missed_allowed tells how many missed cut sites
            are allowed. Default is a very large number.
        """
        ind = 0
        for pchain in self.parent.ProteinList:
            if pchain.chain_name == self.chain_name:
                chain_number = ind
            ind += 1
        self.proteolysis_type = plys_type
        if plys_type == 'trypsin':
            # first identify all of the cut sites
            specificity = ['R','K']
            cut_sites = [0]
            index = 0
            for rez in self.residues:
                if rez.res_type1 in specificity:
                    if index != len(self.residues)-1:
                        if self.residues[index+1].res_type1 != 'P':
                            cut_sites.append(index)
                    else:
                        cut_sites.append(index)
                index = index + 1
            cut_sites.append(len(self.residues))
        elif plys_type == 'V8E':
            # first identify all of the cut sites
            specificity = ['E']
            cut_sites = [0]
            index = 0
            for rez in self.residues:
                if rez.res_type1 in specificity:
                    if index != len(self.residues)-1:
                        cut_sites.append(index)
                index = index + 1
            cut_sites.append(len(self.residues))
        elif plys_type == 'V8' or plys_type == 'V8DE':
            # first identify all of the cut sites
            specificity = ['E', 'D']
            cut_sites = [0]
            index = 0
            for rez in self.residues:
                if rez.res_type1 in specificity:
                    if index != len(self.residues)-1:
                        cut_sites.append(index)
                index = index + 1
            cut_sites.append(len(self.residues))
        print "cut sites at",cut_sites
        # now use them to make a list of possible fragments
        self.fragment_list    = []          # an array of [weight,start,end]
        for site1_index in range(0,len(cut_sites)):
            for site2_index in range(site1_index, len(cut_sites)):
                if site1_index == site2_index or site2_index > site1_index + misses_allowed:
                    continue
                weight = 18.00       # start with the weight of water to cover the termini
                for rez_index in range(cut_sites[site1_index],cut_sites[site2_index]):
                    weight = weight + self.residues[rez_index].mw
                self.fragment_list.append([cut_sites[site1_index],cut_sites[site2_index],weight, chain_number])
        fragdict = {}
        for frag in self.fragment_list:
            fragdict[frag[2]] = frag
        keys = fragdict.keys()
        keys.sort()
        for key in keys:
            if key > 1000 and key < 10000:
                print key, fragdict[key]
            
    # before using this, use fill_proteolysis_fragments
    # outdated by get_nearest_proteolysis_fragment in MS
    def get_nearest_proteolysis_fragment(self, query_weight):
        qw = string.atof(query_weight)
        closest_index = 0
        closest_dist  = 1000000
        weight_index = 0
        for weight in self.frag_weight_list:
            if abs(weight - qw) < closest_dist:
                closest_index = weight_index
                closest_dist = abs(weight - qw)
            weight_index = weight_index + 1                
        return_stuff = [self.frag_start_list[closest_index],
                        self.frag_end_list[closest_index],
                        self.frag_weight_list[closest_index]]
        return return_stuff

    def fill_pseudo_sidechains(self, type=1):
        """ calculate a single coordinate to represent the sidechain.
            type 0 calculates 2 A out from the bisector of a triangle formed between
            three consecutive alpha carbons. If the residue terminates a structural fragment,
            it is assigned to the average of its non-H sidechain atoms, or to the central_pt's
            position if there are no sidechain atoms. type 1 calculates the average coordinate from
            all non-H sidechain atoms. 
        """
        if type == 0:
            for rez_index in range(len(self.residues)):
                if not self.residues[rez_index].has_central_pt:
                    continue
                if self.residues[rez_index].is_Nterm == 1 or self.residues[rez_index].is_Cterm == 1:
                    self.residues[rez_index].pseudo_sidechain = copy.deepcopy(self.residues[rez_index].central_pt)
                    continue
                a1 = self.residues[rez_index-1].central_pt
                a2 = self.residues[rez_index].central_pt
                a3 = self.residues[rez_index+1].central_pt
                l1 = self.residues[rez_index].central_pt.dist(self.residues[rez_index-1].central_pt)
                l2 = self.residues[rez_index].central_pt.dist(self.residues[rez_index+1].central_pt)
                if l1<l2:
                    x = ((l2/l1)*(a1.x-a2.x)) + a2.x
                    y = ((l2/l1)*(a1.y-a2.y)) + a2.y
                    z = ((l2/l1)*(a1.z-a2.z)) + a2.z
                    tp = Point(x,y,z)
                    x = (a3.x + tp.x)/2
                    y = (a3.y + tp.y)/2
                    z = (a3.z + tp.z)/2
                    mp = Point(x,y,z)
                else:
                    x = ((l1/l2)*(a3.x-a2.x)) + a2.x
                    y = ((l1/l2)*(a3.y-a2.y)) + a2.y
                    z = ((l1/l2)*(a3.z-a2.z)) + a2.z
                    tp = Point(x,y,z)
                    x = a1.x + tp.x / 2
                    y = a1.y + tp.y / 2
                    z = a1.z + tp.z / 2
                    mp = Point(x,y,z)
                l2_mp_dist = a2.dist(mp)
                lmp_B_dist = l2_mp_dist + 2.0
                x = ( (lmp_B_dist/l2_mp_dist) * ( a2.x-mp.x) ) + mp.x
                y = ( (lmp_B_dist/l2_mp_dist) * ( a2.y-mp.y) ) + mp.y
                z = ( (lmp_B_dist/l2_mp_dist) * ( a2.z-mp.z) ) + mp.z
                self.residues[rez_index].pseudo_sidechain = Point(x,y,z)
        elif type == 1:
            for rez in self.residues:
                stuff = rez.get_average_sidechain_position()
                rez.pseudo_sidechain = Point(stuff[0],stuff[1],stuff[2])


    def fill_neighbors_lists(self, qscore=0.35, dist_thresh=10.0, force_rewrite=0):
        # first see if the contacts file has been previously generated
        s = len(self.residues)
        contact_list = numpy.zeros([s,s])
        create_new = 0
        filename = self.parent.get_filename_by_extension('.ctc', self.chain_name)
        if force_rewrite:
            create_new = 1
        else:
            try:
                contact_file = open(filename)
            except IOError:
                create_new = 1
        if create_new:
            # initialize the contacts list
            # initialize a 2D array to hold inter-sidechain distances
            distance_list = numpy.zeros([s,s])
            for rex in range(len(self.residues)):
                for rex2 in range(rex,len(self.residues)):
                    if rex != rex2:
                        distance_list[rex][rex2] = self.residues[rex].pseudo_sidechain.dist(self.residues[rex2].pseudo_sidechain)
            print 'sorting distance list'
            # now sort the lists -- first create a new 2D array
            sorted_list = numpy.zeros([s,s])
            for rex in range(len(self.residues)):
                print '.',
                taken = numpy.zeros([s])
                for rex2 in range(len(self.residues)):
                    minDist = 100.0
                    if rex != rex2:
                        for rex3 in range(len(self.residues)):
                            if not taken[rex3] and rex3 != rex:
                                curLen = distance_list[rex][rex3]
                                if (( curLen != 0.00 ) and (curLen <= minDist)):
                                    minDist = curLen
                                    saveK = rex3
                        taken[saveK] = 1
                        sorted_list[rex][rex2] = saveK
            print
            # calculate shielding
            print 'calculating shielding'
            for rex in range(len(self.residues)):
                print '.',
                p1a = self.residues[rex].central_pt
                p1b = self.residues[rex].pseudo_sidechain
                for rex2 in range(len(self.residues)):
                    if distance_list[rex][sorted_list[rex][rex2]] > dist_thresh:
                        contact_list[rex][sorted_list[rex][rex2]] = 10.0
                        break
                    p2a = self.residues[sorted_list[rex][rex2]].central_pt
                    p2b = self.residues[sorted_list[rex][rex2]].pseudo_sidechain
                    dst = p1a.point_line_distance(p1b,p2b)
                    dst2 = p2a.point_line_distance(p1b,p2b)

                    shielding = s = s2 = 1
                    if dst >= 0.00:
                        s = 1-math.exp(-((dst*dst)/4.0))
                    if dst2 >= 0.00:
                        s2 = 1-math.exp(-((dst2*dst2)/4.0))
                    shielding = s;
                    shielding = shielding * s2

                    for rex3 in range(rex2):
                        if ((not self.residues[sorted_list[rex][rex3]]) or (rex3==rex)):
                            continue
                        dst = [0,0]
                        dst[0] = self.residues[sorted_list[rex][rex3]].central_pt.point_line_distance(p1b,p2b)
                        dst[1] = self.residues[sorted_list[rex][rex3]].pseudo_sidechain.point_line_distance(p1b,p2b)
                        s = s2 = 1
                        # calculate the shielding value
                        if dst[0] != -1:
                            s = 1-math.exp(-((dst[0]**2)/4.0))
                        if dst[1] != -1:
                            s2 = 1-math.exp(-((dst[1]**2)/4.0))
                        shielding = shielding * s
                        shielding = shielding * s2
                    # this is the q score
                    contact_list[rex][sorted_list[rex][rex2]] = 20.0 * shielding / (p1b.dist(p2b) + 1.0)
            # write the contacts to a file
            print 'done shielding'
            contact_file = open(filename, 'w')
            for rex in range(len(self.residues)):
                write_string = ""
                for rex2 in range(rex, len(self.residues)):
                    write_string = write_string + "%5.3f, "%(contact_list[rex][rex2])
                write_string = write_string + '\n'
                if write_string != '\n':
                    contact_file.write(write_string)
            contact_file.close()
        else:           # else read the contacts_file to fill the contact_list
            for rex in range(len(self.residues)):
                buffer = contact_file.readline()
                if len(buffer) == 0:
                    break
                tokens = string.split(buffer, ',')
                token_index = 0
                for rex2 in range(rex, len(self.residues)):
                    contact_list[rex][rex2] = 100 * (float(tokens[token_index]))
                    token_index = token_index + 1
            contact_file.close()
        # now fill the neighbors lists
        for rex in range(len(self.residues)):
            self.residues[rex].neighbors = []
            self.residues[rex].data['neighbors'] = []
        for rex in range(len(self.residues)):
            neighborCntr = 0
            for rex2 in range(len(self.residues)):
                if contact_list[rex,rex2] > qscore:
                    # compare the distance between beta carbons
                    inter_beta_dist = self.residues[rex].pseudo_sidechain.dist(self.residues[rex2].pseudo_sidechain)
                    if inter_beta_dist < dist_thresh:
                        if rex2 not in self.residues[rex].neighbors:
                            self.residues[rex].neighbors.append(rex2)
                            self.residues[rex].data['neighbors'].append(rex2)
                        if rex not in self.residues[rex2].neighbors:
                            self.residues[rex2].neighbors.append(rex)
                            self.residues[rex2].data['neighbors'].append(rex)
        # shielding attempts to normalize the number of neighbors into a
        # score similar to percent solvent accessability 
        lengths = []
        for rez in self.residues:
            lengths.append(len(rez.neighbors))
        lengths.sort()
        top_index = int(math.floor(len(self.residues)*0.90))
        bottom_index = int(math.ceil(len(self.residues)*0.10))
        print 'min %s max %s (%s %s)'%(lengths[bottom_index], lengths[top_index], bottom_index, top_index)
        for rez in self.residues:
            if len(rez.neighbors) > lengths[top_index]:
                # if the number of neighbors is > the top 90'th, set at top 90'th
                nb = lengths[top_index]
            elif len(rez.neighbors) < lengths[bottom_index]:
                # if the number of neighbors is < the bottom 10'th, set at the bottom 10'th
                nb = lengths[bottom_index]
            else:
                nb = len(rez.neighbors)
            rez.features['shielding'] = (nb-lengths[bottom_index])/float(lengths[top_index]-lengths[bottom_index])
        # and the average over all neighbors tells about the region
        for res in self.residues:
            sum = res.features['shielding']
            for rez in res.neighbors:
                sum += self.residues[rez].features['shielding']
            res.features['average_shielding'] = sum/(len(res.neighbors)+1.0)
        # calculate and print the average number of neighbors
        sum = 0.0
        for res in self.residues:
            sum += len(res.neighbors)
        sum /= len(self.residues)
        print '%5.3f neighbors per residue, on average'%(sum)

    def detect_domains(self):
        # make sure neighbors have been calculated
        for rez in self.residues:
            try:
                rez.features['shielding']
            except KeyError:
                print 'calculating neighbors'
                self.fill_pseudo_sidechains()
                self.fill_neighbors_lists()

    def assign_asa(self, solvent_radius, point_count, forced_rewrite=0):
        filename = self.parent.get_filename_by_extension('.bsa', self.chain_name)
        create_new = 0
        if forced_rewrite:
            create_new = 1
        else:
            try:
                asa_file = open(filename)
            except IOError:
                create_new = 1
        if create_new:
            # speedup so that areas are not calculated every time
            areas  = {'C':(4.0/3.0) * 3.141592654 * ((solvent_radius + 1.75)**3),
                      'O':(4.0/3.0) * 3.141592654 * ((solvent_radius + 1.4)**3),
                      'N':(4.0/3.0) * 3.141592654 * ((solvent_radius + 1.55)**3),
                      'S':(4.0/3.0) * 3.141592654 * ((solvent_radius + 1.8)**3),
                      'P':(4.0/3.0) * 3.141592654 * ((solvent_radius + 2.0)**3),
                      'H':(4.0/3.0) * 3.141592654 * ((solvent_radius + 1.17)**3),
                      'Z':(4.0/3.0) * 3.141592654 * ((solvent_radius + 3.0)**3)}
            default_area = (4.0/3.0) * 3.141592654 * ((solvent_radius + 1.8)**3)

            sphere_res = 15
            if self.x_table == None:
                self.build_futamura_intersection_table(solvent_radius)
            x_table = self.x_table
            # create spheres for each atom
            for res in self.residues:
                total_points         = 0
                last_total_points    = 0
                last_intra_inaccess  = 0
                total_side           = 0
                intra_inaccessible   = 0
                chain_intra          = 0
                chain_inter          = 0
                side_intra           = 0
                side_inter           = 0
                side_area            = 0.0
                total_area           = 0.0
                total_sidechain_area = 0.0
                for atom in res.atoms:
                    this_atoms_noncovalent_points = 0
                    radius = atom.radius + solvent_radius
                    radius_sq = radius**2
                    # figure out which ones to keep
                    for i in range(point_count):
                        # build the point
                        angle = random.random() * 2 * 3.141592654
                        z = (random.random() * 2 * radius) - radius
                        z_sq = z**2;
                        x_store = math.sqrt(radius_sq - z_sq) * math.cos(angle) + atom.x
                        y_store = math.sqrt(radius_sq - z_sq) * math.sin(angle) + atom.y
                        z_store = z + atom.z
                        externally_broken = 0
                        # see if the point is blocked by any other atoms
                        for second_atom in x_table['%s'%(atom.atom_number)]:
                            # if the point is within range of a second atom from the intersection table,
                            if math.sqrt(pow(x_store-second_atom[0],2) + pow(y_store-second_atom[1],2) + pow(z_store-second_atom[2],2)) <= (second_atom[3]):
                                if atom.atom_type == 'C' and second_atom[6] == 'N' and second_atom[4]-atom.res_number==1:
                                    break
                                elif atom.atom_type == 'N' and second_atom[6] == 'C' and atom.res_number-second_atom[4]==1:
                                    break
                                elif second_atom[4] == res.res_number and second_atom[5] == res.chain_name:
                                    break
                                else:
                                    externally_broken = 1
                        else:
                            if externally_broken:
                                this_atoms_noncovalent_points += 1
                                intra_inaccessible += 1
                                if atom.atom_type not in ['N', 'C', 'O']:
                                    total_side += 1
                                    side_intra += 1
                            else:
                                this_atoms_noncovalent_points += 1
                                if atom.atom_type not in ['N', 'C', 'O']:
                                    total_side += 1
                    if this_atoms_noncovalent_points > 0:
                        total_area += (this_atoms_noncovalent_points/float(point_count)) * areas.get(atom.atom_type[0], default_area)
                        if atom.atom_type not in ['N', 'C', 'O']:
                            total_sidechain_area += (this_atoms_noncovalent_points/float(point_count)) * areas.get(atom.atom_type[0], default_area)
                        #total_area += (this_atoms_noncovalent_points/float(point_count)) * (4.0/3.0) * 3.141592654 * ((atom.radius+solvent_radius)**3)
                    #print '%3s%3s total area: %5.1f exposed area: %5.1f percentage %5.2f'%(atom.atom_type,
                    #                                                                       atom.atom_number,
                    #                                                                       atoms_area,
                    #                                                                       this_area,
                    #                                                                       this_atoms_noncovalent_points/float(point_count))
                    total_points += this_atoms_noncovalent_points
                    last_total_points = this_atoms_noncovalent_points
                    last_intra_inaccess = intra_inaccessible
                            
                total_points = total_points + 0.0
                total_side   = total_side   + 0.0
                res.features['asa'] = (total_points-intra_inaccessible) / float(total_points)
                res.features['sidechain_asa'] = (total_side-side_intra) / float(total_side)
                res.data['exposed_area'] = total_area
                res.data['exposed_sidechain_area'] = total_sidechain_area
                print 'res %s%s - perc. asa - %5.2f, %5.2f side; exposed area - %5.2f, %5.2f side'%(res.res_number, res.res_type, res.features['asa'], res.features['sidechain_asa'], res.features['asa']*res.data['exposed_area'], res.features['sidechain_asa']*res.data['exposed_sidechain_area'])
                
            asa_file = open(filename, 'w')
            for rex in range(len(self.residues)):
                asa_file.write("%s %5.3f %5.3f %5.3f %5.3f\n"%(self.residues[rex].res_number, self.residues[rex].features['asa'], self.residues[rex].features['sidechain_asa'], self.residues[rex].data['exposed_area'], self.residues[rex].data['exposed_sidechain_area']))
            asa_file.close()
        else:           # else read the contacts_file to fill the contact_list
            for rex in range(len(self.residues)):
                buffer = asa_file.readline()
                if len(buffer) == 0:
                    break
                tokens = string.split(buffer)
                self.residue_dict[string.atoi(tokens[0])].features['asa'] = string.atof(tokens[1])
                self.residue_dict[string.atoi(tokens[0])].features['sidechain_asa'] = string.atof(tokens[2])
                self.residue_dict[string.atoi(tokens[0])].data['exposed_area'] = string.atof(tokens[3])
                self.residue_dict[string.atoi(tokens[0])].data['exposed_sidechain_area'] = string.atof(tokens[4])
            asa_file.close()

    def get_core_alpha_carbons(self, core_cutoff=8, neighbor_thresh=3):
        first_pass        = []
        central_atom_list = []
        central_atom_nums = []
        for rez in self.residues:
            if rez.has_central_pt==1 and len(rez.neighbors) >= core_cutoff:
                first_pass.append(rez.res_number)
        for rez in self.residues:
            if rez.has_central_pt==1 and len(rez.neighbors) >= core_cutoff:
                core_neighbors = 0
                for neighbor in rez.neighbors:
                    if self.residues[neighbor].res_number in first_pass:
                        core_neighbors += 1
                if core_neighbors >= neighbor_thresh:
                    central_atom_list.append(rez.central_atom)
                    central_atom_nums.append(rez.res_number)
        for rez in self.residues:
            if rez.res_number not in central_atom_nums:
                if (rez.res_number + 1) in central_atom_nums and (rez.res_number - 1) in central_atom_nums:
                    central_atom_list.append(rez.central_atom)
        return central_atom_list

    def assign_core_atom_distances(self, core_cutoff, neighbor_thresh=3, mindist=0.0, maxdist=12.0):
        f_hash = FutamuraHash(self, 0.0, maxdist, 2, core_cutoff, neighbor_thresh)
        atom_list = self.get_core_alpha_carbons(core_cutoff, neighbor_thresh)
        for cp in atom_list:
            cp.data['nearby_cps'] = f_hash.get_atoms_between(cp.atom_number, mindist, maxdist)

    def fill_secondary(self):
        self.assign_ss()
    
    def assign_ss(self):
        """ assign_ss first looks for header information on secondary structure,
            then if it fails it assigns ss using the alpha_carbon_distances.
        """
        self.assign_ss_from_header()

    def assign_ss_from_alpha_carbon_distances(self):
        averages = {'a1a3s':{'A':5.6,'B':6.6},
                    'a1a4s':{'A':5.8,'B':9.6},
                    'a1a5s':{'A':6.8,'B':12.6},
                    'a2a4s':{'A':5.6,'B':6.6},
                    'a2a5s':{'A':5.6,'B':9.8},
                    'a3a5s':{'A':5.5,'B':6.6}}

        rmsthresh = 3.5
        for r_ind in range(len(self.residues)):
            self.residues[r_ind].features['secondary'] = 'C'
            if self.residues[r_ind].is_Nterm or self.residues[r_ind].is_Cterm:
                    continue
            if self.residues[r_ind-1].is_Nterm or self.residues[r_ind-1].is_Cterm:
                    continue
            if self.residues[r_ind+1].is_Nterm or self.residues[r_ind+1].is_Cterm:
                    continue
            rmsc = 0.0
            rmsc += (averages['a1a3s']['A'] - self.residues[r_ind-2].central_atom.dist(self.residues[r_ind+0].central_atom))**2
            rmsc += (averages['a1a4s']['A'] - self.residues[r_ind-2].central_atom.dist(self.residues[r_ind+1].central_atom))**2
            rmsc += (averages['a1a5s']['A'] - self.residues[r_ind-2].central_atom.dist(self.residues[r_ind+2].central_atom))**2
            rmsc += (averages['a2a4s']['A'] - self.residues[r_ind-1].central_atom.dist(self.residues[r_ind+1].central_atom))**2
            rmsc += (averages['a2a5s']['A'] - self.residues[r_ind-1].central_atom.dist(self.residues[r_ind+2].central_atom))**2
            rmsc += (averages['a3a5s']['A'] - self.residues[r_ind+0].central_atom.dist(self.residues[r_ind+2].central_atom))**2
            rmscA = rmsc / math.sqrt(rmsc/5.0)

            rmsc = 0.0
            rmsc += (averages['a1a3s']['B'] - self.residues[r_ind-2].central_atom.dist(self.residues[r_ind+0].central_atom))**2
            rmsc += (averages['a1a4s']['B'] - self.residues[r_ind-2].central_atom.dist(self.residues[r_ind+1].central_atom))**2
            rmsc += (averages['a1a5s']['B'] - self.residues[r_ind-2].central_atom.dist(self.residues[r_ind+2].central_atom))**2
            rmsc += (averages['a2a4s']['B'] - self.residues[r_ind-1].central_atom.dist(self.residues[r_ind+1].central_atom))**2
            rmsc += (averages['a2a5s']['B'] - self.residues[r_ind-1].central_atom.dist(self.residues[r_ind+2].central_atom))**2
            rmsc += (averages['a3a5s']['B'] - self.residues[r_ind+0].central_atom.dist(self.residues[r_ind+2].central_atom))**2
            rmscB = rmsc / math.sqrt(rmsc/5.0)
            if rmscA > rmsthresh and rmscB > rmsthresh:
                pass
            elif rmscB < rmscA:
                self.residues[r_ind].features['secondary'] = 'B'
            elif rmscA < rmscB:
                self.residues[r_ind].features['secondary'] = 'A'

        
    def assign_ss_from_header(self):
        for res in self.residues:
            res.features['secondary'] = 'C'
        found_secondary = 0
        for line in self.parent.HeaderLines:
            if line[:5] == 'HELIX':
                found_secondary = 1
                if string.strip(line[19:20]) == string.strip(self.chain_name):
                    start_res = int(line[21:25])
                    end_res = int(line[33:37])
                    for i in range(start_res, end_res):
                        try: 
                            self.residue_dict[i]
                        except KeyError:
                            pass
                        else:
                            self.residue_dict[i].features['secondary'] = 'A'
            if line[:5] == 'SHEET':
                found_secondary = 1
                if string.strip(line[21:22]) == string.strip(self.chain_name):
                    start_res = int(line[22:26])
                    end_res = int(line[33:37])
                    for i in range(start_res, end_res):
                        try:
                            self.residue_dict[i]
                        except KeyError:
                            pass
                        else:
                            self.residue_dict[i].features['secondary'] = 'B'
        if found_secondary == 0:
            self.assign_ss_from_alpha_carbon_distances()

    def get_hydrogens (self,chain_name, res_number, atom_name):
        """ Return a list of hydrogens that are associated with an atom """
        hydAtoms=[]
        hydrogen_pat=re.compile("^d*H")
        target_res = self.residues_dict[str(res_number)]
        # Loop through all of the atoms and search for the hydrogens
        for atom in target_res.atoms:
            # If it is a hydrogen atom then find its closest non hydrogen
            # atom and that atom is its owner. If the owner is atom_name
            # then add it to the list of hydAtoms
            if (hydrogen_pat.search(atom.atom_type)):
                hyd=atom
                minDist=-1
                owner=None
                for atom2 in target_res.atoms:
                    if (not hydrogen_pat.search(atom2.atom_type) and
                        (hyd.dist(atom2) < minDist or minDist == -1)):
                        owner = atom2
                        minDist = hyd.dist(atom2)
                # If the closest atom is the atom_name then add it to the hydrogen
                # list
                if (owner.atom_type == target_res.atoms_dict[atom_name].atom_type):
                    hydAtoms.append (hyd)
        return hydAtoms

    def find_atoms_for_protonation (self,protonsInfo):
        """ Find all of the atoms that can be protonated in this protein.
            Returns a list of protons where an element contains
            {'atom','aa'[,'prev_aa'],'protonInfo'}. """
        protons = []
        for j in range(len(self.residues)):
            aa = self.residues[j]
        
            # Add all of the amino acid specific protons
            for i in range(len(aa.atoms)):
                # Find the proton information for this atom if it has any
                protonInfo = find_proton_info (aa.res_type,aa.atoms[i],protonsInfo)
                if (protonInfo != None):
                    protons.append ({'atom': aa.atoms[i], 'aa': aa, 'protonInfo': protonInfo})
  
            # if it is the first amino acid in the chain (N-TERMINUS)
            if (j == 0):
                protonInfo = find_proton_info ('N-TERMINUS',aa.atoms_dict['N'],protonsInfo)
                if (protonInfo != None):
                    protons.append ({'atom': aa.atoms_dict['N'], 'aa': aa, 'protonInfo': protonInfo})
            # Add the backbone amino acid, which is common to all amino acids except
            # the first amino acid in a chain
            elif (aa.res_type != 'PRO'): # Every other residue except PRO
                protonInfo = find_proton_info ('BACKBONE',aa.atoms_dict['N'],protonsInfo)
                if (protonInfo != None):
                    # Store the previous amino acid because we need its carbon
                    protons.append ({'atom': aa.atoms_dict['N'], 'aa': aa, 'prev_aa': self.residues[j-1],
                                     'protonInfo': protonInfo})
        return protons

    def protonate (self,protonFile,redo=False):
        """ Protonate a protein. Add the hydrogens to the protein. i.e. create a list
            of possible donors and all of their hydrogens 'atomsHydList': {'atom','hydAtoms'} """
        self.atomsHydList=[]

        protonsInfo = read_protonation_info (protonFile)
        protons = self.find_atoms_for_protonation (protonsInfo)

        # Position the hydrogens
        for proton in protons:
            D=proton['atom']
            aa=proton['aa']
            hydAtoms=[]
            # Initialize the hydrogen atoms to any that are already in the protein
            hydAtoms = self.get_hydrogens(D.chain_name,D.res_number,D.atom_type)
            if (len(hydAtoms) != 0 and redo == False):
                # Save the list of hyrdogen atoms for this atom 
                self.atomsHydList.append ({'atom': D, 'hydAtoms': hydAtoms})
                continue
            elif (len(hydAtoms) != 0 and redo == True):
                print ("TODO: Add deletion of the current hydrogens")
                continue
      
            # sp2, 1H 2DD
            if (proton['protonInfo']['hyb'] == 'sp2' and
                proton['protonInfo']['bonds'] == '1H 2DD'):
                DD1Name=proton['protonInfo']['DD1Name']
                angle_offset=0 # 0 -> DD1-D-H = DD2-D-H
                if (DD1Name == 'None'):
                    print "ERROR: No DD1 atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD1=aa.atoms_dict[DD1Name]
                if (proton['protonInfo']['angles'] == 'DD1-D-H = DD2-D-H'):
                    DD2Name=proton['protonInfo']['DD2Name']
                    if (DD2Name=='None'):
                        print "ERROR: No DD2 atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                        continue
                    DD2=aa.atoms_dict[DD2Name]
                elif (proton['protonInfo']['angles'] == '(C-N-H)-(CA-N-H)=4; C CA N H are planar'):
                    prev_aa=proton['prev_aa']
                    DD2Name=proton['protonInfo']['DD2Name']
                    if (DD2Name=='None'):
                        print "ERROR: No DD2 atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                        continue
                    DD2=prev_aa.atoms_dict[DD2Name]
                    angle_offset = 4
                # Find the angle C-N-H to place the hydrogen
                Dv=r_[D.x,D.y,D.z]
                # Make Dv the origin
                DD2v=r_[DD2.x,DD2.y,DD2.z]-Dv
                DD1v=r_[DD1.x,DD1.y,DD1.z]-Dv
                Dv=r_[0,0,0]
                theta=acos(dot(DD1v,DD2v)/(mag(DD1v)*mag(DD2v)))
                angle=pi-theta/2+radians(angle_offset)
                hydPos=findPlanarPosition (proton['protonInfo']['D-H'],angle*180/pi,D,DD2,DD1,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos))
          
            # sp2, 1H 1DD
            elif (proton['protonInfo']['hyb'] == 'sp2' and
                  proton['protonInfo']['bonds'] == '1H 1DD'):
                DDName=proton['protonInfo']['DD1Name']
                if (DDName == 'None'):
                    print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD=aa.atoms_dict[DDName]
                DDDName=proton['protonInfo']['DDD1Name']
                if (DDDName=='None'):
                    print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DDD=aa.atoms_dict[DDDName]
                # This configuration has two mutually exclusive hydrogen positions
                # TODO: For now we are placing both hydrogens but we need to add some kind of
                # collision resolution
                # TODO: Change find planar to take an array of angles...
                hydPos1=findPlanarPosition (proton['protonInfo']['D-H'],250,D,DD,DDD,True)
                hydPos2=findPlanarPosition (proton['protonInfo']['D-H'],110,D,DD,DDD,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1))
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos2))

            # sp2, 2H 1DD
            elif (proton['protonInfo']['hyb'] == 'sp2' and
                  proton['protonInfo']['bonds'] == '2H 1DD'):
                DDName=proton['protonInfo']['DD1Name']
                if (DDName == 'None'):
                    print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD=aa.atoms_dict[DDName]
                DDDName=proton['protonInfo']['DDD1Name']
                if (DDDName=='None'):
                    print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DDD=aa.atoms_dict[DDDName]
                # This configuration has two hydrogen positions
                # TODO: Change find planar to take an array of angles...
                hydPos1=findPlanarPosition (proton['protonInfo']['D-H'],120,D,DD,DDD,True)
                hydPos2=findPlanarPosition (proton['protonInfo']['D-H'],360-120,D,DD,DDD,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1)) 
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos2))

            # sp3, 1H 1DD
            elif (proton['protonInfo']['hyb'] == 'sp3' and
                  proton['protonInfo']['bonds'] == '1H 1DD'):
                DDName=proton['protonInfo']['DD1Name']
                if (DDName == 'None'):
                    print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD=aa.atoms_dict[DDName]
                DDDName=proton['protonInfo']['DDD1Name']
                if (DDDName=='None'):
                    print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DDD=aa.atoms_dict[DDDName]
                if (proton['protonInfo']['angles'] == 'DD-D-H=110'):
                    hydPos1=findCirclePosition (proton['protonInfo']['D-H'],110,0,D,DD,DDD,True)
                elif (proton['protonInfo']['angles'] == 'DD-D-H=96'):
                    hydPos1=findCirclePosition (proton['protonInfo']['D-H'],96,0,D,DD,DDD,True)
                elif (proton['protonInfo']['angles'] == 'DD-D-H=110; Circle Angle=240'):
                    hydPos1=findCirclePosition (proton['protonInfo']['D-H'],110,240,D,DD,DDD,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1))

            # sp3, 3H 1DD
            elif (proton['protonInfo']['hyb'] == 'sp3' and
                  proton['protonInfo']['bonds'] == '3H 1DD'):
                DDName=proton['protonInfo']['DD1Name']
                if (DDName == 'None'):
                    print "ERROR: No DD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DD=aa.atoms_dict[DDName]
                DDDName=proton['protonInfo']['DDD1Name']
                #print 'proton %s'%(proton['protonInfo'])
                # special case: N-term glycine doesn't really have a DDDName so continue silently
                if proton['protonInfo']['aminoAcidName'] == 'N-TERMINUS' and aa.res_type1 == 'G':
                    continue
                if (DDDName=='None'):
                    print "ERROR: No DDD atom for proton atom# %d, res# %d"%(D.atom_number,D.res_number)
                    continue
                DDD=aa.atoms_dict[DDDName]
                hydPos1=findCirclePosition (proton['protonInfo']['D-H'],110,0,D,DD,DDD,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos1))
                hydPos2=findCirclePosition (proton['protonInfo']['D-H'],110,120,D,DD,DDD,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos2))
                hydPos3=findCirclePosition (proton['protonInfo']['D-H'],110,240,D,DD,DDD,True)
                hydAtoms.append(self.parent.create_hydrogen (D.chain_name,D.res_number,D.atom_type,hydPos3))

            # Save the list of hyrdogen atoms for this atom
            self.atomsHydList.append ({'atom': D, 'hydAtoms': hydAtoms})

        self.parent.renumber_atoms()

    def get_avail_donors (self,donorFile):
        donorsInfo=read_donor_info (donorFile)

        # Find all of the atoms that have hydrogens and that can donate
        availDonors=[]
        number_list = []
        for atomHyd in self.atomsHydList:
            donorInfo = find_donor_info (atomHyd['atom'].parent.res_type,atomHyd['atom'],donorsInfo)
            if (donorInfo != None):
                availDonors.append ({'donorAtom': atomHyd['atom'], 'hydAtoms': atomHyd['hydAtoms']})

        return availDonors

    def get_avail_acceptors (self,accFile):
        accsInfo=read_acc_info(accFile)

         # Find all of the atoms in the protein that can be acceptors
        availAcceptors = []
        for j in range(len(self.residues)):
            aa = self.residues[j]
            for i in range(len(aa.atoms)):
                # Find the acceptor information if it has any
                accInfo = find_acc_info (aa.res_type,aa.atoms[i],accsInfo)
                if (accInfo != None):
                    availAcceptors.append ({'accAtom': aa.atoms[i],#'accInfo': accInfo, Don't need this?
                                            'AAAtom': aa.atoms_dict[accInfo['AAName']]})
        return availAcceptors

    def find_as_neighborhood (self,as_residue_inxs):
        """ Find the neighborhood of amino acids that has the largest number
            of active site residues
            Return: a list of residues
        """
        maxInx=-1
        maxCnt=-1
        # Find the neighborhood with the maximum number of as residues
        for i in range(len(as_residue_inxs)):
            as_rex = as_residue_inxs[i]
            cnt=0
            for rex in range(len(self.residues[as_rex].neighbors)):
                if (rex in as_residue_inxs): # if this neighbor is an as residue
                    cnt+=1
            if (maxCnt < cnt):
                maxInx=i
                maxCnt=cnt

        if (maxInx < 0):
            print "ERROR: incorrect index (find_as)"

        # Create the list of residues based on the neighborhood selected
        as_residues = [self.residues[as_residue_inxs[maxInx]]]
        for rex in self.residues[as_residue_inxs[maxInx]].neighbors:
            if (rex in as_residue_inxs):
                as_residues.append(self.residues[rex])
        return as_residues

    def find_as_concat (self,as_residue_inxs):
        """ Find the active site of the protein by grouping residues based on whether
            a residue can see at least one of the other residues in the group
            Return: a list of residues
        """
        graph=[]
        as_residues=[]
        # Create the graph representation of the neighborhoods
        for i in range(len(as_residue_inxs)):
            rex_row=as_residue_inxs[i]
            row = []
            for j in range(len(as_residue_inxs)):
                rex_col=as_residue_inxs[j]
                if (j == i): # Of course we can see ourselves
                    row.append (1)
                # Can see each other
                elif (rex_col in self.residues[rex_row].neighbors):
                    row.append(1)
                else: # Cannot see each other
                    row.append(0)
            graph.append(row)

        """
        # Test case
        # Should have three groups:
        # [0, 4, 3, 6, 1]
        # [2, 5, 9]
        # [7, 8]
        graph=[[1,0,0,0,1,0,0,0,0,0],
               [0,1,0,0,0,0,1,0,0,0],
               [0,0,1,0,0,1,0,0,0,0],
               [0,0,0,1,1,0,0,0,0,0],
               [1,0,0,1,1,0,1,0,0,0],
               [0,0,1,0,0,1,0,0,0,1],
               [0,1,0,0,1,0,1,0,0,0],
               [0,0,0,0,0,0,0,1,1,0],
               [0,0,0,0,0,0,0,1,1,0],
               [0,0,0,0,0,1,0,0,0,1]]
        """

        # Show the graph
        print 'Initial Graph'
        for row in graph:
            print row

        # Keep track of which residues have been attached to
        # a group
        attached=[]
        for i in range(len(graph)):
            attached.append(0)
            
        # Perform a BFT to get the possible active sites
        groups = []
        # While there is still a residue that is not attached
        while (0 in attached):
            # Find a residue that has not been attached
            for z in range(len(attached)):
                if (attached[z]==0):
                    break
            group = [z]
            attached[z]=1
            stack=[]
            stack.append(z)

            # Traverse the graph
            while (len(stack) > 0):
                row = graph[stack.pop()]
                for j in range(len(row)):
                    # Connected the current node and not already attached
                    # to a group
                    if (row[j] == 1 and attached[j]==0): 
                        group.append(j)
                        attached[j]=1
                        stack.append(j)
            groups.append (group)

        # Print the resultant groups
        print 'Final Groups'
        for group in groups:
            print group
            
        return as_residues


    def find_as_concat_n (self,n,as_residue_inxs):
        """ Find the active site of the protein by grouping residues based on whether
            a residue can see at least 'n' of the other residues in the group
            Return: a list of residues
        """
        graph=[]
        as_residues=[]
        # Create the graph representation of the neighborhoods
        for i in range(len(as_residue_inxs)):
            rex_row=as_residue_inxs[i]
            row = []
            for j in range(len(as_residue_inxs)):
                rex_col=as_residue_inxs[j]
                if (j == i): # Of course we can see ourselves
                    row.append (1)
                # Can see each other
                elif (rex_col in self.residues[rex_row].neighbors):
                    row.append(1)
                else: # Cannot see each other
                    row.append(0)
            graph.append(row)

        """
        # Test case
        # Should have three groups:
        # [0, 4, 3, 6, 1]
        # [2, 5, 9]
        # [7, 8]
        graph=[[1,0,0,0,1,0,0,0,0,0],
               [0,1,0,0,0,0,1,0,0,0],
               [0,0,1,0,0,1,0,0,0,0],
               [0,0,0,1,1,0,0,0,0,0],
               [1,0,0,1,1,0,1,0,0,0],
               [0,0,1,0,0,1,0,0,0,1],
               [0,1,0,0,1,0,1,0,0,0],
               [0,0,0,0,0,0,0,1,1,0],
               [0,0,0,0,0,0,0,1,1,0],
               [0,0,0,0,0,1,0,0,0,1]]
        """

        # Show the graph
        print 'Initial Graph'
        for row in graph:
            print row

        # Keep track of which residues have been attached to
        # a group
        attached=[]
        for i in range(len(graph)):
            attached.append(0)
            
        # Perform a BFT to get the possible active sites
        groups = []
        # While there is still a residue that is not attached
        while (0 in attached):
            # Find a residue that has not been attached
            for head in range(len(attached)):
                if (attached[head]==0):
                    break
            group = [head]
            attached[head]=1
            stack=[]
            stack.append(head)

            # Traverse the graph
            while (len(stack) > 0):
                row = graph[stack.pop()]
                for j in range(len(row)):
                    # Connected to the current node and not already attached
                    # to a group
                    if (row[j] == 1 and attached[j]==0): 
                        group.append(j)
                        attached[j]=1
                        stack.append(j)
            groups.append (group)

        # Find out if each amino acid in the group can see at least two
        # of the other amino acids by examining the neighborhoods
        for i in range(len(groups)):
            group=groups[i]
            # For each member of the group count the number it can
            # see
            group_tmp=[]
            for inx in group:
                residue=self.residues[as_residue_inxs[inx]]
                cnt=0
                # Search through all of the other members of the group
                for inx2 in group:
                    if (inx == inx2): continue
                    # If it is the neighborhood then increment the count
                    if (as_residue_inxs[inx2] in residue.neighbors):
                        cnt+=1
                if (cnt >= n):
                    group_tmp.append(inx)
            print 'Before:',group
            print 'After: ',group_tmp
            groups[i]=group_tmp
                    

        # Print the resultant groups
        print 'Final Groups'
        for group in groups:
            print group

        # Find the group with the largest number of residues = active site
        maxInx=-1
        maxCnt=-1
        for i in range(len(groups)):
            group = groups[i]
            if (len(group) > maxCnt):
                maxCnt=len(group)
                maxInx=i

        # Using that group get the actual residue numbers
        if (maxInx < 0):
            print "ERROR: incorrect index (find_as)"

        # Create the list of residues based on the group selected
        as_residues = []
        for inx in groups[maxInx]:
            as_residues.append(self.residues[as_residue_inxs[inx]])

        return as_residues

    def fill_densities(self, distance_cutoff=9.0, force_rewrite=0):
        Molecule.fill_densities(self, distance_cutoff, force_rewrite)
        self.fill_pseudo_sidechains(0)
        self.fill_neighbors_lists()
        for res in self.residues:
            sum = 0.0
            for atom in res.atoms:
                sum += atom.features['atomic_density']
            res.features['atomic_density'] = sum/float(len(res.atoms))
        for res in self.residues:
            sum = res.features['atomic_density']
            for rex in res.neighbors:
                rez = self.residues[rex]
                sum += rez.features['atomic_density']
            res.features['average_atomic_density'] = sum / (len(res.neighbors)+1.0)

            
        

