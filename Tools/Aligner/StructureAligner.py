import Numeric
import math
import copy

import sys
sys.path.append('./Tools/Superimpose')
import Superimpose.Superimpose
sys.path.append('./Tools/Aligner')
import Aligner

class StructureAligner(Aligner.Aligner):
    def __init__(self, print_param=0, alignment_type='local', gap_scores=[0.8,0.5]):
        Aligner.Aligner.__init__(self, print_param, alignment_type, gap_scores)

    def _prepare_polymer(self, polymer):
        """ Add StructureAlignerDistances to the CAs .features dictionaries.
            Returns a central atom list.
        """
        polymer.fill_pseudo_sidechains(1)
        polymer.fill_neighbors_lists(0.50,10.0,0)
        res_ind = 0
        for res in polymer.residues:
            res_ind += 1
        
        cind = 0
        cas = polymer.get_central_atom_list()
        terms = []
        for ca in cas:
            if ca.parent.is_Nterm or ca.parent.is_Cterm:
                terms.append(cind)
            cind += 1
        cind = 0
        for ca in cas:
            if cind not in terms and cind-1 not in terms and cind+1 not in terms and cind-2 not in terms and cind+2 not in terms:
                ax0,ay0,az0 = cas[cind-2].x, cas[cind-2].y, cas[cind-2].z
                ax1,ay1,az1 = cas[cind-1].x, cas[cind-1].y, cas[cind-1].z
                ax2,ay2,az2 = cas[cind+0].x, cas[cind+0].y, cas[cind+0].z
                ax3,ay3,az3 = cas[cind+1].x, cas[cind+1].y, cas[cind+1].z
                ax4,ay4,az4 = cas[cind+2].x, cas[cind+2].y, cas[cind+2].z
                a0_a2 = math.sqrt(((ax2-ax0)**2) + ((ay2-ay0)**2) + ((az2-az0)**2))
                a0_a3 = math.sqrt(((ax3-ax0)**2) + ((ay3-ay0)**2) + ((az3-az0)**2))
                a0_a4 = math.sqrt(((ax4-ax0)**2) + ((ay4-ay0)**2) + ((az4-az0)**2))
                a1_a3 = math.sqrt(((ax3-ax1)**2) + ((ay3-ay1)**2) + ((az3-az1)**2))
                a1_a4 = math.sqrt(((ax4-ax1)**2) + ((ay4-ay1)**2) + ((az4-az1)**2))
                a2_a4 = math.sqrt(((ax4-ax2)**2) + ((ay4-ay2)**2) + ((az4-az2)**2))
                polymer.residues[cind].features['StructureAlignerDistances'] = [a0_a2,a0_a3,a0_a4,a1_a3,a1_a4,a2_a4]
            cind += 1
        return polymer
    
    def score(self, r, m, w=[0.0,0.11,0.94,0.0]):
        rres = self.target.residues
        mres = self.template.residues
        rlen = len(rres)
        mlen = len(mres)
        #if rres[r-2].is_Nterm or rres[r-2].is_Cterm or rres[r-1].is_Nterm or rres[r-1].is_Cterm or rres[r].is_Nterm or rres[r].is_Cterm or rres[r+1].is_Nterm or rres[r+1].is_Cterm or rres[r+2].is_Nterm or rres[r+2].is_Cterm: 
            #print '%s%s %s%s 0.8'%(self.target.residues[r].res_type1, self.target.residues[r].res_number,
            #                      self.template.residues[m].res_type1, self.template.residues[m].res_number) 
        try:
            self.target.residues[r].features['StructureAlignerDistances']
        except KeyError:
            return 0.0
        try:
            self.template.residues[m].features['StructureAlignerDistances']
        except KeyError:
            return 0.0
        d1 = self.target.residues[r].features['StructureAlignerDistances']
        d2 = self.template.residues[m].features['StructureAlignerDistances']
        s1 = self.target.residues[r].features['shielding']
        s2 = self.template.residues[m].features['shielding']
        #print r,m,s1,s2
        types = ['difference', 'rms', 'shielding', 'to_neighbors_difference', 'mix']
        run_type = types[4]
        if run_type == 'difference':
            similarity = 0.0
            for i in range(len(d1)):
                if abs(d1[i]-d2[i]) < 1.0:
                    similarity += 1.0 - abs(d1[i]-d2[i])
            score = w[0]*similarity/len(d1)
        elif run_type == 'rms':
            distance = 0.0001
            for i in range(len(d1)):
                distance += (d1[i]-d2[i])**2
            score = 1/math.sqrt(distance/len(d1))
        elif run_type == 'shielding':
            similarity = 1.0 - abs(s1-s2)
            score = w[2]*similarity/len(d1)
        elif run_type == 'to_neighbors_difference':
            vector_list = [[], []]  # vector list; one for target, one for template
            neighbors = [self.target.residues[r].neighbors, self.template.residues[r].neighbors]
            protein = [self.target, self.template]
            central = [self.target.residues[r].central_atom, self.template.residues[r].central_atom]
            if len(neighbors[0]) > 1 and len(neighbors[1]) > 1:
                for i in range(len(vector_list)):       # iterating through residues saves code rewriting
                    for t1 in range(0, len(neighbors[i])):
                        ac1 = protein[i].residues[neighbors[i][t1]].central_atom
                        vector_list[i].append([ac1.dist(central[i])])
            score = w[3]*self.__get_distance_score(vector_list)
        elif run_type == 'mix':
            difference_similarity = 0.0
            w[0] = 0.0
            if w[0] != 0.0:
                # difference
                for i in range(len(d1)):
                    if abs(d1[i]-d2[i]) < 1.0:
                        difference_similarity += 1.0 - abs(d1[i]-d2[i])
                difference_similarity /= len(d1)
            rms_similarity = 0.0
            if w[1] != 0.0:
                # inverted rms
                distance = 0.0001
                for i in range(len(d1)):
                    distance += (d1[i]-d2[i])**2
                rms_similarity = 1/math.sqrt(distance/len(d1))
            shielding_similarity = 0.0
            if w[2] != 0.0:
                similarity = 1.0 - abs(s1-s2)
                shielding_similarity = similarity/len(d1)
            to_neighbors_similarity = 0.0
            w[3] = 0.0
            if w[3] != 0.0:
                # to neighbors distance
                vector_list = [[], []]  # vector list; one for target, one for template
                neighbors = [self.target.residues[r].neighbors, self.template.residues[r].neighbors]
                protein = [self.target, self.template]
                central = [self.target.residues[r].central_atom, self.template.residues[r].central_atom]
                if len(neighbors[0]) > 3 and len(neighbors[1]) > 3:
                    for i in range(len(vector_list)):       # iterating through residues saves code rewriting
                        for t1 in range(0, len(neighbors[i])):
                            ac1 = protein[i].residues[neighbors[i][t1]].central_atom
                            vector_list[i].append([ac1.dist(central[i])])
                to_neighbors_similarity = self.__get_distance_score(vector_list)
            wt_sum = w[0] + w[1] + w[2] + w[3]
            w0 = w[0]/wt_sum
            w1 = w[1]/wt_sum
            w2 = w[2]/wt_sum
            w3 = w[3]/wt_sum
            score = (w0*difference_similarity + w1*rms_similarity + w2*shielding_similarity + w3*to_neighbors_similarity)
        return score
                                
    def __get_distance_score(self, vector_list):
            success = 0.0
            fail    = 0.001
            minlen = min(len(vector_list[0]),len(vector_list[1]))
            if minlen == 0:
                return 0.0
            if len(vector_list[0]) < 100 and len(vector_list[1]) < 100:
                changed = 1
                while len(vector_list[0])>=0 and len(vector_list[1])>=0 and changed==1:
                    changed = 0
                    v1_ind = 0
                    for v1 in vector_list[0]:
                        v2_ind = 0
                        for v2 in vector_list[1]:
                            avg_dist = 0.0
                            for d_ind in range(len(v1)):
                                avg_dist += abs(v1[d_ind]-v2[d_ind])
                            avg_dist /= len(v1)
                            reset = 0
                            if avg_dist < 0.5:
                                success += 1
                                del vector_list[0][v1_ind]
                                del vector_list[1][v2_ind]
                                changed = 1
                                break
                            v2_ind += 1
                        else:
                            fail += 1
                        v1_ind += 1
                        if changed == 1:
                            break
            return success / minlen
        
    def align_structures(self):
        aout, bout = self.align()
        list1 = []
        list2 = []
        for ind in range(len(aout)):
            if aout[ind] != '-' and bout[ind] != '-':
                aout[ind].vtk_arg_list['trace']['color'] = [1.0,1.0,1.0]
                bout[ind].vtk_arg_list['trace']['color'] = [1.0,1.0,1.0]
                list1.append(aout[ind].central_atom)
                list2.append(bout[ind].central_atom)
        ps = [self.target.atoms, self.template.atoms]
        tm, rm = Superimpose.Superimpose._fit_pair_atoms([list1,list2])
        scores = 0.0
        cntr  = 0
        for i in range(len(aout)):
            if aout[i] != '-' and bout[i] != '-':
                scores += self.score_by_key(aout[i].res_number, bout[i].res_number, [0.0,0.11,0.91,0.0])
                cntr += 1
        print 'avg score %s'%(scores/cntr)
        Superimpose.Superimpose._transform(tm, rm, self.template.atoms)
        rms = Superimpose.Superimpose._get_atoms_rms([list1, list2])
        print 'initial rms = %s over %s atoms'%(rms, len(list1))

        for max_dist_for_pairing in [5.0,4.0,3.0]:
            target_central_atoms   = copy.deepcopy(list1)
            template_central_atoms = copy.deepcopy(list2)
            new_target_atoms       = []
            new_template_atoms     = []
            for i in range(len(target_central_atoms)):
                #print '%s %s %s'%(target_central_atoms[i].dist(template_central_atoms[i]), target_central_atoms[i].res_type, template_central_atoms[i].res_type)
                if target_central_atoms[i].dist(template_central_atoms[i]) < max_dist_for_pairing:
                    new_target_atoms.append(target_central_atoms[i])
                    new_template_atoms.append(template_central_atoms[i])
            if len(new_target_atoms) > 5:
                tm, rm = Superimpose.Superimpose._fit_pair_atoms([new_target_atoms,new_template_atoms])
                Superimpose.Superimpose._transform(tm, rm, self.template.atoms)
                rms = Superimpose.Superimpose._fit_atoms_rms([new_target_atoms, new_template_atoms])
                print 'new rms = %s over %s atom pairs < %s distance'%(rms, len(new_target_atoms), max_dist_for_pairing)
            else:
                print 'fewer than 5 atoms fit good in fit'
        
        for res in self.target.residues:
            res.central_atom.vtk_arg_list['atoms']['color'] = [0.0,0.0,1.0]
        for res in self.template.residues:
            res.central_atom.vtk_arg_list['atoms']['color'] = [0.0,0.0,1.0]
            
        target_central_atoms   = copy.deepcopy(list1)
        template_central_atoms = copy.deepcopy(list2)
        max_dist_for_pairing = 2.0
        for i in range(len(list1)):
            for atom in [list1[i], list2[i]]:
                atom.parent.vtk_arg_list['trace']['color'] = [1.0,1.0,1.0]
                atom.vtk_arg_list['atoms']['color'] = [1.0,1.0,1.0]
            if list1[i].dist(list2[i]) < max_dist_for_pairing:
                new_target_atoms.append(target_central_atoms[i])
                new_template_atoms.append(template_central_atoms[i])
                for atom in [list1[i], list2[i]]:
                    atom.parent.vtk_arg_list['trace']['color'] = [1.0,1.0,0.0]
                    atom.vtk_arg_list['atoms']['color'] = [1.0,1.0,0.0]
        if len(new_target_atoms) > 5:
            # now transform the original list of equivalent (non-gap) alpha carbon
            # atoms to get the overall rms
            tm, rm = Superimpose.Superimpose._fit_pair_atoms([new_target_atoms,new_template_atoms])
            Superimpose.Superimpose._transform(tm, rm, self.template.atoms)
            rms = Superimpose.Superimpose._get_atoms_rms([new_target_atoms, new_template_atoms])
            print 'final rms = %s over %s atom pairs < %s distance'%(rms, len(new_target_atoms), max_dist_for_pairing)
        else:
            print 'fewer than 5 atoms fit good in final fit'
        
        return rms, len(new_target_atoms)
                
