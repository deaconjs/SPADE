import copy
import math
import sys
import string
import time
import cPickle

from operator import xor

sys.path.append('../..')
import MolecularSystem
import Numeric

def multi_superimpose(systems):
    # first calculate inter-alpha-carbon distances

    ### PARMS ###
    core_only        = 2           # to limit to core alpha carbons use 2, else 1
    core_cutoff      = 8
    neighbor_thresh  = 8
    
    sort_by_sequence = 1

    shell_collection_start = 7
    shell_collection_end   = 14

    all_ten                = 0
    dsf                    = 0.35                 # distance scaling factor

    do_replicate           = 1
    replicate_thresh       = 0.025*len(systems)
    
    do_rms_filter          = 1
    max_alpha_rms          = 0.1  + 0.15*len(systems)
    max_beta_rms           = 0.5 + 0.25*len(systems)
    max_alpha_distance     = 1.0  + 0.5 *len(systems)
    max_beta_distance      = 2.5  + 1.0 *len(systems)
    
    max_lowres_hash_size = 5*len(systems)-5            # scale this per number of proteins
    
    if len(systems) < 5:
        min_cluster_dist = 1.0 + (0.5*len(systems))
    else:
        min_cluster_dist = 3.5
    
    min_hit_count    = 10

    max_cluster_dist = 0.35
    prefilter_translation_thresh = 5.0

    system_list = systems

    protein_list = []
    for system in systems:
        protein_list.append(system.ProteinList[0])

    for protein in protein_list:
        protein.assign_core_atom_distances(core_cutoff, neighbor_thresh, shell_collection_start, shell_collection_end)
    ref_chain = protein_list[0]
    src_chains = protein_list[1:]

    # create reference and source hash tables with all of the distances
    hash_tables = __create_4_musta_hashes(protein_list, core_cutoff, neighbor_thresh, dsf, do_replicate, replicate_thresh, core_only, sort_by_sequence)
    print 'l_ref-%s, h_ref-%s, l_src-%s, h_src-%s'%(len(hash_tables[0].keys()), len(hash_tables[1].keys()), len(hash_tables[2].keys()), len(hash_tables[3].keys()))
    
    # locate the bins with all structures, enumerate all possible combinations
    print 'creating, filtering complete buckets'
    complete_buckets = __create_complete_buckets(ref_chain, src_chains, hash_tables, max_alpha_distance, max_beta_distance, max_lowres_hash_size, max_alpha_rms, max_beta_rms, do_rms_filter)
    
    record = cluster_by_singles(ref_chain, src_chains, complete_buckets, min_cluster_dist, min_hit_count, max_cluster_dist, prefilter_translation_thresh)

    print '\n\n%s possible cores located\n\n'%(len(record))

    for item in record:
        print 'next item structure %s'%(item['structures'])
        for index in range(len(item['structures'])):
            for anum in item['ref_buckets'][index]:
                print ref_chain.atom_dict[anum].res_number,',',
        print ""
        for index in range(len(item['structures'])):
            for anum in item['src_buckets'][index]:
                print src_chains[item['structures'][index]-1].atom_dict[anum].res_number,',',
        print ""

    complete_solutions = []
    bin_count = 0
    indices = []
    bucket_counter = 0
    # each record represents one key which has had extended hits from multiple structures
    for j in range(len(src_chains)):
        indices.append(-1)
    for item in record:
        atomlist = []
        for anum in item['src_buckets'][0]:
            atom = protein_list[item['structures'][0]].atom_dict[anum]
            atomlist.append(atom)
        complete_solutions.append({'structure':item['structures'][0], 'ref':item['ref_buckets'][0], 'src':item['src_buckets'][0]})

    print '%s complete solutions located'%(len(complete_solutions))

    system_list[0].save_pdb('./ref.pdb')
    solution_count = 0
    for solution in complete_solutions:
        ref_atoms = []
        for anum in solution['ref']:
            atom = protein_list[0].atom_dict[anum]
            ref_atoms.append(atom)
        src_atoms = []
        for anum in solution['src']:
            atom = src_chains[solution['structure']-1].atom_dict[anum]
            src_atoms.append(atom)
        __atoms_to_pdb(ref_atoms, 'atoms1.pdb', "A")
        __atoms_to_pdb(src_atoms, 'atoms2.pdb', "A")

        list1, list2 = [], []
        for atom in ref_atoms:
            list1.append([atom.x, atom.y, atom.z])
        for atom in src_atoms:
            list2.append([atom.x, atom.y, atom.z])
        ta,rm =_fit_pair([list1,list2])
        atoms = copy.deepcopy(protein_list[solution['structure']].atoms)
        _transform(ta,rm, protein_list[solution['structure']].atoms)
        system_list[0].save_pdb('./sln%s_str0.pdb'%(solution_count))
        system_list[solution['structure']].save_pdb('./sln%s_str%s.pdb'%(solution_count, solution['structure']))
        protein_list[solution['structure']].atoms = atoms

        solution_count += 1
        if solution_count == 10:
            break


### TRANSFORMATIONS ###

def _transform(trns, rm, atom_list):
    """ apply this to the template. transforms to center, rotates, then transforms to the target """
    for atom in atom_list:
        x,y,z = atom.x+trns[1], atom.y+trns[3], atom.z+trns[5]
        rx,ry,rz = _rotate_point(x,y,z,rm)
        atom.x, atom.y, atom.z = rx-trns[0], ry-trns[2], rz-trns[4]

def _transform_coordinates(trns, rm, atom_list):
    for atom in atom_list:
        x,y,z = atom[0]+trns[1], atom[1]+trns[3], atom[2]+trns[5]
        rx,ry,rz = _rotate_point(x,y,z,rm)
        atom[0], atom[1], atom[2] = rx-trns[0], ry-trns[2], rz-trns[4]



### PDB SAVES ###

def __atoms_to_pdb(atomlist, filename, chain_name="A"):
    file = open(filename, 'w')
    for atom in atomlist:
        line = 'ATOM %6d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n'%(atom.atom_number,"CA","GLY",chain_name,atom.res_number,atom.x,atom.y,atom.z, " ")
        file.write(line)
    file.close()

def __bucket_to_pdb(bucket, filename, chain_name="A"):
    file = open(filename, 'w')
    for atom_data in bucket:
        x,y,z = atom_data['coordinates'][0], atom_data['coordinates'][1], atom_data['coordinates'][2]
        line = 'ATOM %6d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n'%(atom_data['atom_number'],"CA","GLY",chain_name,atom_data['atom_number'],x,y,z, " ")
        file.write(line)
    file.close()
        
def __coor_list_to_pdb(list, filename, chain_name="A"):
    file = open(filename, 'w')
    cntr = 1
    for sublist in list:
        x,y,z = sublist[0], sublist[1], sublist[2]
        line = 'ATOM %6d%5s %3s %1s%4d %11.3f %7.3f %7.3f%26s\n'%(cntr,"CA","GLY",chain_name,cntr,x,y,z, " ")
        file.write(line)
        cntr += 1
    file.close()
   

### CENTERING ###

def _center(frag, center_atoms=0):
    """ normalize centroids of multiple lists of atoms to the origin, (0,0,0) """
    cx,cy,cz = 0.0,0.0,0.0
    # collect the sum
    for atom in frag:
        cx, cy, cz = cx+atom.x, cy+atom.y, cz+atom.z
    points = len(frag)
    x,y,z = cx/points,cy/points,cz/points
    if center_atoms:
        # decrement by the average
        for atom in frag:
            atom.x, atom.y, atom.z = atom.x-x, atom.y-y, atom.z-z
    return -x, -y, -z
    
def _center_coordinates(frag, center_atoms=0):
    """ normalize centroids of multiple lists of [x,y,z] to the origin, (0,0,0) """
    cx,cy,cz = 0.0,0.0,0.0
    # collect the sum
    for i in range(len(frag)):
        cx, cy, cz = cx+frag[i][0], cy+frag[i][1], cz+frag[i][2]
    points = len(frag)
    x,y,z = cx/points,cy/points,cz/points
    if center_atoms:
        # decrement by the average
        for i in range(len(frag)):
            frag[i][0], frag[i][1], frag[i][2] = frag[i][0]-x, frag[i][1]-y, frag[i][2]-z
    return -x, -y, -z
    

### ROTATIONS ###


def _rotate(fragment, matin):
    """ rotates a single list of atoms by the given rotation matrix
    """
    for atom in fragment:
        x,y,z = atom.x,atom.y,atom.z
        atom.x,atom.y,atom.z = _rotate_point(x,y,z,matin)

def _rotate_coordinates(fragment, matin):
    """ rotates a single list of [x,y,z] lists by the given rotation matrix
    """
    for atom in fragment:
        x,y,z = atom[0],atom[1],atom[2]
        atom[0], atom[1], atom[2] = _rotate_point(x,y,z,matin)

def _rotate_point(x,y,z,matin):
    rx = x * matin[0][0] + y * matin[1][0] + z * matin[2][0]
    ry = x * matin[0][1] + y * matin[1][1] + z * matin[2][1]
    rz = x * matin[0][2] + y * matin[1][2] + z * matin[2][2]
    return rx,ry,rz


### SUPERPOSITIONING ###

def _fit(trns, rm, atom1, atom2):
    """ apply xyz translations and the rotation to atom2. Returns the distance between 
        the two atoms after transformation
    """
    coor_list = [[atom2.x,atom2.y,atom2.z]]
    _transform_coordinates(trns, rm, coor_list)
    return math.sqrt(((atom1.x-coor_list[0][0])**2) + ((atom1.y-coor_list[0][1])**2) + ((atom1.z-coor_list[0][2])**2))
    
def _fit_pair(fraglist, wtl=None):
    """ fit the second atom list to the first. returns a 2 component 
        vector [ref_x_transl, src_x_transl, ref_y_transl, src_y_transl, ref_z_transl, src_z_transl] and rotation_matrix
    """
    new_frags = copy.deepcopy(fraglist)
    ref_x_transl,ref_y_transl,ref_z_transl = _center_coordinates(new_frags[0], 1)
    src_x_transl,src_y_transl,src_z_transl = _center_coordinates(new_frags[1], 1)
    rm = _fit_pair_rotation(new_frags, 0, wtl, 0)
    return [ref_x_transl, src_x_transl, ref_y_transl, src_y_transl, ref_z_transl, src_z_transl], rm

def _fit_pair_atoms(atomlists):
    new_atomlists = copy.deepcopy(atomlists)
    ref_x_transl, ref_y_transl, ref_z_transl = _center(new_atomlists[0], 1)
    src_x_transl, src_y_transl, src_z_transl = _center(new_atomlists[1], 1)
    rm = _fit_pair_rotation(new_atomlists, 0, None, 1)
    return [ref_x_transl, src_x_transl, ref_y_transl, src_y_transl, ref_z_transl, src_z_transl], rm
    

def _fit_multiple(ref, srcs, wtl=None):
    """ fit each of the coordinate lists to the reference. returns a list of 2 component _fit_pair() vectors """
    vectors = []
    for i in range(len(srcs)):
        new_frags = [copy.deepcopy(ref), copy.deepcopy(srcs[i])]
        ref_x_transl,ref_y_transl,ref_z_transl = _center_coordinates(new_frags[0], 1)
        src_x_transl,src_y_transl,src_z_transl = _center_coordinates(new_frags[1], 1)
        rm = _fit_pair_rotation(new_frags, 0, wtl, 0)
        vectors.append([[ref_x_transl, src_x_transl, ref_y_transl, src_y_transl, ref_z_transl, src_z_transl], rm])
    return vectors

def superimpose(systemlist):
    ps = [systemlist[0].ProteinList[0].atoms, systemlist[1].ProteinList[0].atoms]
    _center(systemlist[0].ProteinList[0].atoms)
    _center(systemlist[1].ProteinList[0].atoms)
    rm = _fit_pair_rotation(ps, 0)
    print '%s atoms\n'%(systemlist[1].ProteinList[0].atoms)
    _rotate(systemlist[1].ProteinList[0].atoms, rm)
    rms = _fit_atoms_rms(ps)
    print 'rms = %s'%(rms)

def _fit_atoms_rms(atomlists):
    """ a list of k atoms from each of n structures """
    # create a table for multiple positioning
    table1 = []
    for atomlist in atomlists:
        table1.append([])
    for i in range(len(atomlists)):
        for atom in atomlists[i]:
            table1[i].append([atom.x,atom.y,atom.z])
    vectors = _fit_multiple(table1[0], table1[1:])                                  # calculate the transformation vectors

    for i in range(len(vectors)):
        _transform_coordinates(vectors[i][0], vectors[i][1], table1[i+1])           # transform the atoms by the obtained translation and rotation vectors

    # calculate the average distance
    # invert the table to put equivalent positions into lists within the big list
    table2 = []
    for i in range(len(atomlists[0])):
        table2.append([])
    for j in range(len(atomlists[0])):
        for i in range(len(atomlists)):
            table2[j].append(table1[i][j])
            i += 1
    grand_count   = 0
    sum = 0.0
    cnt = 0
    distances = []
    for column in table2:                                                           # for each set of equivalent positions
        for j in range(0,len(column)-1):
            for k in range(j+1,len(column)):                                        # calculate the distance between the two atoms
                distance = math.sqrt(((column[k][0]-column[j][0])**2) + ((column[k][1]-column[j][1])**2) + ((column[k][2]-column[j][2])**2))
                distances.append(distance)
                sum += distance
                cnt += 1

    average_dist = sum/(cnt+0.0)

    sum = 0.0
    for distance in distances:
        sum += (distance-average_dist)**2

    return math.sqrt(sum/(cnt-1.0))


def _get_atoms_rms(atomlists):
    """ a list of k atoms from each of n structures """
    # create a table for multiple positioning
    table1 = []
    for atomlist in atomlists:
        table1.append([])
    for i in range(len(atomlists)):
        for atom in atomlists[i]:
            table1[i].append([atom.x,atom.y,atom.z])
    # calculate the average distance
    # invert the table to put equivalent positions into lists within the big list
    table2 = []
    for i in range(len(atomlists[0])):
        table2.append([])
    for j in range(len(atomlists[0])):
        for i in range(len(atomlists)):
            table2[j].append(table1[i][j])
            i += 1
    grand_count   = 0
    sum = 0.0
    cnt = 0
    distances = []
    for column in table2:                                                           # for each set of equivalent positions
        for j in range(0,len(column)-1):
            for k in range(j+1,len(column)):                                        # calculate the distance between the two atoms
                distance = math.sqrt(((column[k][0]-column[j][0])**2) + ((column[k][1]-column[j][1])**2) + ((column[k][2]-column[j][2])**2))
                distances.append(distance)
                sum += distance
                cnt += 1

    average_dist = sum/(cnt+0.0)

    sum = 0.0
    for distance in distances:
        sum += (distance-average_dist)**2

    return math.sqrt(sum/(cnt-1.0))
    
def _fit_pair_rotation(fraglist, column=1, wt1=None, use_atoms=1):
    """ _fit_rotation finds the best fitting rotation
    """
    if len(fraglist[0]) != len(fraglist[1]):
        print 'frags different lengths %s, %s'%(len(fraglist[0]), len(fraglist[1]))
        return -1
    rm   = [[0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0]]
    umat = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    n = len(fraglist[0])
    if n<2:
        return 0
    if use_atoms == 1:
        for i in range(3):
            for j in range(n):
                dict1 = {0:fraglist[0][j].x, 1:fraglist[0][j].y, 2:fraglist[0][j].z}
                if wt1:
                    w = wt1[j] * dict1[i]
                    umat[i][0] += w * fraglist[1][j].x
                    umat[i][1] += w * fraglist[1][j].y
                    umat[i][2] += w * fraglist[1][j].z
                else:
                    umat[i][0] += dict1[i] * fraglist[1][j].x
                    umat[i][1] += dict1[i] * fraglist[1][j].y
                    umat[i][2] += dict1[i] * fraglist[1][j].z
    else:
        for i in range(3):
            for j in range(n):
                f = fraglist[0][j][i]
                if wt1:
                    w = wt1[j] * f
                    umat[i][0] += w * fraglist[1][j][0]
                    umat[i][1] += w * fraglist[1][j][1]
                    umat[i][2] += w * fraglist[1][j][2]
                else:
                    umat[i][0] += f * fraglist[1][j][0]
                    umat[i][1] += f * fraglist[1][j][1]
                    umat[i][2] += f * fraglist[1][j][2]

    # now fit
    rotcnt    = 0
    rot       = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    turmat    = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    c         = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
    coup      = [0.0,0.0,0.0]
    dir       = [0.0,0.0,0.0]
    step      = [0.0,0.0,0.0]
    v         = [0.0,0.0,0.0]
    jmax      = 30
    rtsum     = umat[0][0] + umat[1][1] + umat[2][2]
    delta     = 0.0
    minimum_d = 0.00000000000000000001
    gRotMat   = []
    for i in range(100):
        gRotMat.append([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    
    for ncyc in range(jmax):
        nsteep = 3
        nrem = ncyc-nsteep*round(ncyc/nsteep)
        if nrem == 0:
            step = [0.0,0.0,0.0]
            clep = 1.0
        # couple
        coup[0] = umat[1][2]-umat[2][1]
        coup[1] = umat[2][0]-umat[0][2]
        coup[2] = umat[0][1]-umat[1][0]
        cle     = math.sqrt(coup[0]**2 + coup[1]**2 + coup[2]**2)
        # gradient vector is now -coup
        gfac    = (cle/clep)**2
        rtsump  = rtsum
        deltap  = delta
        clep    = cle
        if cle < minimum_d:
            break
        # Step vector conjugate to previous
        stp = 0.0
        for i in range(3):
            step[i] = coup[i]+step[i]*gfac
            stp    += step[i]**2
        stp = 1.0/math.sqrt(stp)
        # Normalized step
        for i in range(3):
            dir[i] = stp * step[i]
        # couple resolved along step direction
        stcoup = coup[0] * dir[0] + coup[1]*dir[1] + coup[2]*dir[2]
        # Component of UMAT along direction
        ud = 0.0
        for i in range(3):
            for j in range(3):
                ud += umat[i][j]*dir[i]*dir[j]
        tr = umat[0][0] + umat[1][1] + umat[2][2] - ud
        ta = math.sqrt(tr**2 + stcoup**2)
        cs = tr/ta
        sn = stcoup/ta
        # if cs< 0 then position is unstable, so dont stop
        if cs < 0.0 and abs(sn)<minimum_d:
            break
        # Turn matrix for correcting rotation
        #    symmetric part first
        ac = 1.0-cs
        for i in range(3):
            v[i] = ac*dir[i]
            for j in range(3):
                turmat[i][j] = v[i] * dir[j]
            turmat[i][i] += cs
            v[i] = dir[i] * sn
        #    asymmetric part now
        turmat[0][1] -= v[2]
        turmat[1][2] -= v[0]
        turmat[2][0] -= v[1]
        turmat[1][0] += v[2]
        turmat[2][1] += v[0]
        turmat[0][2] += v[1]
        # update total rotation matrix
        for i in range(3):
            for j in range(3):
                c[i][j] = 0.0
                for k in range(3):
                    c[i][j] += turmat[i][k]*rot[k][j]
        rot = copy.deepcopy(c)
        for i in range(3):
            for j in range(3):
                c[i][j] = 0.0
                for k in range(3):
                    c[i][j] += turmat[i][k]*umat[k][j]
        umat = copy.deepcopy(c)
        rtsum = umat[0][0] + umat[1][1] + umat[2][2]
        delta = rtsum - rtsump
        if abs(delta) < minimum_d:
            break
        # store intermediate copies of rotation matrix
        if rotcnt >= 100:
            print 'out of bounds'
            sys.exit()
        gRotMat[rotcnt] = copy.deepcopy(rot)
        rotcnt += 1
    rsum = rtsum
    # copy rotation matrix for output
    if (column):
        for i in range(3):
            for j in range(3):
                rm[j][i] = rot[i][j]
    else:
        rm[0][0],rm[0][1],rm[0][2] = rot[0][0], rot[0][1], rot[0][2]
        rm[1][0],rm[1][1],rm[1][2] = rot[1][0], rot[1][1], rot[1][2]
        rm[2][0],rm[2][1],rm[2][2] = rot[2][0], rot[2][1], rot[2][2]
    return rm



def __create_4_musta_hashes(protein_list, core_cutoff, neighbor_thresh, dsf, do_replicate, replicate_thresh, core_only, sort_by_sequence):
    print 'collecting hash keys'
    l_src_hash = {}
    h_src_hash = {}
    l_ref_hash = {}
    h_ref_hash = {}
    b_src_hash = {}
    b_ref_hash = {}
    src_ind = -1
    lores_counter = 0
    hires_counter = 0

    # make sure pseudo-sidechains are there
    for protein in protein_list:
        try:
            protein.residues[0].pseudo_sidechain
        except:
            protein.fill_pseudo_sidechains(self, type=0)

    # create the hash table storing sets of coordinates for multidimensional transformation generation
    # the keys are strings of interatomic distances and symmetries
    for src in protein_list:
        start_time = time.clock()
        src_ind += 1
        if src_ind == 0:
            l_hash = l_ref_hash
            h_hash = h_ref_hash
        else:
            l_hash = l_src_hash
            h_hash = h_src_hash

        if core_only == 2:
            atom_list = src.get_core_alpha_carbons(core_cutoff, neighbor_thresh)
            sec_list  = src.get_central_atom_list()
            print '%s of %s are core atoms'%(len(atom_list), len(sec_list))
            for al in atom_list:
                print '%s,'%(al.res_number),
            print ""
        elif core_only == 1:
            atom_list = src.get_central_atom_list()
        # calculate hash keys and store
        # note that the coordinates are all normalized such that the query atom is on the origin
        for central_atom in atom_list:
            #print 'next Ca'
            x1,y1,z1 = central_atom.x, central_atom.y, central_atom.z
            pb1      = central_atom.parent.pseudo_sidechain
            ncps = central_atom.data['nearby_cps']
            ncp_count = len(ncps)
            good_its = 0
            bad_its  = 0
            if ncp_count >= 4:
                for a2 in range(0,ncp_count-3):
                    x2,y2,z2 = ncps[a2].x, ncps[a2].y, ncps[a2].z
                    pb2 = ncps[a2].parent.pseudo_sidechain
                    r_a2_dist = math.sqrt(((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))
                    r_b2_dist = math.sqrt(((pb2.x-pb1.x)**2) + ((pb2.y-pb1.y)**2) + ((pb2.z-pb1.z)**2))
                    for a3 in range(a2+1, ncp_count-2):
                        #if ncps[a2] not in ncps[a3].data['nearby_cps']:
                        #    continue
                        x3,y3,z3 = ncps[a3].x, ncps[a3].y, ncps[a3].z
                        pb3 = ncps[a3].parent.pseudo_sidechain
                        r_a3_dist = math.sqrt(((x3-x1)**2) + ((y3-y1)**2) + ((z3-z1)**2))
                        r_b3_dist = math.sqrt(((pb3.x-pb1.x)**2) + ((pb3.y-pb1.y)**2) + ((pb3.z-pb1.z)**2))
                        for a4 in range(a3+1, ncp_count-1):
                            #if ncps[a2] not in ncps[a4].data['nearby_cps']:
                            #    continue
                            #if ncps[a3] not in ncps[a4].data['nearby_cps']:
                            #    continue
                            x4,y4,z4 = ncps[a4].x, ncps[a4].y, ncps[a4].z
                            pb4 = ncps[a4].parent.pseudo_sidechain
                            r_a4_dist = math.sqrt(((x4-x1)**2) + ((y4-y1)**2) + ((z4-z1)**2))
                            r_b4_dist = math.sqrt(((pb4.x-pb1.x)**2) + ((pb4.y-pb1.y)**2) + ((pb4.z-pb1.z)**2))
                            # start by sorting the four added atoms by their distance to the query
                            dist_list = [r_a2_dist, r_a3_dist, r_a4_dist]
                            beta_dist = [r_b2_dist, r_b3_dist, r_b4_dist]
                            atom_list = [ncps[a2], ncps[a3], ncps[a4]]
                            beta_list = [pb2, pb3, pb4]
                            atom_num_list = [atom_list[0].atom_number, atom_list[1].atom_number, atom_list[2].atom_number]
                            # sort just this one
                            if sort_by_sequence:
                                sorted_list = [atom_list[0].atom_number, atom_list[1].atom_number, atom_list[2].atom_number]
                            else:
                                sorted_list = [r_a2_dist, r_a3_dist, r_a4_dist]
                            sorted_list.sort()
                            # and catch the indices so as to use all of the above lists
                            f = [0,0,0]
                            for i in range(len(sorted_list)):
                                for j in range(len(dist_list)):
                                    if sort_by_sequence:
                                        if atom_num_list[j] == sorted_list[i]:
                                            f[i] = j
                                    else:
                                        if dist_list[j] == sorted_list[i]:
                                            f[i] = j

                            # collect the coordinates for faster access and readability
                            xs = [atom_list[f[0]].x, atom_list[f[1]].x, atom_list[f[2]].x]
                            ys = [atom_list[f[0]].y, atom_list[f[1]].y, atom_list[f[2]].y]
                            zs = [atom_list[f[0]].z, atom_list[f[1]].z, atom_list[f[2]].z]

                            xb = [beta_list[f[0]].x, beta_list[f[1]].x, beta_list[f[2]].x]
                            yb = [beta_list[f[0]].y, beta_list[f[1]].y, beta_list[f[2]].y]
                            zb = [beta_list[f[0]].z, beta_list[f[1]].z, beta_list[f[2]].z]

                            new_distance_list  = [dist_list[f[0]],
                                                  dist_list[f[1]],
                                                  dist_list[f[2]],
                                                  math.sqrt(((xs[0]-xs[1])**2) + ((ys[0]-ys[1])**2) + ((zs[0]-zs[1])**2)),
                                                  math.sqrt(((xs[0]-xs[2])**2) + ((ys[0]-ys[2])**2) + ((zs[0]-zs[2])**2)),
                                                  math.sqrt(((xs[1]-xs[2])**2) + ((ys[1]-ys[2])**2) + ((zs[1]-zs[2])**2))]

                            beta_distances     = [beta_dist[f[0]],
                                                  beta_dist[f[1]],
                                                  beta_dist[f[2]],
                                                  math.sqrt(((xb[0]-xb[1])**2) + ((yb[0]-yb[1])**2) + ((zb[0]-zb[1])**2)),
                                                  math.sqrt(((xb[0]-xb[2])**2) + ((yb[0]-yb[2])**2) + ((zb[0]-zb[2])**2)),
                                                  math.sqrt(((xb[1]-xb[2])**2) + ((yb[1]-yb[2])**2) + ((zb[1]-zb[2])**2))]

                            hires_dl_string  = '%2.1f_%2.1f_%2.1f_%2.1f_%2.1f_%2.1f_'%(new_distance_list[0],
                                                                                       new_distance_list[1],
                                                                                       new_distance_list[2],
                                                                                       new_distance_list[3],
                                                                                       new_distance_list[4],
                                                                                       new_distance_list[5])

                            # store the bin distances to replicate
                            lowres_dl_bins   = [[],[],[],[],[],[]]
                            lowres_dlstrings = []
                            for i in range(len(new_distance_list)):
                                lowres_dl_bins[i].append(math.floor(dsf*new_distance_list[i]))
                                if do_replicate:
                                    if new_distance_list[i] % dsf <= replicate_thresh:              # if the distance is just over an integer change
                                        lowres_dl_bins[i].append((math.floor(dsf*new_distance_list[i]))-1)
                                    elif new_distance_list[i]%dsf >= (1.0-replicate_thresh):
                                        lowres_dl_bins[i].append((math.floor(dsf*new_distance_list[i]))+1)

                            if do_replicate:
                                for i0 in lowres_dl_bins[0]:
                                    for i1 in lowres_dl_bins[1]:
                                        for i2 in lowres_dl_bins[2]:
                                            for i3 in lowres_dl_bins[3]:
                                                for i4 in lowres_dl_bins[4]:
                                                    for i5 in lowres_dl_bins[5]:
                                                        lowres_dlstrings.append('%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_'%(i0,i1,i2,i3,i4,i5))
                            else:
                                lowres_dlstrings.append('%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_'%(lowres_dl_bins[0][0],
                                                                                                lowres_dl_bins[1][0],
                                                                                                lowres_dl_bins[2][0],
                                                                                                lowres_dl_bins[3][0],
                                                                                                lowres_dl_bins[4][0],
                                                                                                lowres_dl_bins[5][0]))
                            normal_vector = [ys[0]*zs[1] - zs[0]*ys[1],
                                             zs[0]*xs[1] - xs[0]*zs[1],
                                             xs[0]*ys[1] - ys[0]*xs[1]]

                            dot_4 = xs[2]*normal_vector[0] + ys[2]*normal_vector[1] + zs[2]*normal_vector[2]

                            # add a marker for symmetry (above or below the plane defined by the other three atoms)
                            if dot_4>=0:
                                app = '1_'
                            else:
                                app = '0_'

                            hires_dl_string += app

                            lores_counter += len(lowres_dlstrings)
                            hires_counter += 1

                            for dlstring in lowres_dlstrings:
                                dlstring += app
                                try:
                                    l_hash[dlstring].append(hires_dl_string)
                                except KeyError:
                                    l_hash[dlstring] = [hires_dl_string]

                            try:
                                h_hash[hires_dl_string].append({'structure':src_ind, 'beta_distances':beta_distances, 'atoms':[central_atom.atom_number, atom_list[f[0]].atom_number, atom_list[f[1]].atom_number, atom_list[f[2]].atom_number]})
                            except KeyError:
                                h_hash[hires_dl_string] =     [{'structure':src_ind, 'beta_distances':beta_distances, 'atoms':[central_atom.atom_number, atom_list[f[0]].atom_number, atom_list[f[1]].atom_number, atom_list[f[2]].atom_number]}]


        end_time = time.clock()
        # record how many items have been added to the coloring lists of 
        print 'structure done in %5.3f seconds'%(end_time - start_time)
    sum = 0
    cnt = 0
    for l_key in l_hash.keys():
        sum += len(l_hash[l_key])
        cnt += 1
    print 'average lowres storage size = %s'%(sum/(cnt+0.0))

    print 'done key generation. %s lowres, %s hires items stored'%(lores_counter, hires_counter)
    return [l_ref_hash, h_ref_hash, l_src_hash, h_src_hash]




def __create_3_musta_hashes(protein_list, core_cutoff, neighbor_thresh, dsf, do_replicate, replicate_thresh, core_only, sort_by_sequence):
    print 'collecting hash keys'
    l_src_hash = {}
    h_src_hash = {}
    l_ref_hash = {}
    h_ref_hash = {}
    b_src_hash = {}
    b_ref_hash = {}
    src_ind = -1
    lores_counter = 0
    hires_counter = 0

    # make sure pseudo-sidechains are there
    for protein in protein_list:
        try:
            protein.residues[0].pseudo_sidechain
        except:
            protein.fill_pseudo_sidechains(0)

    # create the hash table storing sets of coordinates for multidimensional transformation generation
    # the keys are strings of interatomic distances and symmetries
    for src in protein_list:
        start_time = time.clock()
        src_ind += 1
        if src_ind == 0:
            l_hash = l_ref_hash
            h_hash = h_ref_hash
        else:
            l_hash = l_src_hash
            h_hash = h_src_hash

        if core_only == 2:
            atom_list = src.get_core_alpha_carbons(core_cutoff, neighbor_thresh)
            sec_list  = src.get_central_atom_list()
            print '%s of %s are core atoms'%(len(atom_list), len(sec_list))
            for al in atom_list:
                print '%s,'%(al.res_number),
            print ""
        elif core_only == 1:
            atom_list = src.get_central_atom_list()
        # calculate hash keys and store
        # note that the coordinates are all normalized such that the query atom is on the origin
        for central_atom in atom_list:
            #print 'next Ca'
            x1,y1,z1 = central_atom.x, central_atom.y, central_atom.z
            pb1      = central_atom.parent.pseudo_sidechain
            ncps = central_atom.data['nearby_cps']
            ncp_count = len(ncps)
            good_its = 0
            bad_its  = 0
            if ncp_count >= 3:
                for a2 in range(0,ncp_count-2):
                    x2,y2,z2 = ncps[a2].x, ncps[a2].y, ncps[a2].z
                    pb2 = ncps[a2].parent.pseudo_sidechain
                    r_a2_dist = math.sqrt(((x2-x1)**2) + ((y2-y1)**2) + ((z2-z1)**2))
                    r_b2_dist = math.sqrt(((pb2.x-pb1.x)**2) + ((pb2.y-pb1.y)**2) + ((pb2.z-pb1.z)**2))
                    for a3 in range(a2+1, ncp_count-1):
                        if ncps[a2] not in ncps[a3].data['nearby_cps']:
                            continue
                        x3,y3,z3 = ncps[a3].x, ncps[a3].y, ncps[a3].z
                        pb3 = ncps[a3].parent.pseudo_sidechain
                        r_a3_dist = math.sqrt(((x3-x1)**2) + ((y3-y1)**2) + ((z3-z1)**2))
                        r_b3_dist = math.sqrt(((pb3.x-pb1.x)**2) + ((pb3.y-pb1.y)**2) + ((pb3.z-pb1.z)**2))

                        # start by sorting the four added atoms by their distance to the query
                        dist_list = [r_a2_dist, r_a3_dist]
                        beta_dist = [r_b2_dist, r_b3_dist]
                        atom_list = [ncps[a2], ncps[a3]]
                        beta_list = [pb2, pb3]
                        atom_num_list = [atom_list[0].atom_number, atom_list[1].atom_number]
                        # sort just this one
                        if sort_by_sequence:
                            sorted_list = [atom_list[0].atom_number, atom_list[1].atom_number]
                        else:
                            sorted_list = [r_a2_dist, r_a3_dist]
                        sorted_list.sort()
                        # and catch the indices so as to use all of the above lists
                        f = [0,0]
                        for i in range(len(sorted_list)):
                            for j in range(len(dist_list)):
                                if sort_by_sequence:
                                    if atom_num_list[j] == sorted_list[i]:
                                        f[i] = j
                                else:
                                    if dist_list[j] == sorted_list[i]:
                                        f[i] = j

                        # collect the coordinates for faster access and readability
                        xs = [atom_list[f[0]].x, atom_list[f[1]].x]
                        ys = [atom_list[f[0]].y, atom_list[f[1]].y]
                        zs = [atom_list[f[0]].z, atom_list[f[1]].z]

                        xb = [beta_list[f[0]].x, beta_list[f[1]].x]
                        yb = [beta_list[f[0]].y, beta_list[f[1]].y]
                        zb = [beta_list[f[0]].z, beta_list[f[1]].z]

                        new_distance_list  = [dist_list[f[0]],
                                              dist_list[f[1]],
                                              math.sqrt(((xs[0]-xs[1])**2) + ((ys[0]-ys[1])**2) + ((zs[0]-zs[1])**2))]

                        beta_distances     = [beta_dist[f[0]],
                                              beta_dist[f[1]],
                                              math.sqrt(((xb[0]-xb[1])**2) + ((yb[0]-yb[1])**2) + ((zb[0]-zb[1])**2))]

                        hires_dl_string  = '%2.1f_%2.1f_%2.1f_'%(new_distance_list[0],
                                                                 new_distance_list[1],
                                                                 new_distance_list[2])

                        # store the bin distances to replicate
                        lowres_dl_bins   = [[],[],[]]
                        lowres_dlstrings = []
                        for i in range(len(new_distance_list)):
                            lowres_dl_bins[i].append(math.floor(dsf*new_distance_list[i]))
                            if do_replicate:
                                if new_distance_list[i] % dsf <= replicate_thresh:              # if the distance is just over an integer change
                                    lowres_dl_bins[i].append((math.floor(dsf*new_distance_list[i]))-1)
                                elif new_distance_list[i]%dsf >= (1.0-replicate_thresh):
                                    lowres_dl_bins[i].append((math.floor(dsf*new_distance_list[i]))+1)

                        if do_replicate:
                            for i0 in lowres_dl_bins[0]:
                                for i1 in lowres_dl_bins[1]:
                                    for i2 in lowres_dl_bins[2]:
                                        lowres_dlstrings.append('%2.0f_%2.0f_%2.0f_'%(i0,i1,i2))
                        else:
                            lowres_dlstrings.append('%2.0f_%2.0f_%2.0f_'%(lowres_dl_bins[0][0],
                                                                          lowres_dl_bins[1][0],
                                                                          lowres_dl_bins[2][0]))

                        lores_counter += len(lowres_dlstrings)
                        hires_counter += 1

                        for dlstring in lowres_dlstrings:
                            try:
                                l_hash[dlstring].append(hires_dl_string)
                            except KeyError:
                                l_hash[dlstring] = [hires_dl_string]

                        try:
                            h_hash[hires_dl_string].append({'structure':src_ind, 'beta_distances':beta_distances, 'atoms':[central_atom.atom_number, atom_list[f[0]].atom_number, atom_list[f[1]].atom_number]})
                        except KeyError:
                            h_hash[hires_dl_string] =     [{'structure':src_ind, 'beta_distances':beta_distances, 'atoms':[central_atom.atom_number, atom_list[f[0]].atom_number, atom_list[f[1]].atom_number]}]


        end_time = time.clock()
        print 'structure done in %5.3f seconds'%(end_time - start_time)
    sum = 0
    cnt = 0
    for l_key in l_hash.keys():
        sum += len(l_hash[l_key])
        cnt += 1
    print 'average lowres storage size = %s'%(sum/(cnt+0.0))

    print 'done key generation. %s lowres, %s hires items stored'%(lores_counter, hires_counter)
    return [l_ref_hash, h_ref_hash, l_src_hash, h_src_hash]




def __create_complete_buckets(ref_chain, src_chains, hash_tables, max_alpha_distance, max_beta_distance, max_lowres_hash_size, max_alpha_rms=0.5, max_beta_rms=0.5, do_rms_filter=1):
    # query with the reference keys -- then put bins with representatives from all structures into final_hash
    # when structures have multiple representatives, enumerate all possible combinations
    # the ref_hash data is in the following format, where a 5 atom hash hit is shown
    #    try:      # check for KeyError
    #        ref_hash[dlstring].append({'structure':src_ind,   'atoms':[central_atom, atom_list[f[0]], atom_list[f[1]], atom_list[f[2]]], atom_list[f[3]]]})
    #    except:
    #        ref_hash[dlstring] = [{'structure':src_ind,       'atoms':[central_atom, atom_list[f[0]], atom_list[f[1]], atom_list[f[2]]], atom_list[f[3]]]}]
    l_ref_hash = hash_tables[0]
    h_ref_hash = hash_tables[1]
    l_musta_hash = hash_tables[2]
    h_musta_hash = hash_tables[3]
    bin_count = 0
    indices = []
    for j in range(len(src_chains)):
        indices.append(0)

    # filter out missing ref keys
    hit = 0
    miss = 0
    for key in l_ref_hash.keys():
        try:
            l_musta_hash[key]
        except KeyError:
            del l_ref_hash[key]
            miss += 1
        else:
            hit += 1
    print 'of the remaining %s ref keys, %s match the %s src keys'%(hit + miss, hit, len(l_musta_hash.keys()))

    # filter out source keys with too few hashed items
    hit = 0
    miss = 0
    for key in l_ref_hash.keys():
        for i in range(len(indices)):
            indices[i] = 0
        for h_ref_key in l_musta_hash[key]:
            for set in h_musta_hash[h_ref_key]:
                indices[set['structure']-1] = 1
        if 0 in indices:
            hit += 1
            del l_ref_hash[key]
        else:
            miss += 1
    print 'of the remaining %s ref keys, %s match src keys with all structures present'%(hit + miss, miss)

    # filter out source keys with too many hash entries
    hit = 0
    miss = 0
    for key in l_ref_hash.keys():
        if len(l_musta_hash[key]) >= max_lowres_hash_size or len(l_ref_hash[key]) >= max_lowres_hash_size:
            hit += 1
            del l_ref_hash[key]
        else:
            miss += 1
    print 'of the remaining %s ref keys, %s match src keys without > %s hash entries'%(hit + miss, miss, max_lowres_hash_size)

    complete_buckets = []
    for j in range(len(src_chains)):
        indices[j] = -1
    full_hits = 0
    too_many  = 0
    alpha_keep  = 0
    beta_keep   = 0
    alpha_throw = 0
    beta_throw  = 0
    gross_throw = 0

    iteration_counter = 0
    for l_ref_key in l_ref_hash.keys():
        for h_ref_key in l_ref_hash[l_ref_key]:                                 # for each high-res key in each lo-res reference bin
            h_ref_tokens = string.split(h_ref_key[:-1], '_')
            for ref_hit in h_ref_hash[h_ref_key]:                               # for each set of atoms in the reference hash
                decent_hashes = []                                              # collect sets of atoms from source structures that are close to this one
                for h_mst_key in l_musta_hash[l_ref_key]:                       # for each hi-res source key by the indicated lo-res signature
                    h_mst_tokens = string.split(h_mst_key[:-1], '_')            # make sure the hi-res keys are compatible
                    diff_distance = 0.0
                    for i in range(len(h_mst_tokens)):                          # count the distance; these are single-decimal precision distances
                        diff_distance += abs(string.atof(h_mst_tokens[i]) - string.atof(h_ref_tokens[i]))**2
                    diff_distance = math.sqrt(diff_distance)
                    if diff_distance < max_alpha_distance:
                        for set in h_musta_hash[h_mst_key]:                     # for each set under this hi-res key
                            iteration_counter += 1
    print '%s reference to source comparisons to consider'%(iteration_counter)
    start_time = time.clock()

    cntr = 0

    # iterate first through the lores ref keys. 
    for l_ref_key in l_ref_hash.keys():
        for h_ref_key in l_ref_hash[l_ref_key]:                                 # for each high-res key in each lo-res reference bin
            h_ref_tokens = string.split(h_ref_key[:-1], '_')
            for ref_hit in h_ref_hash[h_ref_key]:                               # for each set of atoms in the reference hash
                decent_hashes = []                                              # collect sets of atoms from source structures that are close to this one
                for h_mst_key in l_musta_hash[l_ref_key]:                       # for each hi-res source key by the indicated lo-res signature
                    h_mst_tokens = string.split(h_mst_key[:-1], '_')            # make sure the hi-res keys are compatible
                    diff_distance = 0.0
                    for i in range(len(h_mst_tokens)):                          # count the distance; these are single-decimal precision distances
                        diff_distance += (string.atof(h_mst_tokens[i]) - string.atof(h_ref_tokens[i]))**2
                    diff_distance = math.sqrt(diff_distance)
                    #print 'a dis %4.2f'%(diff_distance)
                    if diff_distance < max_alpha_distance:
                        for set in h_musta_hash[h_mst_key]:                     # if the alphas are compatible, check the betas
                            cntr += 1
                            if cntr % 10000 == 0:
                                print '10k'
                            diff_distance = 0.0
                            for i in range(len(ref_hit['beta_distances'])):     # count the distance
                                diff_distance += (ref_hit['beta_distances'][i]-set['beta_distances'][i])**2
                            diff_distance = math.sqrt(diff_distance)
                            #print '         b      dis %4.2f'%(diff_distance)
                            if diff_distance <= max_beta_distance:
                                decent_hashes.append(set)                       # if the betas are good too, then store the set
                if len(decent_hashes) < len(src_chains):                            # pass if not enough structures
                    continue
                present_call, present_buckets = [],[]                               # create some bins to watch present calls, and to collect hash bins that have each structure
                for i in range(len(src_chains)):
                    present_call.append(0)
                    present_buckets.append([])
                for i in range(len(decent_hashes)):                                 # for each of the source motifs falling in this reference's key
                    present_call[decent_hashes[i]['structure']-1] = 1               # mark the present call for that structure
                    present_buckets[decent_hashes[i]['structure']-1].append(i)      # and collect the index of this element in the hashes' bin
                if 0 not in present_call:
                    for ref_hit in h_ref_hash[h_ref_key]:                           # go through any duplicates for this hi-res index
                        full_hits += 1
                        bins_to_create = 1
                        for i in range(len(present_buckets)):                       # first count the number of buckets to create
                            bins_to_create *= len(present_buckets[i])
                        bin_count += bins_to_create                                 # record the count here
                        for i in range(bins_to_create):                             # enumerate all of the combinations iteratively
                            mod = i
                            for j in range(len(src_chains)):
                                bins_below = 1
                                for k in range(j+1,len(src_chains)):
                                    bins_below *= len(present_buckets[k])
                                indices[j] = int(math.floor(mod/bins_below))
                                mod = math.floor(mod%bins_below)
                            atomlist1 = [[]]
                            atomlist2 = [[]]
                            for atomnum in ref_hit['atoms']:
                                atomlist1[0].append(ref_chain.atom_dict[atomnum])
                                atomlist2[0].append(ref_chain.atom_dict[atomnum].parent.pseudo_sidechain)
                            for j in range(len(src_chains)):
                                atomlist1.append([])
                                atomlist2.append([])
                                for atomnum in decent_hashes[present_buckets[j][indices[j]]]['atoms']:
                                    atomlist1[-1].append(src_chains[j].atom_dict[atomnum])
                                    atomlist2[-1].append(src_chains[j].atom_dict[atomnum].parent.pseudo_sidechain)
                            if do_rms_filter:
                                rms = _fit_atoms_rms(atomlist1)
                                #print 'a rms %4.2f'%(rms)
                                if rms <= max_alpha_rms:
                                    rms = _fit_atoms_rms(atomlist2)
                                    #print '    b rms %4.2f'%(rms)
                                    if rms <= max_beta_rms:
                                        complete_buckets.append({'data':[copy.deepcopy(ref_hit)], 'key':l_ref_key})
                                        for j in range(len(src_chains)):
                                            complete_buckets[-1]['data'].append(copy.deepcopy(decent_hashes[present_buckets[j][indices[j]]]))
                            else:
                                complete_buckets.append({'data':[copy.deepcopy(ref_hit)], 'key':l_ref_key})
                                for j in range(len(src_chains)):
                                    complete_buckets[-1]['data'].append(copy.deepcopy(decent_hashes[present_buckets[j][indices[j]]]))

    end_time = time.clock()
    # record how many items have been added to the coloring lists of 
    print 'done complete bins %5.3f seconds'%(end_time - start_time)
    #print 'alpha - kept %s threw %s %s'%(alpha_keep, alpha_throw, alpha_keep/(alpha_throw+0.0))
    #print 'beta  - kept %s threw %s %s'%(beta_keep, beta_throw, beta_keep/(beta_throw+0.0))
    print '%s ref hits to %s complete bins (%s under %s %s RMS)'%(full_hits, bin_count, len(complete_buckets), max_alpha_rms, max_beta_rms)
    return complete_buckets


def cluster_by_singles(ref_chain, src_chains, complete_buckets, min_cluster_dist, min_hit_count, max_cluster_dist, prefilter_translation_thresh):
    # use the complete buckets as seeds to cluster buckets. If a cluster exceeds min_hit_count, 
    # store the structure index (i), and the reference and source buckets in record. Index record by the
    # key of the seed bucket
    record = []
    # do each chain separately
    for bucket in complete_buckets:
        bucket['available'] = 1

    print '%s buckets'%(len(complete_buckets))

    bucket_backup = copy.deepcopy(complete_buckets)

    for i in range(1,len(src_chains)+1):
        print 'structure %s %s'%(i,src_chains[i-1].parent.filename)
        complete_buckets = copy.deepcopy(bucket_backup)
        # collect the total match list for this src-ref combo. these are simply all nonredundant pairs of atoms that appear
        total_match_list = []
        total_match_index = 0
        for j in range(len(complete_buckets)):
            complete_buckets[j]['match_indices'] = []                         # record the matches that are consistent with the transform
            for k in range(len(complete_buckets[j]['data'][i]['atoms'])):
                current_match_index = 0
                for match in total_match_list:
                    # skip if the pair is already present
                    if match['src'].atom_number == complete_buckets[j]['data'][i]['atoms'][k] and match['ref'].atom_number == complete_buckets[j]['data'][0]['atoms'][k]:
                        complete_buckets[j]['match_indices'].append(current_match_index)
                        break
                    current_match_index += 1
                else:
                    # collect ref-src atom pairs
                    complete_buckets[j]['match_indices'].append(total_match_index)
                    total_match_list.append({'src':src_chains[i-1].atom_dict[complete_buckets[j]['data'][i]['atoms'][k]],'ref':ref_chain.atom_dict[complete_buckets[j]['data'][0]['atoms'][k]]})
                    total_match_index += 1
        print '%s in total_match_list'%(len(total_match_list))

        added_a_new_one  = 1
        while added_a_new_one == 1:                                               # iterate until no more clustering is available
            print 'iterating'
            start_time = time.clock()

            for bucket in complete_buckets:                                       # reset cluster centers to available for clustering and pickup
                if bucket['available'] == 2:
                    bucket['available'] = 1
            added_a_new_one = 0

            for bucket in complete_buckets:
                if len(bucket['match_indices']) < 3:
                    bucket['available'] = 0

            # take the consensus of the used indices from the total_match_list
            consensus_transform_indices = []
            for bucket in complete_buckets:
                if bucket['available'] == 1:
                    for match_index in bucket['match_indices']:
                        if match_index not in consensus_transform_indices:
                            consensus_transform_indices.append(match_index)
            print '%s transform indices left - calculating masks'%(len(consensus_transform_indices))

            # first update the match_indices; the transformation mask
            added_matches = 0
            failed_matches = 0
            for bucket in complete_buckets:                                       # for each bucket
                if bucket['available'] == 1:
                    list1,list2 = [],[]
                    for transform_index in bucket['match_indices']:
                        atom = total_match_list[transform_index]['ref']
                        list1.append([atom.x,atom.y,atom.z])
                        atom = total_match_list[transform_index]['src']
                        list2.append([atom.x,atom.y,atom.z])
                
                    ta, rm = _fit_pair([list1, list2])                                # calculate the transformation vector
                    new_ta = [ta[0],ta[1],ta[2],ta[3],ta[4],ta[5]]                    # make a copy of it
                    new_rm = [[rm[0][0],rm[0][1],rm[0][2]],[rm[1][0],rm[1][1],rm[1][2]],[rm[2][0],rm[2][1],rm[2][2]]]
                    bucket['ta'] = [ta[0],ta[1],ta[2],ta[3],ta[4],ta[5]]
                    bucket['rm'] = [[rm[0][0],rm[0][1],rm[0][2]],[rm[1][0],rm[1][1],rm[1][2]],[rm[2][0],rm[2][1],rm[2][2]]]
                    transformation_vector = {'translation':new_ta, 'rotation':new_rm}
                    fitting_indices = []
                    for match_index in consensus_transform_indices:                  # test each match for consistency with this bucket
                        if match_index not in bucket['match_indices']:
                            dist = _fit(new_ta, new_rm, total_match_list[match_index]['ref'], total_match_list[match_index]['src'])
                            if dist < min_cluster_dist:
                                added_matches += 1
                                bucket['match_indices'].append(match_index)               # store the match_indices list for this bucket
                            else:
                                failed_matches += 1
                                
            left_avail_count = 0
            for bucket in complete_buckets:
                if bucket['available'] == 1:
                    left_avail_count += 1
            print 'done transformation masks - merged %s of %s clusters - begin clustering'%(added_matches, added_matches+failed_matches)
            
            still_available = 0                                                   # count how many are left
            for bucket in complete_buckets:
                if bucket['available'] == 1:
                    still_available += 1

            match_count      = 0
            last_match_count = 1
            while still_available > 0:                                            # cluster all of the groups by iteratively passing through
                first_bucket = None
                bi1 = -1
                for bucket in complete_buckets:                                   # get the first available bucket for cluster initialization
                    bi1 += 1
                    if bucket['available'] == 1:
                        first_bucket = bucket
                        bucket['available'] = 2                                   # mark it as unavailable for extending itself
                        still_available -= 1
                        break
                else:
                    break
                if still_available == 0:
                    break
                bi2 = -1
                sec_bucket = None

                for bucket in complete_buckets:                                   # compare to all of the remaining available buckets
                    # pre-filter by translational distance of the source molecule
                    sec_bucket = bucket
                    distance_sum = 0.0
                    for k in [1,3,5]:
                        distance_sum += (first_bucket['ta'][k] - sec_bucket['ta'][k])**2
                    distance_sum = math.sqrt(distance_sum)
                    if distance_sum > prefilter_translation_thresh:
                        continue
                    bi2 += 1
                    if sec_bucket['available'] == 1:                              # if available and not self,
                        found_it = 0
                        for ind1 in first_bucket['match_indices']:                # count the number shared
                            for ind2 in sec_bucket['match_indices']:
                                if ind1 == ind2:
                                    found_it += 1
                        if len(first_bucket['match_indices']) < len(sec_bucket['match_indices']):
                            inconsistent_pairs = len(first_bucket['match_indices']) - found_it
                        else:
                            inconsistent_pairs = len(sec_bucket['match_indices']) - found_it
                        if len(first_bucket) > len(sec_bucket):
                            lesserval = len(sec_bucket)
                        else:
                            lesserval = len(first_bucket)
                        if inconsistent_pairs <= max_cluster_dist * (lesserval + 0.0):
                            added_a_new_one = 1
                            sec_bucket['available'] = 0
                            still_available -= 1
                            for new_index in sec_bucket['match_indices']:
                                if new_index not in first_bucket['match_indices']:
                                    dist = _fit(first_bucket['ta'], first_bucket['rm'], total_match_list[new_index]['ref'], total_match_list[new_index]['src'])
                                    if dist > min_cluster_dist:
                                        break
                            else:
                                first_bucket['match_indices'].append(new_index)               # store the match_indices list for this bucket
                            if still_available == 0:
                                break

                # reduce any ambiguous match_indices by taking only the best of the possiblities
                for bucket in complete_buckets:
                    beforelen = len(bucket['match_indices'])
                    fitting_indices = first_bucket['match_indices']
                    first_bucket['match_indices'] = []
                    good_indices = []
                    while len(fitting_indices) > 0:
                        continue_test = 0
                        first_index = fitting_indices[0]
                        for f_ind in range(1,len(fitting_indices)):
                            if xor(total_match_list[first_index]['ref'] == total_match_list[fitting_indices[f_ind]]['ref'], total_match_list[first_index]['src'] == total_match_list[fitting_indices[f_ind]]['src']):
                                d1 = _fit(first_bucket['ta'], first_bucket['rm'], total_match_list[first_index]['ref'], total_match_list[first_index]['src'])
                                d2 = _fit(first_bucket['ta'], first_bucket['rm'], total_match_list[fitting_indices[f_ind]]['ref'], total_match_list[fitting_indices[f_ind]]['src'])
                                if d1 > d2:
                                    fitting_indices = fitting_indices[1:]
                                    break
                                else:
                                    fitting_indices = fitting_indices[:f_ind] + fitting_indices[f_ind+1:]
                                    break
                        else:
                            good_indices.append(fitting_indices[0])
                            fitting_indices = fitting_indices[1:]
                    for good_index in good_indices:
                        first_bucket['match_indices'].append(good_index)               # store in a the match_indices list for this bucket
                    afterlen = len(bucket['match_indices'])

                beforelen = len(first_bucket['match_indices'])
 
                for match_index in first_bucket['match_indices']:                 # append located atoms
                    for anum in first_bucket['data'][0]['atoms']:
                        atom = ref_chain.atom_dict[anum]
                        if total_match_list[match_index]['ref'].atom_number == atom.atom_number:
                            break
                    else:

                        for anum in first_bucket['data'][0]['atoms']:
                            atom = ref_chain.atom_dict[anum]
                            if total_match_list[match_index]['ref'].atom_number == atom.atom_number:
                                break
                        else:
                            first_bucket['data'][0]['atoms'].append(total_match_list[match_index]['ref'].atom_number)
                            first_bucket['data'][i]['atoms'].append(total_match_list[match_index]['src'].atom_number)

                afterlen = len(first_bucket['match_indices'])
                #if beforelen != afterlen:
                #    print 'added atoms %s %s'
                #if len(first_bucket['data'][0]['atoms']) != len(first_bucket['data'][i]['atoms']):
                #    print 'unequal numbers of atoms between ref and src... memory leak ! ! ! ! ! ! ! ! ! !'
            match_count    = 0
            for bucket in complete_buckets:
                if len(bucket['data'][0]['atoms']) > match_count:
                    match_count = len(bucket['data'][0]['atoms'])
            end_time = time.clock()
            print 'iteration done in %5.3f seconds'%(end_time - start_time)
            print 'largest match %s residues'%(match_count)

        for bucket in complete_buckets:
            if len(bucket['data'][0]['atoms']) >= min_hit_count:
                record.append({'structures':[i],'ref_buckets':[copy.deepcopy(bucket['data'][0]['atoms'])],'src_buckets':[copy.deepcopy(bucket['data'][i]['atoms'])]})
    print '%s records found'%(len(record))
    return record


if __name__ == '__main__':
    start_time = time.clock()
    systems = []
    #system_list = ['1L1OB.pdb','1L1OA.pdb']
    system_list = ['1L1OB.pdb','1KAW.pdb','1L1OA.pdb','1L1OC.pdb']
    for system in system_list:
        systems.append(MolecularSystem.System(None))
        systems[-1].load_pdb(system)
    multi = 1
    if multi:
        for system in systems:
            system.ProteinList[0].fill_pseudo_sidechains(1)
            system.ProteinList[0].fill_neighbors_lists(0.35,15.0)
        multi_superimpose(systems)
    else:
        superimpose(systems)
        systems[0].save_pdb('atoms3.pdb')
        systems[1].save_pdb('atoms4.pdb')

    end_time = time.clock()
    print 'program done in %5.3f seconds'%(end_time - start_time)



