import math
import random

def auto_decompose(pchain, viewer=None):
    """
        this function identifies closed networks of residues. 
        Identification occurs by first detecting potential core residues.
        From those find tightly packed networks of residues. The constraints
        are very high so as to not initiate crystals at protein-protein 
        interfaces, which are more loosely packed than residues participating
        in hydrophobic cores. Initial crystals have a minimum number of residues
        and must be within a close distance to be considered. Extention constraints
        are less filtering, but should not allow walk across domains or around barriers.
        The closed domains are analyzed for entanglement, and unentangled
        networks are differentially displayed.
    """
    pchain.fill_pseudo_sidechains(0)
    pchain.fill_neighbors_lists(0.35,15.0)
    pchain.fill_secondary()
    if viewer != None:
        viewer.color_trace_by_residue_feature(pchain, 'shielding')
    secondary_structure_filter = 0    # allows only residues in secondary structure units to be used in the core
    merging = 1
  
    # once the initial core crystal is located, do not allow extention outside of the range of this core.
    disallow_ext_beyond_crystal_range = 1

    # used to select initial core fragments
    min_core_sa                                       = 0.7
    Cavg_factor                                       = 0.7   # neighbors necessary for inclusion in core selection
    Havg_factor                                       = 0.5   # avg regional hydrophobicity for inclusion in core selection
    minimum_crystal_start                             = 5
    network_start_distance                            = 12.0
    minimum_crystal_ext                               = 3
    network_ext_distance                              = 11.0
    blocked_range_percentage                          = 0.8   # denotes the total range of the region to consider
    min_alpha_for_merging                             = 1.0
    selection_choice                                  = 8
    min_core_merging_interactions                     = 1000 # the minimum number of interactions that must be observed to merge cores in the first merging
    core_ext_distance                                 = 8.0  # the minimum distance for scanning the core merging interactions
    min_entanglement_for_merging_captured_domains     = 0.8  # the percent of residues in one domain that must be within the range of another domain to consider it for merging cores
    max_entanglement_for_not_merging_captured_domains = 0.2  # the max percent of residues from the capturing domain to consider it separate
    min_entanglement_for_merging_entangled_domains    = 0.35 # each domain of a pair must be entangled to this degree to be merged
    min_contacts_for_merging_entangled_domains        = 29   # each domain of a pair must be entangled to this degree to be merged
    minimum_core_rez_for_ranges                       = 10   # the number of contiguous residues needed for including a segment of the protein in one cores range.
    alpha_merge_distance                              = 5.0  # the threshold distance for merging domains based on many close alpha carbons
    alpha_merge_min_cnt                               = 0.12 # the number of alpha carbons that must be within the threshold distance
    # sheets are currently disabled
    adj_beta_merg_cnt                                 = 1.85 # the number of adjacent, alternating residues located in strands needed to merge two cores
    interface_deletion_threshold                      = 0.35  # the threshold percentage of core residues that must be remaining after extending all other cores, in order to delete interfaces that are identified as cores

    multi_domain     = 0
    contained_domain = 0

    #avail_core stuff takes care of the core that is made available in the beginning of execution
    avail_core_count = 0
    temp_core_count  = 0
    perm_core_count  = 0
    avail_core_list = []
    temp_core_list = []
    perm_core_list = []
    mstln            = len(pchain.residues)
    for j in range(mstln):
        avail_core_list.append(0)             # avail_core_list keeps track of all the residues that have been used
        temp_core_list.append(0)             # temp_core_list keeps track of the residues that have been used by the current core
        perm_core_list.append(0)             # perm_core_list keeps permentant track of all core residues so that temp_core_list can be reset

    # all_core_list stuff takes care of 
    all_core_list = []
    for j in range(mstln):
        all_core_list.append(0)
    all_core_count = 0

    # full_core_list stuff keeps all of the chains that have been found, separated by -1's
    full_core_list = []
    for k in range(2*mstln):
        full_core_list.append(0)
    full_core_count = 0
    
    full_core_list = []
    for k in range(2*mstln):
        full_core_list.append(0)
    
    number_of_cores = 0

    # current_core stuff takes care of the crystallizing core
    current_core_count = 0
    current_core_list = []
    for k in range(mstln):
        current_core_list.append(0)

    crystal_watcher = 0
    last_avail_core_count = 0
    # first identify the average number of contacts, and the average hydrophobicity
    hscore = 0
    cscore = 0
    hb_dict = {'A': 8,'R': 0,'N': 0,'D': 0,'C': 4,'Q': 0,'E': 0,'G': 3,'H': 1,'I': 12,'L': 11,'K': 2,'M': 9,'F': 10,'P': 0,'S': 0,'T': 5,'W': 7,'Y': 6,'V': 13,'-': 0}
    for res in pchain.residues:
        Hscore = Cscore = 0.0
        
        for n in res.neighbors:
            Cscore += len(pchain.residues[n].neighbors)
            Cscore += len(res.neighbors)
            Hscore += hb_dict[pchain.residues[n].res_type1]
            Hscore += hb_dict[res.res_type1]
        cscore += Cscore / (len(res.neighbors)*2)
        hscore += Hscore / (len(res.neighbors)*2)

    Havg = hscore / len(pchain.residues)
    Cavg = cscore / len(pchain.residues)

    # Use these standards to identify potential core residues. Make a list of residues above average for both
    j = 0
    for res in pchain.residues:
        Cscore = 0.0
        Hscore = 0.0
        for p in res.neighbors:
            Cscore += len(pchain.residues[p].neighbors)
            Cscore += len(res.neighbors)
            Hscore += hb_dict[pchain.residues[p].res_type1]
            Hscore += hb_dict[res.res_type1]
        cscore = Cscore / (len(res.neighbors)*2)
        hscore = Hscore / (len(res.neighbors)*2)
        print '%s > %s and %s > %s and %s >= %s'%(len(res.neighbors), (Cavg_factor*Cavg), res.features['shielding'], min_core_sa, hscore, Havg_factor*Havg)
        if len(res.neighbors) > (Cavg_factor * Cavg) and res.features['shielding'] > min_core_sa and hscore >= Havg_factor * Havg:
            if secondary_structure_filter and res.features['secondary'] == 'C':
                continue
            avail_core_list[j] = 1
            avail_core_count  += 1
            perm_core_list[j]  = 1
            perm_core_count   += 1
        else:
            for atom in res.atoms:
                atom.selected = 0
        j += 1
    if viewer != None:
        viewer.color_atoms_by_residue_feature(pchain, 'shielding')
        viewer.atoms_display_on()        
    
    win_counter = 0
    # Pt pt1, pt2
    temp_start_size = 0.0

    core_count = avail_core_count

    crystal_range_start = 0
    crystal_range_end = 0

    # now locate domains
    while 1:
        start_size = network_start_distance  # 12
        for n in range(avail_core_count):
            pass_counter = 0
            current_core_count = 0
            start_holder = 0
            for k in range(mstln):
                if avail_core_list[k]:
                    pass_counter += 1
                    if pass_counter == n+1:
                        start_holder = k
                        break
            print 'start holder %s'%(start_holder)
            for j in range(start_holder, mstln):
                j_res = pchain.residues[j]
                if avail_core_list[j]:
                    if current_core_count == 0:                                      # if its the first to be added to the current crystal
                        current_core_list[current_core_count] = j                    # initiate the core list
                        current_core_count += 1
                        crystal_range_start = j
                    else:                                                            # else add it to the core list if its networked with it
                        win_counter = 0
                        for k in range(current_core_count):
                            for partner in pchain.residues[current_core_list[k]].neighbors:
                                if pchain.residues[partner] == j_res:
                                    p1x = pchain.residues[partner].pseudo_sidechain.x
                                    p1y = pchain.residues[partner].pseudo_sidechain.y
                                    p1z = pchain.residues[partner].pseudo_sidechain.z
                                    p2x = j_res.x
                                    p2y = j_res.y
                                    p2z = j_res.z
                                    if start_size > math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2):
                                        crystal_range_end = j                        # collect the points
                                        win_counter += 1                             # then tally the observation
                        if win_counter == current_core_count:                      # if its in each member of the current core
                            current_core_list[current_core_count] = j                # then add it to the core
                            current_core_count += 1
            if current_core_count >= minimum_crystal_start:
                break
        print 'found %s in current core, %s required'%(current_core_count, minimum_crystal_start)
        if current_core_count < minimum_crystal_start:                              # if a minimum crystal was not found in the search, 
            break
        else:
            temp_core_count = 0
            for h in range(mstln):                                                   # first reset the temp_core_list to reflect the perm_core_list
                temp_core_list[h] = perm_core_list[h]
                if perm_core_list[h]:
                    temp_core_count += 1
            for h in range(current_core_count):
                avail_core_list[current_core_list[h]] = 0
                avail_core_count -= 1
                temp_core_list[current_core_list[h]] = 0
                temp_core_count -= 1
            # color the current crystal
            if viewer != None:
                for rind in range(len(pchain.residues)):
                    if rind in current_core_list:
                        for atom in pchain.residues[rind].atoms:
                            atom.vtk_arg_list['atoms']['color'] = [1.0,1.0,0.0]
                viewer.update_view()
            # finally, extend the crystal
            while 1:                                                                 # loop until no new residues are added
                last_avail_core_count = temp_core_count                              # use to monitor whether a residue was added - break later if not
                for h in range(crystal_range_start, crystal_range_end):              # dont extend beyond the current crystal range
                    hres = pchain.residues[h]
                    if temp_core_list[h]:
                        win_counter = 0
                        for k in range(current_core_count):
                            kres = pchain.residues[current_core_list[k]]
                            for m in range(len(kres.neighbors)):
                                if kres.neighbors[m] == hres:
                                    p1x = kres.pseudo_sidechain.x    # also limit by the distance between beta carbons
                                    p1y = kres.pseudo_sidechain.y
                                    p1z = kres.pseudo_sidechain.z
                                    p2x = hres.pseudo_sidechain.x
                                    p2y = hres.pseudo_sidechain.y
                                    p2z = hres.pseudo_sidechain.z
                                    if math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2) < network_ext_distance:
                                        win_counter += 1                                # then tally the observation
                        if win_counter >= minimum_crystal_ext:                     # if its in so many members of the current core
                            avail_core_list[h] = 0
                            avail_core_count  -= 1
                            temp_core_list[h]  = 0
                            temp_core_count   -= 1
                            current_core_list[current_core_count] = h                # then add it to the core
                            current_core_count+= 1
                if temp_core_count == last_avail_core_count:                       # break when none are added in a full loop
                    break
            crystal_watcher = 0

            # ok, done with the core - now add it to the full_list
            for m in range(current_core_count):
                full_core_list[full_core_count+number_of_cores+m] = current_core_list[m]
            full_core_count += current_core_count
            number_of_cores += 1
            full_core_list[full_core_count+number_of_cores-1] = -1
            for m in range(current_core_count):
                avail_core_list[current_core_list[m]] = 0
            for m in range(all_core_count, all_core_count+current_core_count):
                all_core_list[m] = current_core_list[m - all_core_count]
            all_core_count += current_core_count
            for m in range(current_core_count):
                current_core_list[m] =0
            current_core_count = 0
        if avail_core_count < minimum_crystal_start:                               # if there are not enough residues left in the avail_core_list, break
            break
    return    
    # first rearrange the data so that its easier to handle.
    # store each core in its own index of a two-dimensional array
    cores = []
    for j in range(number_of_cores):
        cores_count.append(0)
        cores.append([])
        for k in range(mstln):
          cores[-1].append([])
    fcr_cntr = 0
    for j in range(number_of_cores):
        ncr_cntr = 0
        while 1:
            cores[j][ncr_cntr] = full_core_list[fcr_cntr]
            fcr_cntr += 1
            ncr_cntr += 1
            if full_core_list[fcr_cntr] == -1:
                fcr_cntr += 1
                cores_count[j] = ncr_cntr
                break

    # now try to delete any interfaces
    # for each core, try to delete it by extending all of the rest of the interfaces one residue in each direction.
    # if the targeted core gets exempted, then it is an interface.

    # delete the most removed first,then re-evaluate the rest

    remaining_core_percent = []
    for i in range(number_of_cores):
        remaining_core_percent.append(0.0)
    remaining_core_watcher = 0

    while (1):
        lowest_remaining_core = 1.5
        for k in range(number_of_cores):
            temp_core_size = cores_count[k]
            for m in range(cores_count[k]):
                for n in range(number_of_cores):
                    if n == k:
                        continue
                    else:
                        for p in range(cores_count[n]):
                            if cores[k][m] == cores[n][p]+1 or cores[k][m] == cores[n][p] or cores[k][m] == cores[n][p]-1:
                                temp_core_size -= 1
                                break
                            elif cores[k][m] < 15 or cores[k][m] > len(pchain.residues)-15:    # if its in the first or last 15 amino acids, delete it
                                temp_core_size -= 1
                                break
            remaining_core_percent[k] = temp_core_size/float(cores_count[k])
            if remaining_core_percent[k] < lowest_remaining_core:
                remaining_core_watcher = k
                lowest_remaining_core = remaining_core_percent[k]
        if lowest_remaining_core > interface_deletion_threshold:             # now, if the least remaining core is above the threshold percentage, break, else delete that core
            break
        else:                                                                # delete the cor
            for n in range(remaining_core_watcher+1, number_of_cores):
                for p in range(mstln):
                    cores[n-1][p] = cores[n][p]
                cores_count[n-1] = cores_count[n]                            # move cores_count info
            number_of_cores -= 1

    ranges = []
    for i in range(number_of_cores):
        ranges.append(0)
    for j in range(number_of_cores):
        ranges[-1].append([0,0])
    for j in range(number_of_cores):
        ranges[j][0] = 0
        ranges[j][1] = 0
   
    # now calculate their ranges
    for j in range(number_of_cores):
        start_searching = math.ceil(cores_count[j] * ((1.0 - blocked_range_percentage) / 2.0))
        blocked = 0                                                             # first find the lower end
        for m in range(len(pchain.residues)):
            for n in range(cores_count[j]):                                 # for all of the residues in this core
                if cores[j][n] == m:                                      # if this one is the mth one
                    if blocked >= start_searching:                        # if were already past the blocked residues
                        ranges[j][0] = m                                    # this one starts at m
                        break                                               # finished with the mth residue
                    blocked += 1
            if ranges[j][0] != 0:                                           # if the lower range has been found, break
                break

        blocked = 0
        # now find the upper end
        for m in range(mstln-1,-1,-1):
            for n in range(cores_count[j]):
                if cores[j][n] == m:
                    if blocked >= start_searching:
                        ranges[j][1] = m
                        break
                    blocked += 1
            if ranges[j][1] != 0:
                break


    # now that we have the cores, refine them.

    # now, for any cores  that entangle with two other cores that are otherwise unentangled, 
    #      merge the entangled one with the one it is most entangled with.
    #      but only merge those residues that are in the range of the latter.
    #      this is a search to remove domain interfaces that might be tightly packed
    # write out the number of core residues from each region are also in the region of another core.
    entanglement_matrix = []
    for j in range(number_of_cores):
        entanglement_matrix.append([])
        for k in range(number_of_cores):
            entanglement_matrix[-1].append(0)
    for j in range(number_of_cores):
        for k in range(number_of_cores):
            if k == j:
                continue
            else:
                for m in range(ranges[j][0],ranges[j][1]+1):
                    if m > ranges[k][0] and m < ranges[k][1]:                  # if that residue is in the range of the kth core
                        for p in range(cores_count[j]):
                            if cores[j][p] == m:
                                entanglement_matrix[j][k] += 1                           # then add it to the entanglement matrix


    # write out the number of residues from each region are in contact with another core.
    complex_entanglement_matrix = []
    for j in range(number_of_cores):
        complex_entanglement_matrix.append([])
        for k in range(number_of_cores):
            complex_entanglement_matrix[-1].append(0)

    for j in range(number_of_cores):
        for k in range(number_of_cores):
            if k == j:
                continue
            else:
                for m in range(ranges[j][0], ranges[j][1]+1):
                    for n in range(cores_count[k]):
                        for p in range(len(pchain.residues[cores[k][n]].neighbors)):
                            if pchain.residues[cores[k][n]].neighbors[p] == m:
                                complex_entanglement_matrix[j][k] += 1

    # now, delete any domains that entangle any other two, otherwise unentangled domains
    tmp_numcr = 0
    deletion_made = 0
    while (1):
        tmp_numcr = number_of_cores
        deletion_made = 0
        for j in range(number_of_cores):
            for k in range(j+1,number_of_cores):                   # for each pair of cores
                for m in range(number_of_cores):                   # for each third core
                    if m == j or m == k:
                        continue
                    else:
                        # if the j and k domains are not entangled
                        if  ( ((float(entanglement_matrix[j][k]) / float(cores_count[j])) < min_entanglement_for_merging_entangled_domains) 
                           and ((float(entanglement_matrix[k][j]) / float(cores_count[k])) < min_entanglement_for_merging_entangled_domains) ):
                              # if the j and m domains are entangled though
                            if   ( ((((float(entanglement_matrix[j][m]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains) )
                                 and (((float(entanglement_matrix[m][j]) / float(cores_count[m])) >= min_entanglement_for_merging_entangled_domains) ) )
                              # and the k and m domains are also entangled
                               and  ((((float(entanglement_matrix[m][k]) / float(cores_count[m])) >= min_entanglement_for_merging_entangled_domains) )
                                 and (((float(entanglement_matrix[k][m]) / float(cores_count[k])) >= min_entanglement_for_merging_entangled_domains) ) ) ):
                                for n in range(m+1,number_of_cores):      # then delete the mth domain
                                    for p in range(mstln):                # move cores info
                                        cores[n-1][p] = cores[n][p]
                                    ranges[n-1][0] = ranges[n][0]         # move ranges info
                                    ranges[n-1][1] = ranges[n][1]
                                    cores_count[n-1] = cores_count[n]     # move cores_count info
                                number_of_cores -= 1
                                deletion_made = 1
                                break
                if deletion_made:
                    break
            if deletion_made:
                break

        # now recalculate the entanglement and complex entanglement matrices
        for j in range(number_of_cores):
            for k in range(number_of_cores):
                entanglement_matrix[j][k] = 0

        for j in range(number_of_cores):
            for k in range(number_of_cores):
                if k == j:
                    continue
                else:
                    for m in range(ranges[j][0],ranges[j][1]+1):
                        if m > ranges[k][0] and m < ranges[k][1]:
                            for p in range(cores_count[j]):
                                if cores[j][p] == m:
                                    entanglement_matrix[j][k] += 1
 
        # write out the number of residues from each region are in contact with another core.
        for j in range(number_of_cores):
            for k in range(number_of_cores):
                complex_entanglement_matrix[j][k] = 0
     
        for j in range(number_of_cores):
            for k in range(number_of_cores):
                if k == j:
                    continue
                else:
                    for m in range(ranges[j][0], ranges[j][1]+1):
                        mres = pchain.residues[m]
                        for n in range(cores_count[k]):
                            for partner in pchain.residues[cores[k][n]].neighbors:
                                if partner == mres:
                                    complex_entanglement_matrix[j][k] += 1
        if number_of_cores == tmp_numcr:
            break

    # first, merge any cores that are highly contacted with each other - extention by multiple sets of contacts
    #  less figuratively, count the number of contacts between the two cores. if above a threshold, merge
    # the threshold should be pretty high; we'll recalculate this later, once we've removed any interfaces that might be problematic
    # for each core
    interaction = 0
    reset_button = 0
    interactions_matrix = []
    for j in range(number_of_cores):
        interaction_matrix.append([])
        for k in range(number_of_cores):
            interaction_matrix.append(0)
    
    for j in range(number_of_cores):
        for k in range(number_of_cores):
            interaction = 0
            for m in range(cores_count[j]):                             # for each contact of the first
                for partner in pchain.residues[cores[j][m]].neighbors:
                    for p in range(cores_count[k]):
                        if partner == cores[k][p]:                      # if it contacts a residue of the second core
                            p1x = pchain.residues[cores[j][m]].pseudo_sidechain.x
                            p1y = pchain.residues[cores[j][m]].pseudo_sidechain.y
                            p1z = pchain.residues[cores[j][m]].pseudo_sidechain.z
                            p2x = pchain.residues[cores[k][p]].pseudo_sidechain.x
                            p2y = pchain.residues[cores[k][p]].pseudo_sidechain.y
                            p2z = pchain.residues[cores[k][p]].pseudo_sidechain.z
                            if math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2) < core_ext_distance:
                                interaction  += 1                           # then tally the observation
                                interactions_matrix[j][k] += 1
                                interactions_matrix[k][j] += 1

    min_alpha_distances = []
    for j in range(number_of_cores):
        min_alpha_distances.append([])
        for k in range(number_of_cores):
            min_alpha_distances.append(10000.0)

    # now calculate the distance of the closest alpha carbons between the two cores.
    for j in range(number_of_cores):
        for k in range(number_of_cores):
            if k == j:
                continue
            else:
                for m in range(cores_count):
                    for n in range(cores_count[x]):
                        p1x = pchain.residues[cores[j][m]].x
                        p1y = pchain.residues[cores[j][m]].y
                        p1z = pchain.residues[cores[j][m]].z
                        p2x = pchain.residues[cores[k][n]].x
                        p2y = pchain.residues[cores[k][n]].y
                        p2z = pchain.residues[cores[k][n]].z
                        d = math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2)
                        if d < min_alpha_distances[j][k]:
                            min_alpha_distances[j][k] = d
                            min_alpha_distances[k][j] = d

    # now calculate the number of alpha carbon atoms within a threshold distance
    close_alphas = []
    for j in range(number_of_cores):
        close_alphas.append([])
        for k in range(number_of_cores):
            close_alphas.append(0)
    
    for j in range(number_of_cores):
        for k in range(number_of_cores):
            if k == j:
                continue
            else:
                for m in range(cores_count[j]):
                    for n in range(cores_count[k]):
                        p1x = pchain.residues[cores[j][m]].x
                        p1y = pchain.residues[cores[j][m]].y
                        p1z = pchain.residues[cores[j][m]].z
                        p2x = pchain.residues[cores[k][n]].x
                        p2y = pchain.residues[cores[k][n]].y
                        p2z = pchain.residues[cores[k][n]].z
                        if math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2) < alpha_merge_distance:
                            close_alphas[j][k] += 1

    # now count the number of residues that are in beta strands from different cores, but belong to the same beta sheet
    same_sheet_rez = []
    for j in range(number_of_cores):
        same_sheet_rez.append([])
        for k in range(number_of_cores):
            same_sheet_rez.append(0)
    
    """
    for j in range(number_of_cores):
        for k in range(number_of_cores):
            if k == j:
                continue
            else:
                # for each core residue of the first domain, find out which sheet its in.
                # then step through the core residues of the second domain. if they are in the same sheet, increment same_sheet_rez
                owner_sheet = 0
                for m in range(cores_count[j]):
                    if self.residues[cores[j][m]].features['secondary'] == 'B':                  # if its in a strand
                        for n in range(pchain.sheets.p[i].FragmentSetCnt):
                            for p in range(pchain.sheets.p[i].s[n].FragCnt):
                                for q in range(pchain.sheets.p[i].s[n].f[p].PosCnt):
                                    if cores[j][m] == pchain.sheets.p[i].s[n].f[p].p[q]:
                                        owner_sheet = n
                        for n in range(cores_count[k]):
                            if pchain.residues[cores[k][n]].features['secondary'] == 'B':
                                for p in range(pchain.sheets.p[i].FragmentSetCnt):
                                    for q in range(pchain.sheets.p[i].s[p].FragCnt):
                                        for r in range(pchain.sheets.p[i].s[p].f[q].PosCnt):
                                            if cores[k][n] == pchain.sheets.p[i].s[p].f[q].p[r]:
                                                if owner_sheet == p:
                                                    same_sheet_rez[j][k] += 1
    """
    for q in range(mstln):
        full_core_list[q] = 0
    k_ctr = -1
    for q in range(number_of_cores):
        for r in range(cores_count[q]):
            k_ctr += 1
            full_core_list[k_ctr] = cores[q][r]
        k_ctr += 1
        full_core_list[k_ctr] = -1


    # now, using the contacts matrix and the entanglement matrix, 
    #   decipher whether one core contained within another core is a continuation, or is actually two domains.
    #   if its a continuation, merge the cores, else leave it alone
    #   also decipher whether or not two independent domains should be entangled
    selection = 0
    while 1:
        if not merging:
            break
        reset_button = 0
        for j in range(number_of_cores):
            for k in range(number_of_cores):
                if j == k:
                    continue
                else:
                    # this switch statement allows control over which domains are merged. It can be set up to merge different types of domains at different times.
                    case = selection_choice
                    if case == 0:                         # merge cases where one domain contains the other or where two domains are entangled
                        selection = ((((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= min_entanglement_for_merging_captured_domains)
                                 and ((float(entanglement_matrix[k][j]) / float(cores_count[k])) <= max_entanglement_for_not_merging_captured_domains)
                                 and (interactions_matrix[j][k] >= min_core_merging_interactions))
                                 or (((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains)
                                  and ((float(entanglement_matrix[k][j]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains)))
                    elif case in [1,3,6]:                         # merge cases where one domain contains another
                        selection = (((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= min_entanglement_for_merging_captured_domains)
                                and ((float(entanglement_matrix[k][j]) / float(cores_count[k])) <= max_entanglement_for_not_merging_captured_domains))
                                            #and (interactions_matrix[j][k] >= min_core_merging_interactions))
                    elif case == 2:                         # merge cases where two domains entangle each other
                        selection = (((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains)
                                 and ((float(entanglement_matrix[k][j]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains))
                    elif case in [4,5]:
                        selection = ((complex_entanglement_matrix[j][k] >= min_contacts_for_merging_entangled_domains)
                                  or (complex_entanglement_matrix[k][j] >= min_contacts_for_merging_entangled_domains))
                    elif case == 7:
                        selection = ((((complex_entanglement_matrix[j][k] >= min_contacts_for_merging_entangled_domains)
                                    or ((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains))
                                  and ((complex_entanglement_matrix[k][j] >= min_contacts_for_merging_entangled_domains)
                                    or ((float(entanglement_matrix[k][j]) / float(cores_count[k])) >= min_entanglement_for_merging_entangled_domains)))
                                and (min_alpha_distances[j][k] < min_alpha_for_merging))
                    elif case == 8:
                        # if both domains are entangled, by one or the other measure, and have the minimum alphas within a distance,
                        # or if both domains are three times the normal degree of entanglement
                        selection = ((((((complex_entanglement_matrix[j][k] >= min_contacts_for_merging_entangled_domains)
                                      or ((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= min_entanglement_for_merging_entangled_domains))
                                    and ((complex_entanglement_matrix[k][j] >= min_contacts_for_merging_entangled_domains)
                                      or ((float(entanglement_matrix[k][j]) / float(cores_count[k])) >= min_entanglement_for_merging_entangled_domains)))
                                    or ((interactions_matrix[j][k] >= min_core_merging_interactions)
                                     or (interactions_matrix[k][j] >= min_core_merging_interactions)))
                                 and (((float(close_alphas[j][k])/float(cores_count[j])) >= alpha_merge_min_cnt)
                                   or ((float(close_alphas[k][j])/float(cores_count[k])) >= alpha_merge_min_cnt)))
                                 or (((complex_entanglement_matrix[j][k] >= 3 * min_contacts_for_merging_entangled_domains)
                                   or ((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= 3 * min_entanglement_for_merging_entangled_domains))
                                 and ((complex_entanglement_matrix[k][j] >= 3 * min_contacts_for_merging_entangled_domains)
                                   or ((float(entanglement_matrix[k][j]) / float(cores_count[k])) >= 3 * min_entanglement_for_merging_entangled_domains))))
                    elif case == 9:
                        # if they are on opposite sides of a beta sheet
                        selection = ((float(same_sheet_rez[j][k]) > (adj_beta_merg_cnt * float(cores_count[j]))) or (float(same_sheet_rez[k][j]) > (adj_beta_merg_cnt * float(cores_count[k]))))
                        if not selection:
                            selection = (( (((complex_entanglement_matrix[j][k] >= 2 * min_contacts_for_merging_entangled_domains)
                                          or ((float(entanglement_matrix[j][k]) / float(cores_count[j])) >= 2 * min_entanglement_for_merging_entangled_domains))
                                         or ((complex_entanglement_matrix[k][j] >= 2 * min_contacts_for_merging_entangled_domains)
                                          or ((float(entanglement_matrix[k][j]) / float(cores_count[k])) >= 2 * min_entanglement_for_merging_entangled_domains))))
                                     and ((float(same_sheet_rez[j][k]) > float(cores_count[j])) or (float(same_sheet_rez[k][j]) > float(cores_count[k]))))
                    if selection:
                        # one domain is contained within another, but they are highly contacted. merge them
                        reset_button = 1

                        # first, merge these two cores into the index of the first
                        if j<k:
                            for m in range(cores_count[k]):
                                cores[j][cores_count[j]+m] = cores[k][m]
                                cores[k][m] = 0
                            # join their counts
                            cores_count[j] += cores_count[k]
                            cores_count[k] = 0
                        else:
                            for m in range(cores_count[j]):
                                cores[k][cores_count[k]+m] = cores[j][m]
                                cores[j][m] = 0
                            # join their counts
                            cores_count[k] += cores_count[j]
                            cores_count[j] = 0
                        if j<k:
                            # next, copy the rest of the cores back to fill in the gap
                            for m in range(number_of_cores):
                                for n in range(cores_count[m]):
                                    cores[m-1][n] = cores[m][n]
                                    cores[m][n] = 0
                                cores_count[m-1] = cores_count[m]
                                cores_count[m] = 0
                            cores_count[number_of_cores-1] = 0
                            number_of_cores -= 1
                        else:
                            # next, copy the rest of the cores back to fill in the gap
                            for m in range(j+1, number_of_cores):
                                for n in range(cores_count[m]):
                                    cores[m-1][n] = cores[m][n]
                                    cores[m][n] = 0
                                cores_count[m-1] = cores_count[m]
                                cores_count[m] = 0
                            cores_count[number_of_cores-1] = 0
                            number_of_cores -= 1
                        # first recalculate ranges for the new cores
                        for p in range(number_of_cores):
                            ranges[p][0] = 0
                            ranges[p][1] = 0
                            start_searching = int(math.ceil(float(cores_count[p]) * ((1.0 - blocked_range_percentage) / 2.0)))
                            # first find the lower end
                            blocked = 0
                            for m in range(mstln):
                                for n in range(cores_count[p]):
                                    if cores[p][n] == m:                              # if this one is the mth one
                                        if blocked >= start_searching:                  # if were already past the blocked residues
                                            ranges[p][0] = m                            # this one starts at m
                                            break                                       # finished with the mth residue
                                        blocked += 1
                                if ranges[p][0]:                                      # if the lower range has been found, break
                                    break
                            blocked = 0
                            # now find the upper end
                            for m in range(mstln-1, -1, -1):
                                for n in range(cores_count[p]):
                                    if cores[p][n] == m:
                                        if blocked >= start_searching:
                                            ranges[p][1] = m
                                            break
                                        blocked += 1
                                if ranges[p][1]:
                                    break

                        # now recalculate the minimum distances between alpha carbons of the remaining cores
                        for m in range(number_of_cores):
                            for n in range(number_of_cores):
                                min_alpha_distances[m][n] = 10000.0
                        # now calculate the distance of the closest alpha carbons between the two cores.
                        for q in range(number_of_cores):
                            for r in range(number_of_cores):
                                if r == q:
                                    continue
                                else:
                                    for m in range(cores_count[q]):
                                        for n in range(cores_count[r]):
                                            p1x = pchain.residues[cores[q][m]].x
                                            p1y = pchain.residues[cores[q][m]].y
                                            p1z = pchain.residues[cores[q][m]].z
                                            p2x = pchain.residues[cores[r][n]].x
                                            p2y = pchain.residues[cores[r][n]].y
                                            p2z = pchain.residues[cores[r][n]].z
                                            d = math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2) 
                                            if d < min_alpha_distances[q][r]:
                                                min_alpha_distances[q][r] = d
                                                min_alpha_distances[r][q] = min_alpha_distances[q][r]

                        # now recalculate the number of alpha carbon atoms within a threshold distance
                        # for (int m=0; m<number_of_cores; m++){
                        #   close_alphas[m] = new int[number_of_cores]
                        # }
                        for m in range(number_of_cores):
                            for n in range(number_of_cores):
                                close_alphas[m][n] = 0

                        for q in range(number_of_cores):
                            for r in range(number_of_cores):
                                if r == q:
                                    continue
                                else:
                                    for m in range(cores_count[q]):
                                        for n in range(cores_count[r]):
                                            p1x = pchain.residues[cores[q][m]].x
                                            p1y = pchain.residues[cores[q][m]].y
                                            p1z = pchain.residues[cores[q][m]].z
                                            p2x = pchain.residues[cores[r][n]].x
                                            p2y = pchain.residues[cores[r][n]].y
                                            p2z = pchain.residues[cores[r][n]].z
                                            if math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2) < alpha_merge_distance:
                                                close_alphas[q][r] += 1
                        for z in range(number_of_cores):
                            for y in range(number_of_cores):
                                same_sheet_rez[z][y] = 0
                        """
                        for z in range(number_of_cores):
                            for y in range(number_of_cores):
                                if y==z:
                                    continue
                                else:
                                    # for each core residue of the first domain, find out which sheet its in.
                                    # then step through the core residues of the second domain. if they are in the same sheet, increment same_sheet_rez
                                    owner_sheet = 0
                                    for m in range(cores_count[z]):
                                        if pchain.residues.features['secondary'] == 'B':
                                            for n in range(pchain.sheets.p[i].FragmentSetCnt):
                                                for p in range(pchain.sheets.p[i].s[n].FragCnt):        # for each strand in this sheet
                                                    for q in range(pchain.sheets.p[i].s[n].f[p].PosCnt):
                                                        if cores[z][m] == pchain.sheets.p[i].s[n].f[p].p[q]:
                                                            owner_sheet = n
                                            for n in range(cores_count[y]):
                                                if pchain.residues[cores[y][n]].features['secondary'] == 'B':
                                                    for p in range(pchain.sheets.p[i].FragmentSetCnt):
                                                        for q in range(pchain.sheets.p[i].s[p].FragCnt):
                                                            for r in range(pchain.sheets.p[i].s[p].f[q].PosCnt):
                                                                if cores[y][n] == pchain.sheets.p[i].s[p].f[q].p[r]:
                                                                    if owner_sheet == p:
                                                                        same_sheet_rez[z][y] += 1
                        """
                        # now recalculate the interactions and entanglements
                        for m in range(number_of_cores):
                            for n in range(number_of_cores):
                                interactions_matrix[m][n] = 0
                        for q in range(number_of_cores):
                            for r in range(q+1, number_of_cores):
                                interaction = 0
                                for m in range(cores_count[q]):                            # for each contact of the first
                                    for n in range(len(pchain.residues[cores[q][m]].neighbors)):
                                        for p in range(cores_count[r]):
                                            if pchain.residues[cores[q][m]].neighbors[n] == cores[r][p]:
                                                # check the distances
                                                p1x = pchain.residues[cores[q][m]].pseudo_sidechain.x
                                                p1y = pchain.residues[cores[q][m]].pseudo_sidechain.y
                                                p1z = pchain.residues[cores[q][m]].pseudo_sidechain.z
                                                p2x = pchain.residues[cores[r][p]].pseudo_sidechain.x
                                                p2y = pchain.residues[cores[r][p]].pseudo_sidechain.y
                                                p2z = pchain.residues[cores[r][p]].pseudo_sidechain.z
                                                if math.sqrt((p1x-p2x)**2+(p1y-p2y)**2+(p1z-p2z)**2) < core_ext_distance:
                                                    interaction += 1                                          # then tally the observation
                                                    interactions_matrix[q][r] += 1
                                                    interactions_matrix[q][r] += 1
                        for q in range(number_of_cores):
                            for r in range(number_of_cores):
                                entanglement_matrix[q][r] = 0
                        for q in range(number_of_cores):
                            for r in range(number_of_cores):
                                if r == q:
                                    continue
                                else:
                                    for m in range(ranges[q][0], ranges[q][1]+1):
                                        if m>ranges[r][0] and m<ranges[r][1]:
                                            for p in range(cores_count[q]):
                                                if cores[q][p] == m:
                                                    entanglement_matrix[q][r] += 1
                        for q in range(number_of_cores):
                            for r in range(number_of_cores):
                                complex_entanglement_matrix[q][r] = 0
                        for q in range(number_of_cores):
                            for r in range(number_of_cores):
                                if r == q:
                                    continue
                                else:
                                    for m in range(ranges[q][0], ranges[q][1]+1):
                                        for n in range(cores_count[r]):
                                            for p in range(len(pchain.residues[cores[r][n]].neighbors)):
                                                if pchain.residues[cores[r][n]].neighbors[p] == m:
                                                    complex_entanglement_matrix[q][r] += 1
                if reset_button:
                    break

            # now re-write the full_core_list array to pass the display over to dino
            for q in range(mstln):
                full_core_list[q] = 0
            k_ctr = -1
            for q in range(number_of_cores):
                for r in range(cores_count[q]):
                    k_ctr += 1
                    full_core_list[k_ctr] = cores[q][r]
                k_ctr += 1
                full_core_list[k_ctr] = -1
  
            if reset_button:
                break
        
        if not reset_button and selection_choice in [1,5]:
            selection_choice = 2
        elif not reset_button and selection_choice == 6:
            selection_choice = 7
        elif not reset_button and selection_choice == 8:
            selection_choice = 9
        elif not reset_button and selection_choice == 9:
            break
            selection_choice = 3
        elif not reset_button and selection_choice in [0,2,3,4,7]:
            break


    # now re-write the full_core_list array to pass the display over to dino
    for j in range(mstln):
        full_core_list[j] = 0
    k_ctr = -1
    for j in range(number_of_cores):
        for k in range(cores_count[j]):
            k_ctr += 1
            full_core_list[k_ctr] = cores[j][k]
        k_ctr += 1
        full_core_list[k_ctr] = -1
  
    # now print out the cores using ranges. a range is here defined for each segment of the protein that 
    #  contains more than minimum_core_rez_for_ranges residues that are contiguously of one core.

    # first make an array of bools for each residue in the mstln; set all non-gap positions to 1
    segments_available_cnt = 0
    core_types, segments_available = [],[]
    for n in range(mstln):
        core_types.append(0)
        segments_available.append(0)
    for j in range(mstln):
        # if this residue is in one of the cores set to 1, else 0
        segments_available[j] = 0
        core_types[j] = 0
        for m in range(number_of_cores):
            for n in range(cores_count[m]):
                if cores[m][n] == j:
                    segments_available[j] = 1
                    segments_available_cnt += 1
                    core_types[j] = m

    # now evaluate the selections
    # first compare the number of cores predicted vs. the number of cores observed.
    if multi_domain > 1:
        domain_counter = 0
        for j in range(multi_domain):               # for each multi_domain - the number of domains observed for this protein
            for k in range(contained_domain+1):
                if gold[0][domain_counter] != 0:
                    domain_counter += 1                            # increment the domain_counter

    # now assign the rest of the residues
    # use a window to get a majority vote
    # determine to which core the closest 5 core residues belong, and assign 'this' residue to the majority

    # first make an array of ints with core assignments, to make it all easier
    window_size = 7
    linear_core_assignments = []
    final_core_assignments  = []
    for n in range(mstln):
        linear_core_assignments.append(-1)                    # -1 represents default - not in a core
        final_core_assignments.append(-1)
        
    for j in range(mstln):
        for m in range(number_of_cores):
            for n in range(cores_count[m]):
                if cores[m][n] == j:
                    linear_core_assignments[j] = m
                    final_core_assignments[j] = m

    # now make the actual assignments
    my_closest_cores = []
    for n in range(number_of_cores):
        my_closest_cores.append(0)
    if number_of_cores > 1:
        for j in range(mstln):
          got_this_many = 0
          if linear_core_assignments[j] == -1:                   # if -1, then it needs assigned
              for k in range(number_of_cores):
                  my_closest_cores[k] = -1
              increment_search = 1
              top_of_search = j
              bot_of_search = j
              while 1:                                             # add residues to the neighbors
                  while 1:                                           # first find the next forward residue that belongs to a core
                      top_of_search += 1
                      if top_of_search >= mstln:                       # but break if past the end of the array
                          break
                      elif linear_core_assignments[top_of_search] > -1:
                          break
                  if top_of_search < mstln:                          # if < mstln, then it broke on a core residue
                      my_closest_cores[linear_core_assignments[top_of_search]] += 1    # add one to this observation
                      got_this_many += 1
                  if got_this_many == window_size:
                      break
                  while 1:                                           # first find the next forward residue that belongs to a core
                      bot_of_search -= 1
                      if bot_of_search < 0:                            # but break if past the end of the array
                          break
                      elif linear_core_assignments[bot_of_search] > -1:
                          break
                  if bot_of_search >= 0:                             # if >= 0, then it broke on a core residue
                      my_closest_cores[linear_core_assignments[bot_of_search]] += 1
                      got_this_many += 1
                  if got_this_many == window_size:
                      break
              # now assign it to the majority core
              majority_core = -1
              majority_core_count = 0
              second_core = -1
              second_core_count = 0
              for k in range(number_of_cores):
                  if my_closest_cores[k] > majority_core_count:
                      second_core_count = majority_core_count
                      second_core = majority_core
                      majority_core_count = my_closest_cores[k]
                      majority_core = k
              if second_core_count == majority_core_count:          # if its a tie, randomly choose one
                  if math.random()*2:
                      majority_core = second_core
              final_core_assignments[j] = majority_core
    else:
        for j in range(mstln):
            if final_core_assignments[j] != -2:
                final_core_assignments[j] = 0

    # now do another round of smoothing, where each residue is assigned to the core of the surrounding 2 residues
    first_back_index = -1
    second_back_index = -1
    for j in range(mstln):
        if second_back_index != -1 and first_back_index != -1:
            if (final_core_assignments[j] == final_core_assignments[second_back_index]) or (final_core_assignments[j] == final_core_assignments[first_back_index]):
                final_core_assignments[first_back_index] = final_core_assignments[j]
            elif final_core_assignments[first_back_index] == final_core_assignments[second_back_index]:
                final_core_assignments[first_back_index] = final_core_assignments[first_back_index]
            else:
                if random.random()*2:
                    final_core_assignments[first_back_index] = final_core_assignments[second_back_index]
                else:
                    final_core_assignments[first_back_index] = final_core_assignments[j]
        second_back_index = first_back_index
        first_back_index = j

    return float(number_of_cores)

