import string
import sys
sys.path.append('./Tools/Aligner')
import SequenceAligner
import math

# normalized PET91, gap score 6
scoring_matrix =  {'A-':6/15.0,
                   'AA':15/15.0, 'R-':6/15.0,
                   'AR': 7/15.0, 'RR':15/15.0, 'N-':6/15.0,
                   'AN': 9/15.0, 'RN': 8/15.0, 'NN':15/15.0, 'D-':6/15.0,
                   'AD': 8/15.0, 'RD': 6/15.0, 'ND':13/15.0, 'DD':15/15.0, 'C-': 6/15.0,
                   'AC': 7/15.0, 'RC': 8/15.0, 'NC': 7/15.0, 'DC': 3/15.0, 'CC':15/15.0, 'Q-': 6/15.0,
                   'AQ': 7/15.0, 'RQ':12/15.0, 'NQ': 9/15.0, 'DQ': 9/15.0, 'CQ': 4/15.0, 'QQ':15/15.0, 'E-': 6/15.0,
                   'AE': 8/15.0, 'RE': 7/15.0, 'NE': 9/15.0, 'DE':15/15.0, 'CE': 2/15.0, 'QE':12/15.0, 'EE':15/15.0, 'G-': 6/15.0,
                   'AG':10/15.0, 'RG': 9/15.0, 'NG': 9/15.0, 'DG':10/15.0, 'CG': 7/15.0, 'QG': 6/15.0, 'EG': 9/15.0, 'GG':15/15.0, 'H-': 6/15.0,
                   'AH': 6/15.0, 'RH':12/15.0, 'NH':12/15.0, 'DH': 9/15.0, 'CH': 8/15.0, 'QH':14/15.0, 'EH': 7/15.0, 'GH': 6/15.0, 'HH':15/15.0, 'I-': 6/15.0,
                   'AI': 9/15.0, 'RI': 4/15.0, 'NI': 6/15.0, 'DI': 3/15.0, 'CI': 5/15.0, 'QI': 4/15.0, 'EI': 3/15.0, 'GI': 4/15.0, 'HI': 4/15.0, 'II':15/15.0, 'L-': 6/15.0,
                   'AL': 6/15.0, 'RL': 5/15.0, 'NL': 4/15.0, 'DL': 2/15.0, 'CL': 5/15.0, 'QL': 7/15.0, 'EL': 2/15.0, 'GL': 2/15.0, 'HL': 6/15.0, 'IL':12/15.0, 'LL':15/15.0, 'K-': 6/15.0,
                   'AK': 6/15.0, 'RK':14/15.0, 'NK':11/15.0, 'DK':12/15.0, 'CK': 4/15.0, 'QK':12/15.0, 'EK':10/15.0, 'GK': 6/15.0, 'HK': 9/15.0, 'IK': 4/15.0, 'LK': 4/15.0, 'KK':15/15.0, 'M-': 6/15.0,
                   'AM': 8/15.0, 'RM': 6/15.0, 'NM': 5/15.0, 'DM': 3/15.0, 'CM': 5/15.0, 'QM': 6/15.0, 'EM': 4/15.0, 'GM': 4/15.0, 'HM': 5/15.0, 'IM':13/15.0, 'LM':13/15.0, 'KM': 6/15.0, 'MM':15/15.0, 'F-': 6/15.0,
                   'AF': 4/15.0, 'RF': 2/15.0, 'NF': 4/15.0, 'DF': 1/15.0, 'CF': 8/15.0, 'QF': 3/15.0, 'EF': 0/15.0, 'GF': 1/15.0, 'HF': 8/15.0, 'IF': 9/15.0, 'LF':12/15.0, 'KF': 1/15.0, 'MF': 8/15.0, 'FF':15/15.0, 'P-': 6/15.0,
                   'AP':10/15.0, 'RP': 8/15.0, 'NP': 7/15.0, 'DP': 5/15.0, 'CP': 5/15.0, 'QP':10/15.0, 'EP': 5/15.0, 'GP': 6/15.0, 'HP': 9/15.0, 'IP': 5/15.0, 'LP': 8/15.0, 'KP': 6/15.0, 'MP': 5/15.0, 'FP': 5/15.0, 'PP':15/15.0, 'S-': 6/15.0,
                   'AS':12/15.0, 'RS': 8/15.0, 'NS':12/15.0, 'DS': 8/15.0, 'CS':10/15.0, 'QS': 7/15.0, 'ES': 7/15.0, 'GS':10/15.0, 'HS': 8/15.0, 'IS': 7/15.0, 'LS': 7/15.0, 'KS': 7/15.0, 'MS': 7/15.0, 'FS': 7/15.0, 'PS':11/15.0, 'SS':15/15.0, 'T-': 6/15.0,
                   'AT':12/15.0, 'RT': 7/15.0, 'NT':11/15.0, 'DT': 7/15.0, 'CT': 7/15.0, 'QT': 7/15.0, 'ET': 6/15.0, 'GT': 8/15.0, 'HT': 7/15.0, 'IT':10/15.0, 'LT': 6/15.0, 'KT': 8/15.0, 'MT':10/15.0, 'FT': 4/15.0, 'PT':10/15.0, 'ST':12/15.0, 'TT':15/15.0, 'W-': 6/15.0,
                   'AW': 2/15.0, 'RW': 9/15.0, 'NW': 2/15.0, 'DW': 1/15.0, 'CW':10/15.0, 'QW': 4/15.0, 'EW': 1/15.0, 'GW': 7/15.0, 'HW': 4/15.0, 'IW': 3/15.0, 'LW': 6/15.0, 'KW': 4/15.0, 'MW': 4/15.0, 'FW': 7/15.0, 'PW': 2/15.0, 'SW': 5/15.0, 'TW': 3/15.0, 'WW':15/15.0, 'Y-': 6/15.0,
                   'AY': 3/15.0, 'RY': 5/15.0, 'NY': 8/15.0, 'DY': 6/15.0, 'CY':12/15.0, 'QY': 6/15.0, 'EY': 3/15.0, 'GY': 2/15.0, 'HY':14/15.0, 'IY': 5/15.0, 'LY': 6/15.0, 'KY': 3/15.0, 'MY': 4/15.0, 'FY':15/15.0, 'PY': 3/15.0, 'SY': 7/15.0, 'TY': 4/15.0, 'WY': 8/15.0, 'YY':15/15.0, 'V-': 6/15.0,
                   'AV':11/15.0, 'RV': 4/15.0, 'NV': 5/15.0, 'DV': 5/15.0, 'CV': 7/15.0, 'QV': 4/15.0, 'EV': 5/15.0, 'GV': 6/15.0, 'HV': 4/15.0, 'IV':15/15.0, 'LV':11/15.0, 'KV': 4/15.0, 'MV':12/15.0, 'FV': 8/15.0, 'PV': 7/15.0, 'SV': 7/15.0, 'TV':10/15.0, 'WV': 4/15.0, 'YV': 4/15.0, 'VV':15/15.0, '--':11/15.0}
for key in scoring_matrix.keys():
    if key[1]+key[0] not in scoring_matrix.keys():
        scoring_matrix[key[1]+key[0]] = scoring_matrix[key[0]+key[1]]

def fetch_msq_conservation(pchain):
    """ shortcut to apply .msq files as alignments: probably not a permenant function """
    filename = pchain.parent.get_filename_by_extension('.msq', pchain.chain_name)
    try:
        contact_file = open(filename)
    except IOError:
        print 'no .msq file present for chain %s'%(pchain.chain_name)
        return 0

    lines = contact_file.readlines()
    firstline = ""
    for line in lines:
        if line[0] != '-':
            firstline += line[0]

    a = SequenceAligner.SequenceAligner(1, 'global')
    a.add_target(pchain)
    a.add_template_sequence(firstline)
    target_out, template_out = a.get_alignment()

    template_index = 0
    target_index = 0
    for i in range(len(target_out)):
        if target_out[i] == '-' and template_out[i] != '-':
            # gap in the target... advance template
            print "\n\n!!!  MISALIGNMENT  !!!\n\n"
            template_index += 1
        elif target_out[i] != '-' and template_out[i] == '-':
            # gap in the template... advance target
            print "\n\n!!!  MISALIGNMENT  !!!\n\n"
            pchain.residues[target_index].data['conservation'] = ""
            target_index += 1
        elif target_out[i] != '-' and template_out[i] != '-':
            pchain.residues[target_index].data['conservation'] = string.strip(lines[template_index])
            template_index += 1
            target_index += 1
        else:
            print 'error'
    return 1 # success -- found an msq file

def apply_sequence_alignment(system, sequences):
    """ given a set of aligned sequences, where the top corresponds to the query,
        fill conservation features with strings from columns of the input table.
    """
    first_sequence = string.upper(sequences[0])
    match_sequence = ""
    for c in first_sequence:
        if c in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
            match_sequence += c
            
    for pchain in system.ProteinList:
        s = pchain.get_sequence()
        if s != match_sequence:
            print 'no match for chain %s'%(pchain.chain_name)
            continue
        else:
            print 'match ok for chain %s'%(pchain.chain_name)
            
        rs_index = 0
        #print '***' + len(first_sequence), first_sequence
        #for seq in sequences:
        #    print len(seq), seq
        
        for fs_ind in range(len(first_sequence)):
            k = string.upper(first_sequence[fs_ind])
            #print fs_ind, len(first_sequence), k
            if k in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
                conservation = ''
                for seq in sequences:
                    if seq[fs_ind] not in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
                        conservation += '-'
                    else:
                        conservation += seq[fs_ind]
                pchain.residues[rs_index].data['conservation'] = conservation
                rs_index += 1
        
        alignment_sequences = []
        for c in pchain.residues[0].data['conservation']:
            alignment_sequences.append("")
        for res in pchain.residues:
            ind = 0
            for c in res.data['conservation']:
                alignment_sequences[ind] += c
                ind += 1
        pchain.data['sequence_alignment'] = alignment_sequences
    
def print_sequence_alignment(system):
    for pchain in system.ProteinList:
        try:
            pchain.data['sequence_alignment']
        except KeyError:
            pass
        else:
            for sequence in pchain.data['sequence_alignment']:
                print sequence
            print pchain.data['sequence_alignment']
            for res in pchain.residues:
                print res.data['conservation']

def calculate_conservation(system, asa_style='sidechain_asa', rewrite=0, viewer=None):
    system.calculate_differential_system_asa(1.4, 1000, rewrite)

    for pchain in system.ProteinList:
        try:
            pchain.residues[0].data['conservation']
        except KeyError:
            continue

        use_asa = 0
        for res in pchain.residues:
            try:
                res.features[asa_style]
            except KeyError:
                print 'no solvent accessibilities calculated, aborting conservation calculations'
                break
        else:
            use_asa = 1

        cons_test = 1
        for res in pchain.residues:
            try:
                res.data['conservation']
            except KeyError:
                print 'fetching conservation from msqs'
                cons_test = fetch_msq_conservation(pchain)
                break

        if not cons_test:
            print 'no alignment file found for chain %s'%(pchain.chain_name)
            continue

        for res in pchain.residues:             # test for shielding calculations
            try:
                res.data['neighbors']
            except KeyError:
                pchain.fill_pseudo_sidechains(0)
                pchain.fill_neighbors_lists(0.35, 6.0)
                break

        # calculate 0D conservation
        for res in pchain.residues:
            cons_string = res.data['conservation']
            if len(cons_string) > 0:
                score = 0.0
                count = 0.0
                for i1 in range(len(cons_string)-1):
                    for i2 in range(i1+1, len(cons_string)):
                        key = '%s%s'%(cons_string[i1],cons_string[i2])
                        try:
                            score += scoring_matrix[key]
                        except KeyError:
                            score += scoring_matrix['A-']
                        count += 1
                res.data['0D_conservation'] = score / count
                res.data['noactsit_0D_conservation'] = score / count
            else:
                res.data['0D_conservation'] = None

        # calculate the 0D normalized scores here so they can be used for filters later
        minval,maxval=10.0,0.0
        for res in pchain.residues:
            if use_asa and res.features[asa_style] < 0.05:
                continue
            f = res.data['0D_conservation']
            if f < minval:
                minval = f
            elif f > maxval:
                maxval = f
        for res in pchain.residues:
            if use_asa and res.features[asa_style] < 0.05:
                res.features['normalized_0D_conservation'] = -1
                continue
            res.features['normalized_0D_conservation'] = (res.data['0D_conservation']-minval)/(maxval-minval)

        # now that conservation has been calculated
        _identify_low_conservation_loops(pchain)

        # calculate 1D conservation
        for i in range(len(pchain.residues)):
            sum, count = _get_1D_sum_and_count(i, pchain)
            pchain.residues[i].data['1D_conservation'] = sum/count

        # 3D conservation
        for i in range(len(pchain.residues)):
            sum, count = _get_3D_sum_and_count(i, pchain, asa_style)
            pchain.residues[i].data['3D_conservation'] = sum/count
            
        # calculate 3D shielding
        i = 0
        for res in pchain.residues:
            sum, count = _get_1D_sum_and_count(i, pchain)
            sum2,count2= _get_s3D_sum_and_count(i, pchain, asa_style, use_asa)
            res.data['ms3D_conservation'] = (sum+sum2)/(count+count2)
            i += 1

        # calculate filtered 3D shielding
        # create a list of atoms with ligands to filter out active sites
        ligatoms = []
        for lig in system.LigandList:
            for atom in lig.atoms:
                ligatoms.append(atom)

        # calculate 3D shielding
        i = 0
        for res in pchain.residues:
            sum, count = _get_1D_sum_and_count(i, pchain)
            sum2,count2= _get_noactsit_s3D_sum_and_count(i, pchain, ligatoms, asa_style, use_asa)
            sum3,count3= _get_nobadloop_s3D_sum_and_count(i, pchain, asa_style, use_asa)
            res.data['noactsit_ms3D_conservation'] = (sum+sum2)/(count+count2)
            res.data['nobadloop_ms3D_conservation'] = (sum+sum3)/(count+count3)
            i += 1
                
        # offer normalized versions for all scores
        # calculate 3D shielding with and without active sites                
        #types = ['0D_conservation', '1D_conservation', '3D_conservation', 'ms3D_conservation']
        types = ['noactsit_0D_conservation', '1D_conservation', '3D_conservation', 'ms3D_conservation', 'noactsit_ms3D_conservation', 'nobadloop_ms3D_conservation']
        for type in types:
            print 'calculating conservation %s'%(type)
            minval,maxval=10.0,0.0
            # first find the range
            for res in pchain.residues:
                if use_asa and res.features[asa_style] < 0.05:
                    continue
                if type in ['noactsit_ms3D_conservation', 'noactsit_0D_conservation']:
                    lighit = 0
                    for resatom in res.atoms:
                        for ligatom in ligatoms:
                            if ligatom.dist(resatom) < 5.0:
                                lighit=1
                    if lighit:
                        continue
                if type == 'nobadloop_ms3D_conservation':
                    if res.data['loop_status']:
                        conssum = 0.0
                        conscnt = 0.0
                        for rex in res.data['loop_status']:
                            conssum += pchain.residues[rex].features['normalized_0D_conservation']
                            conscnt += 1.0
                        if (conssum / conscnt) < 0.25:
                            continue
                    
                f = res.data[type]
                if f < minval:
                    minval = f
                elif f > maxval:
                    maxval = f
            str = 'normalized_' + type
            res_index = 0
            
            # now normalize by it
            for res in pchain.residues:
                if use_asa and res.features[asa_style] < 0.05:
                    res.features[str] = -1
                    continue
                if type in ['noactsit_ms3D_conservation', 'noactsit_0D_conservation']:
                    lighit = 0
                    for resatom in res.atoms:
                        for ligatom in ligatoms:
                            if ligatom.dist(resatom) < 5.0:
                                lighit = 1
                    if lighit:
                        res.features[str] = -1
                        continue
                if type == 'nobadloop_ms3D_conservation':
                    if res.data['loop_status']:
                        conssum = 0.0
                        conscnt = 0.0
                        for rex in res.data['loop_status']:
                            conssum += pchain.residues[rex].features['normalized_0D_conservation']
                            conscnt += 1.0
                        if (conssum / conscnt) < 0.25:
                            res.features[str] = -1
                            continue
                res.features[str] = (res.data[type]-minval)/(maxval-minval)
                res_index += 1

        """
        # the following block creates Z-scores for the conservation values and is not necessary
        types = ['0D_conservation', 'noactsit_0D_conservation', '1D_conservation', '3D_conservation', 'ms3D_conservation', 'noactsit_ms3D_conservation', 'nobadloop_ms3D_conservation']
        for type in types:
            str = 'normalized_' + type
            # calculate the average
            sum = 0.0
            count = 0.0
            for res in pchain.residues:
                if res.features[str] and res.features[str] != -1:
                    sum += res.features[str]
                    count += 1.0
            ave = sum/count
            # and the SD
            sum = 0.0
            for res in pchain.residues:
                if res.features[str] and res.features[str] != -1:
                    sum += (ave-res.features[str])**2
            sd = math.sqrt(sum/(count-1.0))
            # now transform each of the values
            for res in pchain.residues:
                if res.features[str] and res.features[str] != -1:
                    t = res.features[str]
                    res.features[str] = (t-ave)/sd
            """

        if viewer != None:
            viewer.rebuild_color_menu()

def _identify_low_conservation_loops(pchain):
    pchain.fill_secondary()
    # look for loops at least three amino acids long
    len_so_far = 0
    start_res   = 0
    for rex in range(len(pchain.residues)):
        res = pchain.residues[rex]
        res.data['loop_status'] = 0
        if len_so_far == 0:
            start_res = rex
        if res.features['secondary'] == 'C':
            len_so_far += 1
        else:
            if len_so_far >= 3:
                loop_list = []
                for i in range(start_res, rex):
                    loop_list.append(i)
                for i in range(start_res, rex):
                    pchain.residues[i].data['loop_status'] = loop_list
            len_so_far = 0

def _get_1D_sum_and_count(i, pchain):
    calset = {0:[0,1,2,3],
              1:[-1,0,1,2,3],
              2:[-2,-1,0,1,2,3],
              (len(pchain.residues)-3):[-3,-2,-1,0,1,2],
              (len(pchain.residues)-2):[-3,-2,-1,0,1],
              (len(pchain.residues)-1):[-3,-2,-1,0]}
    if i in calset.keys():
        set = calset[i]
    else:
        set = [-3,-2,-1,0,1,2,3]
    sum = 0.0
    count = 0
    for offset in set:
        count += 1
        sum += pchain.residues[i+offset].data['0D_conservation']
    return sum, count

def _get_3D_sum_and_count(i, pchain, asa_style='sidechain_asa', exposure=None):
    sum = pchain.residues[i].data['0D_conservation']
    count = 1
    c1 = pchain.residues[i].central_atom
    for ind in range(len(pchain.residues)):
        if i == ind:
            continue
        c2 = pchain.residues[ind].central_atom
        if c1.dist(c2) > 6.0:
            continue
        if exposure and pchain.residues[ind].features[asa_style] < 0.05:
            continue
        else:
            count += 1
            sum += pchain.residues[ind].data['0D_conservation']
    return sum, count
    
def _get_s3D_sum_and_count(i, pchain, asa_style='sidechain_asa', exposure=None):
    sum = pchain.residues[i].data['0D_conservation']
    count = 1
    for ind in pchain.residues[i].neighbors:
        if exposure and pchain.residues[ind].features[asa_style] < 0.05:
            continue
        else:
            count += 1
            sum += pchain.residues[ind].data['0D_conservation']
    return sum, count

def _get_noactsit_s3D_sum_and_count(i, pchain, ligatoms, asa_style='sidechain_asa', exposure=None):
    sum = pchain.residues[i].data['0D_conservation']
    count = 1
    for ind in pchain.residues[i].neighbors:
        if exposure and pchain.residues[ind].features[asa_style] < 0.05:
            continue
        lighit = 0
        for atom in pchain.residues[ind].atoms:
            for ligatom in ligatoms:
                if ligatom.dist(atom) < 5.0:
                    lighit = 1
        if lighit:
            continue
            
        count += 1
        sum += pchain.residues[ind].data['0D_conservation']
    return sum, count


def _get_nobadloop_s3D_sum_and_count(i, pchain, ligatoms, asa_style='sidechain_asa', exposure=None):
    sum = pchain.residues[i].data['0D_conservation']
    count = 1
    for ind in pchain.residues[i].neighbors:
        if exposure and pchain.residues[ind].features[asa_style] < 0.05:
            continue
        res = pchain.residues[i]
        # if is in a loop that is low conservation, skip
        if res.data['loop_status']:
            conssum = 0.0
            conscnt = 0.0
            for rex in res.data['loop_status']:
                conssum += pchain.residues[rex].features['normalized_0D_conservation']
                conscnt += 1.0
            if (conssum / conscnt) < 0.25:
                continue
        count += 1
        sum += pchain.residues[ind].data['0D_conservation']
    return sum, count

