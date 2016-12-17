import Numeric
import math
import sys
import copy
sys.path.append('./Tools/Superimpose')
import Superimpose
sys.path.append('./Tools/Aligner')
import Aligner

class SequenceAligner(Aligner.Aligner):
    def __init__(self, print_param=0, alignment_type='local', gap_scores=[0.2,0.05]):
        Aligner.Aligner.__init__(self, print_param, alignment_type, gap_scores)
        # PET91 
        self.scoring_matrix =  {'AA':15/15.0, 'A-':6/15.0,
                                'AR': 7/15.0, 'RR':15/15.0, 'R-':6/15.0,
                                'AN': 9/15.0, 'RN': 8/15.0, 'NN':15/15.0, 'N-':6/15.0,
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
                                'AV':11/15.0, 'RV': 4/15.0, 'NV': 5/15.0, 'DV': 5/15.0, 'CV': 7/15.0, 'QV': 4/15.0, 'EV': 5/15.0, 'GV': 6/15.0, 'HV': 4/15.0, 'IV':15/15.0, 'LV':11/15.0, 'KV': 4/15.0, 'MV':12/15.0, 'FV': 8/15.0, 'PV': 7/15.0, 'SV': 7/15.0, 'TV':10/15.0, 'WV': 4/15.0, 'YV': 4/15.0, 'VV':15/15.0}
        for key in self.scoring_matrix.keys():
            if key[1]+key[0] not in self.scoring_matrix.keys():
                self.scoring_matrix[key[1]+key[0]] = self.scoring_matrix[key[0]+key[1]]

    def score(self, r, m):
        return self.scoring_matrix[self.target.residues[r].res_type1+self.template.residues[m].res_type1]
        
    def align_sequences(self):
        aout, bout = self.align()
        
        identical_pairs = 0
        total_pairs     = 0
        for ind in range(len(aout)):
            if aout[ind] != '-' and bout[ind] != '-':
                if aout[ind].res_type1 == bout[ind].res_type1:
                    identical_pairs += 1
            total_pairs += 1
        return (identical_pairs + 0.0)/total_pairs

    def get_alignment(self):
        aout, bout = self.align()
        
        identical_pairs = 0
        total_pairs     = 0
        for ind in range(len(aout)):
            if aout[ind] != '-' and bout[ind] != '-':
                if aout[ind].res_type1 == bout[ind].res_type1:
                    identical_pairs += 1
            total_pairs += 1
        print '%5.3f percent identity'%((identical_pairs + 0.0)/total_pairs)
        return aout, bout
        

