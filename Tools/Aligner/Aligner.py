import numpy
import math

class Aligner:
    """ Base Class Only: Implements a versitile dynamic programming algorithm """
    def __init__(self, print_param=0, alignment_type='local', gap_scores=[0.8,0.05]):
        """ May be alignment_types 'Needleman-Wunsch', 'Semiglobal', 'Smith-Waterman'
            alternatively,         'global',           'semiglobal', 'local'
        """
        # use global, semiglobal, local
        if alignment_type == 'Smith-Waterman':
            alignment_type = 'local'
        elif alignment_type == 'Needleman-Wunsch':
            alignment_type = 'global'
        elif alignment_type == 'Semiglobal':
            alignment_type = 'semiglobal'
        self.alignment_type = alignment_type
        self.print_param = print_param
        self.gap_score = gap_scores[0]
        self.gap_extension_score = gap_scores[1]
        
    def add_target(self, polymer):
        """ aligns .central_atoms from passed in polymer objects """
        # this style of function should keep flexibility for doing multiple alignments 
        # turn .target into a list and append polymers to avoid changing the interface
        self.target = self._prepare_polymer(polymer)
        
    def add_template(self, polymer):
        """ aligns .central_atoms from passed in polymer objects """
        self.template = self._prepare_polymer(polymer)

    def add_target_sequence(self, sequence):
        self.add_target(self._prepare_fake_polymer(sequence))
    
    def add_template_sequence(self, sequence):
        self.add_template(self._prepare_fake_polymer(sequence))
    
    def _prepare_polymer(self, polymer):
        return polymer

    def _prepare_fake_polymer(self, sequence):
        # need a class that holds .residues that support .res_type1
        class fake_residue:
            def __init__(self, type):
                self.res_type1 = type
        class fake_protein:
            def __init__(self, sequence):
                self.residues = []
                for type in sequence:
                    self.residues.append(fake_residue(type))
        return fake_protein(sequence)
    
    def score_by_key(self, key1, key2, w):
        """ This is slow, and should only be used in testing. Use .score instead.
        """
        i,r = 0,0
        j,m = 0,0
        for res in self.target.residues:
            if res.res_number == key1:
                r = i
                break
            i += 1
        for res in self.template.residues:
            if res.res_number == key2:
                m = j
                break
            j += 1
        return self.score(r,m,w)
            
    def score(self, r, m):
        return 0.0

    def align(self):
        """ returns two lists of residues, of equal length. These residues are aligned.
        """
        print_param = self.print_param
        if len(self.template.residues) > len(self.target.residues):
            N = len(self.template.residues)+2
        else:
            N = len(self.target.residues)+2
            
        Max=xMax=yMax=0
        x_traceback = numpy.full(((N, N),-1), dytpe=numpy.int)         # initialize to -1
        y_traceback = numpy.full(((N, N),-1), dtype=numpy.int)
        score_table = numpy.zeros((N, N))
        aout = []
        bout = []
        for i in range(N*2):
            aout.append('-')
            bout.append('-')
        for i in range(1, len(self.target.residues)+1):
            lastgap = 'n'
            for j in range(1, len(self.template.residues)+1):
                # handle gap score
                down_gap_score = self.gap_score
                right_gap_score = self.gap_score
                if lastgap == 'y':
                    down_gap_score = self.gap_extension_score
                elif lastgap == 'x':
                    right_gap_score = self.gap_extension_score
                if self.alignment_type == 'semiglobal':
                    if i == len(self.target.residues):
                        down_gap_score = 0            # no gap score at last column and row for semiglobal
                    if j == len(self.template.residues):
                        right_gap_score = 0
                # calculate best direction
                diag   = score_table[i-1,j-1] + self.score(i-1,j-1)
                down   = score_table[i-1,j]   - down_gap_score
                right  = score_table[i,j-1]   - right_gap_score
                maxval = max(diag,down,right)
                if maxval <= 0 and self.alignment_type == 'local':
                    # set anything below 0 to 0
                    score_table[i,j] = 0.0
                elif maxval == diag:
                    lastgap = 'n'
                    score_table[i,j] = diag
                    x_traceback[i,j] = i-1
                    y_traceback[i,j] = j-1
                elif maxval == down:
                    lastgap = 'y'
                    score_table[i,j] = down
                    x_traceback[i,j] = i-1
                    y_traceback[i,j] = j
                else:
                    lastgap = 'x'
                    score_table[i,j] = right
                    x_traceback[i,j] = i
                    y_traceback[i,j] = j-1
                if self.alignment_type != 'global':
                    # search for the highest possible position
                    if maxval > Max:
                        Max=maxval
                        xMax=i
                        yMax=j
                #print score_table[i,j]
        if self.alignment_type == 'global':
            # want to start traceback from the bottom right corner
            i = len(self.target.residues)
            j = len(self.template.residues)
        else:
            # start traceback at highest starting point
            i=xMax; j=yMax
        target_index = i
        template_index = j
        x=y=1
        left_go=bottom_go=1
        while (target_index>0 and template_index>0):
            if self.alignment_type != 'global':
                # break if alignment score drops below zero
                if score_table[target_index,template_index] <= 0:
                    break
            if left_go and bottom_go:
                aout[x]=self.target.residues[target_index-1]
                bout[y]=self.template.residues[template_index-1]
                x += 1
                y += 1
            elif left_go:
                bout[y] = self.template.residues[template_index-1]
                y += 1
                x += 1
            elif bottom_go:
                aout[x] = self.target.residues[target_index-1]
                x += 1
                y += 1

            a, b = target_index-1   , x_traceback[target_index-1, template_index-1]
            c, d = template_index-1 , y_traceback[target_index-1, template_index-1]

            bottom_go = a > b
            left_go   = c > d
            tempi     = target_index-1
            tempj     = template_index-1
            target_index   = x_traceback[tempi,tempj]+1
            template_index = y_traceback[tempi,tempj]+1
                
        # delete any elements that have gaps in both positions
        delset = []
        for i in range(len(aout)):
            if aout[i] == '-' and bout[i] == '-':
                delset.append(i)
        delset.sort()
        delset.reverse()
        for delone in delset:  # last first so that the recorded indices remain valid
            del aout[delone]
            del bout[delone]

        a = ""
        b = ""
        for ind in range(len(aout)-1, -1, -1):
            if aout[ind] == '-':
                a += '-'
            else:
                a += aout[ind].res_type1
            if bout[ind] == '-':
                b += '-'
            else:
                b += bout[ind].res_type1
        
        if print_param:
            print a
            print b
            print '%s positions printed'%(len(a))

        aout.reverse()
        bout.reverse()
        return aout, bout
