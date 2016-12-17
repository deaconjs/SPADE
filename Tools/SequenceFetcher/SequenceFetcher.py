import os.path
import os
import sys
sys.path.append(os.path.normpath('../..'))
sys.path.append(os.path.normpath('..'))
sys.path.append(os.path.normpath('.'))
import MolecularSystem
sys.path.append('Tools/Aligner')
import SequenceAligner
import pickle
import string
import urllib
import copy
import re

try:
    from Bio.WWW import NCBI
except ImportError:
    pass
else:
    from Bio.Blast import NCBIWWW
    from Bio import Fasta
    from Bio.Clustalw import MultipleAlignCL
    import Bio.Clustalw

import global_functions

class SequenceFetcher:
    def __init__(self, system, pfam=1):
        """ pfam =1 gives pfam domain assignments
            pfam =0 builds a fasta-clustalw alignment file and returns the handle
        """
        self.system = system
        if pfam == 1:
            self.fetch_from = 'pfam'
            self.pfq = PFAMQuery()
        elif pfam == 0:
            self.fetch_from = 'fasta'

    def normalize_domain_labels(self):
        for pchain in self.system.ProteinList:
            try:
                pchain.residues[0].features['domain']
            except KeyError:
                continue
            max_domain = 0
            for res in pchain.residues:
                d = res.features['domain']
                if d > max_domain:
                    max_domain = d
            for res in pchain.residues:
                if res.features['domain'] == 0:
                    res.features['domain'] = -1
                else:
                    res.features['domain'] = res.features['domain'] / (max_domain+0.0)

    def fetch(self, forced_rewrite=0):
        if fetch_from == 'pfam':
            self.fetch_from_pfam(forced_rewrite)
        elif fetch_from == 'fasta':
            for pchain in self.system.ProteinList:
                sequence = pchain.get_sequence()
                cwaligned_tag = self.fetch_from_fasta(sequence, forced_rewrite)
                print 'tag = %s'%(cwaligned_tag)

    def fetch_from_fasta(self, sequence, forced_rewrite=0):            
        """ do the entire sequence """
        """
        b_tag = self.system.get_filename_by_extension('.fst')
        cwin_tag = self.system.get_filename_by_extension('.cin')
        dos_cwin_tag = global_functions.translate_filename_to_DOS_8_3_format(cwin_tag)
        cwout_tag = self.system.get_filename_by_extension('.clu')
        dos_cwout_tag = global_functions.translate_filename_to_DOS_8_3_format(cwout_tag)
        cwaligned_tag = self.system.get_filename_by_extension('.aln')
        dos_cwaligned_tag = global_functions.translate_filename_to_DOS_8_3_format(cwaligned_tag)
        #the cwprot_tag can only differ from the tag given to Clustalw_Protein by its extension
        print "BLAST\n"
        self.Blastp_( sequence, b_tag )
        print "BLAST Complete, Fetching Sequences\n"
        self.Create_Clustalw_Input( sequence, b_tag, dos_cwin_tag, 'ON')
        print "Fetch Complete, Clustalw Aligned Sequences\n"
        self.Clustalw_Align( dos_cwin_tag, dos_cwout_tag )
        print "Matching Fasta headers to alignments\n"
        self.Bridge_Fasta_Header_And_Aligned( dos_cwin_tag, dos_cwout_tag, dos_cwaligned_tag )
        return cwaligned_tag
        """
        pass
        

    def fetch_from_pfam(self, forced_rewrite=0):        
        bpm_filename = self.system.get_filename_by_extension('.bpm')   # bliss-pfam mapping
        # first get the alignments, whether from pfam or from file
        create_new = 0
        if forced_rewrite:
            create_new = 1
        else:
            try:
                bpm_file = open(bpm_filename, 'r')
            except IOError:
                create_new = 1
        if create_new:
            mapping = self.pfq.lookupMapping(self.system)
            if mapping == None:
                return
            bpm_file = open(bpm_filename, 'w')
            pickle.dump(mapping, bpm_file)
        else:
            mapping = pickle.load(bpm_file)
        # parse the mapping html to obtain the alignments
        names   = []
        accs    = []
        starts  = []
        ends    = []
        snames  = []
        subseqs = []
        chains  = []
        line_index = 0
        while 1:
            if line_index == len(mapping):
                break
            line = mapping[line_index]
            if string.find(line, 'name=name6') != -1:
                rex = re.compile('value=.*>')
                token = rex.findall(line)[0]
                names.append(token[6:-1])
            elif string.find(line, 'name=acc6') != -1:
                rex = re.compile('value=.*>')
                token = rex.findall(line)[0]
                accs.append(token[6:-1])
            elif string.find(line, 'name=start6') != -1:
                rex = re.compile('value=.*>')
                token = rex.findall(line)[0]
                starts.append(token[6:-1])
            elif string.find(line, 'name=end6') != -1:
                rex = re.compile('value=.*>')
                token = rex.findall(line)[0]
                ends.append(token[6:-1])
            elif string.find(line, 'name=sname6') != -1:
                rex = re.compile('value=.*>')
                token = rex.findall(line)[0]
                snames.append(token[6:-1])
            elif string.find(line, 'name=subseq6') != -1:
                rex = re.compile('value=.*>')
                token = rex.findall(line)[0]
                subseqs.append(token[6:-1])
            if string.find(line, 'matches PDB identifier') != -1:
                token = string.split(line,'font')[-2]
                token = string.split(token,'>',1)[1][:-2]
                chains.append(string.split(token, ":")[1])
            line_index += 1
        print accs
        print starts
        print ends

        # now do the same thing for the profiles: store as an array
        bpp_filename = self.system.get_filename_by_extension('.bpp')   # bliss-pfam profile
        # get the pfam profiles, whether from pfam or from file
        create_new = 0
        if forced_rewrite:
            create_new = 1
        else:
            try:
                bpp_file = open(bpp_filename, 'r')
            except IOError:
                create_new = 1
        if create_new:
            profiles = self.pfq.lookupProfiles(accs)
            if profiles == None:
                return
            sequences = []
            this_key = None
            for profile in profiles:
                sequences.append({})
                for seq in profile:
                    if seq[0] == '>':
                        sequences[-1][seq] = ['','']
                        this_key = seq
                    else:
                        gapped_sequence = string.replace(seq,'X','.')
                        sequences[-1][this_key][0] += gapped_sequence
                        # and delete any gaps. store the gapless one in the second string
                        ungapped_sequence = string.replace(gapped_sequence,'.','')
                        ungapped_sequence = string.replace(ungapped_sequence, '-','')
                        sequences[-1][this_key][1] += ungapped_sequence

            bpp_file = open(bpp_filename, 'w')
            profiles = sequences
            pickle.dump(profiles, bpp_file)
        else:
            profiles = pickle.load(bpp_file)
        print '%s profiles'%(len(profiles))

        # now label the domains
        for pchain in self.system.ProteinList:
            for res in pchain.residues:
                res.features['domain'] = 0
        for j in range(len(profiles)):
            aligner = Aligner.SequenceAligner()
            try:
                pchain = self.system.ProteinDict[chains[j]]
            except KeyError:
                continue
            aligner.clear_profile('template')
            aligner.clear_profile('target')
            aligner.add_sequence_to_profile(subseqs[j], 'template')
            # add the sequences to the aligner
            for key in profiles[j].keys():
                aligner.add_sequence_to_profile(profiles[j][key][0], 'target')
            # perform the alignment
            identity, aout, bout = aligner.align()
            print 'identity = %s'%(identity)
            # aout and bout contain indices to the two sequences. map domain labels back to residues
            for i in range(len(aout)):
                if aout[i] != -1:
                    if bout[i] != -1:
                        pchain = self.system.ProteinDict[chains[j]]
                        try:
                            pchain.residue_dict[bout[i]+string.atoi(starts[j])].features['domain'] = j+1
                        except:
                            print 'chain %s missing res %s'%(pchain.chain_name, bout[i]+string.atoi(starts[j]))

        last_label = -1
        last_resnum = -1
        for pchain in self.system.ProteinList:
            for res in pchain.residues:
                if res.features['domain'] != last_label:
                    if last_resnum == -1:
                        print '%s%s'%(res.res_type1, res.res_number),
                    else:
                        print '- %s%s %s\n%s%s'%(last_restype, last_resnum, last_label, res.res_type1, res.res_number),
                last_label = res.features['domain']
                last_resnum = res.res_number
                last_restype = res.res_type1
        print '- %s%s %s'%(last_restype, last_resnum, last_label)

#

class PFAMQuery:
    """PFAMQuery - Queries a PFAM server for a sequence alignment for a given PDB ID.    """
    def __init__ (self):
        """Initializes the PFAMQuery with the pfam_id, pdb_id mapping """
        # the web address for mapping pfam families to pdbs by pdb accession key
        self.mapping_addr='http://www.sanger.ac.uk/cgi-bin/Pfam/structural_view.pl?'
        # the web address for retrieving pfam profiles by pfam or uniprot key
        self.profile_addr='http://www.sanger.ac.uk/cgi-bin/Pfam/getalignment.pl?%s'

    def lookupMapping (self,system):
        bit = string.split(system.filename, os.path.normpath('/'))[-1]
        bit2 = string.split(bit, '.')[0]
        # first retrieve the pdb to pfam mapping and alignment
        parm_string = 'pdb='+string.lower(bit2[0:4])
        print 'connecting to %s'%(self.mapping_addr + parm_string)
        try:
            f = urllib.urlopen(self.mapping_addr + parm_string)
        except IOError:
            print 'please connect to the network first'
            sys.exit()
        return string.split(f.read(), '\n')
    
    def lookupProfiles(self, profile_keys):
        profiles = []
        last_profile_key = None
        for key in profile_keys:
            params=urllib.urlencode({'acc': key,'type': 'seed', 'format': 'fal'})
            try:
                print 'attempting request %s'%(self.profile_addr%params)
                f=urllib.urlopen(self.profile_addr % params)
            except IOError:
                print 'please connect to the network first'
                sys.exit()
            alignment = f.read()
            alignment=string.replace(alignment,"<pre>\n","")
            alignment=string.replace(alignment,"</pre>","")
            alignment=string.strip(alignment)
            profiles.append(string.split(alignment, '\n'))
        return profiles


if __name__ == '__main__':
    ms = MolecularSystem.System('None')
    database_dir = os.path.normpath('../Databases/OB/')
    subdirs = os.listdir(database_dir)
    #subdirs = ['1FK8_A']
    for subdir in subdirs:
        item = os.path.join(database_dir, subdir, subdir+'.pdb')
        ms.load_pdb(os.path.abspath(item))
        fetcher = SequenceFetcher(ms)
        fetcher.fetch()
        ms.save_system(os.path.join(database_dir, subdir, subdir+'.bps'))

