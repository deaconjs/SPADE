# python imports
import string
import math
import os.path
import pickle
import sys
import copy
# dependency imports
sys.path.append(os.path.abspath('./Dependencies'))
sys.path.append(os.path.abspath('./Dependencies/pybwidget-0.1.2_1.7.0'))
sys.path.append(os.path.abspath('./Dependencies/pybwidget/bwidget'))
import bwidget 
from Tkinter import *
from tkFileDialog import *
import Pmw


class Viewer:
    def __init__(self, parent, experiment, peak_index):
        # build the widget
        self.parent = parent
        self.rtop = Toplevel(self.parent)
        self.experiment = experiment
        self.my_filter_resolution = experiment.get_filter_resolution()
        self.rtop.title('Rationalize Peak')
        # top frame carries mutations options
        self.rtop.optionsframe = Frame(self.rtop, bd=8)
        self.rtop.mw_string = StringVar()
        self.rtop.mw_string.set('%s amu'%(self.parent.parent.parent.max_plot.x[peak_index]))
        self.rtop.mw_label = Label(self.rtop.optionsframe, bd=4, textvariable=self.rtop.mw_string)
        self.rtop.mutantsframe = Frame(self.rtop.optionsframe, relief=GROOVE, bd=2)
        self.rtop.mutantcheck_int = IntVar()
        self.rtop.mutantcheck = Checkbutton(self.rtop.mutantsframe,
                                            text="Check for Mutants",
                                            variable=self.rtop.mutantcheck_int)
        self.rtop.nonnatural_int = IntVar()
        self.rtop.include_nonnatural = Checkbutton(self.rtop.mutantsframe,
                                                   text="include non-natural",
                                                   variable = self.rtop.nonnatural_int,
                                                   padx=8)
        self.rtop.indel_int = IntVar()
        self.rtop.include_indels = Checkbutton(self.rtop.mutantsframe,
                                               text="include insertions/deletions",
                                               padx=12,
                                               variable = self.rtop.indel_int)
        self.rtop.mutant_range_counter = Pmw.Counter(self.rtop.mutantsframe,
                                                     labelpos = 'w',
                                                     label_text = 'maximum count:',
                                                     orient = 'horizontal',
                                                     entry_width = 2,
                                                     entryfield_value = 1,
                                                     datatype = {'counter':'numeric'},
                                                     entryfield_validate = {'validator':'numeric'},
                                                     padx=20,
                                                     increment=1)
        self.rtop.mutantcheck.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.include_nonnatural.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.include_indels.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.mutant_range_counter.pack(side=TOP, expand=NO, fill=NONE)
        self.rtop.mw_label.pack(anchor='w', side=TOP, expand=YES, fill=BOTH)
        self.rtop.mutantsframe.pack(side=TOP, expand=YES, fill=BOTH, pady=8)

        # Create and pack a RadioSelect widget, with radiobuttons.
        self.rtop.modifymutants_radio = Pmw.RadioSelect(self.rtop.optionsframe, buttontype = 'radiobutton', labelpos = 'w', command = None, label_text = 'Modify mutants?')
        for text in ('Yes', 'No'):
            self.rtop.modifymutants_radio.add(text)
        self.rtop.modifymutants_radio.invoke('No')
        self.rtop.modifymutants_radio.pack(side = TOP, expand = YES, fill=BOTH, pady=10)

        # bottom frame
        self.rtop.modsframe = Frame(self.rtop.optionsframe,  relief=GROOVE, bd=2)
        self.rtop.modscheck_int = IntVar()
        self.rtop.modscheck = Checkbutton(self.rtop.modsframe,
                                          text='Check for Modifications',
                                          variable = self.rtop.modscheck_int)
        self.rtop.modscheck.select()
        self.rtop.uncommon_int = IntVar()
        self.rtop.include_uncommon = Checkbutton(self.rtop.modsframe,
                                                 text='include uncommon',
                                                 padx=12,
                                                 variable = self.rtop.uncommon_int)
        self.rtop.nonspecific_int = IntVar()
        self.rtop.include_nonspecific = Checkbutton(self.rtop.modsframe,
                                                    text='include nonspecific',
                                                    padx=12,
                                                    variable = self.rtop.nonspecific_int)
        self.rtop.crosslinks_int = IntVar()
        self.rtop.include_crosslinks = Checkbutton(self.rtop.modsframe,
                                                   text='include crosslinks',
                                                   padx=12,
                                                   variable = self.rtop.crosslinks_int)
        self.rtop.mods_range_counter = Pmw.Counter(self.rtop.modsframe,
                                                   labelpos = 'w',
                                                   label_text='maximum count: ',
                                                   orient = 'horizontal',
                                                   entry_width = 2,
                                                   entryfield_value = 1,
                                                   datatype = {'counter':'numeric'},
                                                   entryfield_validate = {'validator':'numeric'},
                                                   padx=20,
                                                   increment=1)
        
        self.rtop.modscheck.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.include_uncommon.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.include_nonspecific.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.include_crosslinks.pack(anchor='w', side=TOP, expand=NO, fill=None)
        self.rtop.mods_range_counter.pack(anchor='w', side=TOP, expand=NO, fill=NONE)
        self.rtop.modsframe.pack(side=TOP, expand=YES, fill=BOTH, pady=8)

        self.rtop.limitsframe = Frame(self.rtop.optionsframe, relief=GROOVE, bd=2)
        
        self.rtop.fragweightlimit_counter = Pmw.Counter(self.rtop.limitsframe,
                                                        labelpos='w',
                                                        label_text='examine fragments within',
                                                        orient='horizontal',
                                                        entry_width=5,
                                                        entryfield_value=1000,
                                                        datatype={'counter':'numeric'},
                                                        entryfield_validate={'validator':'numeric'},
                                                        increment=50)
        self.rtop.ignorechangesbelow_counter = Pmw.Counter(self.rtop.limitsframe,
                                                           labelpos='w',
                                                           label_text='ignore weight changes below:',
                                                           orient='horizontal',
                                                           entry_width=3,
                                                           entryfield_value=3,
                                                           datatype={'counter':'numeric'},
                                                           entryfield_validate={'validator':'numeric'},
                                                           increment=1)
        self.rtop.my_filter_resolution_counter = Pmw.Counter(self.rtop.limitsframe,
                                                             labelpos='w',
                                                             label_text='apply filter resolution:',
                                                             orient='horizontal',
                                                             entry_width=6,
                                                             entryfield_value=self.my_filter_resolution,
                                                             datatype={'counter':'real'},
                                                             entryfield_validate={'validator':'real'},
                                                             increment=0.001)
        self.rtop.limitsframe.pack(side=TOP, expand=NO, fill=X, pady=10)
        self.rtop.fragweightlimit_counter.pack(side=TOP, expand=NO, fill=X, pady=10)
        self.rtop.ignorechangesbelow_counter.pack(side=TOP, expand=NO, fill=X,pady=10)
        self.rtop.my_filter_resolution_counter.pack(side=TOP, expand=NO, fill=X, pady=10)
        # button box
        self.rtop.rationalize_buttonbox = Pmw.ButtonBox(self.rtop.optionsframe, orient='horizontal')
        self.rtop.rationalize_buttonbox.add('Calculate', command=self._rerun_rationalize_peak)
        self.rtop.rationalize_buttonbox.add('Close', command=self.rtop.destroy)
        self.rtop.rationalize_buttonbox.pack(side=BOTTOM, expand=NO, fill=X)

        # now add a scrolled frame for holding fragment options and final output trees
        self.rtop.outframe = Frame(self.rtop, width=600, height=600)
        self.rtop.fragment_options_frame = Frame(self.rtop.outframe, width=600, height=200)
        self.rtop.output_frame = Frame(self.rtop.outframe, width=600, height=400)
        self.rtop.fragment_options_frame.pack(side=TOP, expand=NO, fill=X)
        self.rtop.output_frame.pack(side=TOP, expand=YES, fill=BOTH)
        self.rtop.optionsframe.pack(side=LEFT, expand=NO, fill=BOTH)
        self.rtop.outframe.pack(side=LEFT, expand=YES, fill=BOTH)

        self.starting_fragments = {}
        self.examination_fragments = {}

    def _rerun_rationalize_peak(self):
        self.starting_fragments = {}
        self.examination_fragments = {}
        self._run_rationalize_peak()

    def _run_rationalize_peak(self, calculation_fragments={}):
        aa_weights = {'A':71.09,  'C':103.15, 'D':115.09, 'E':129.12, 'F':147.18,
                      'G':57.05,  'H':137.14, 'I':113.16, 'K':128.17, 'L':113.16,
                      'M':131.19, 'N':114.11, 'P':97.12,  'Q':128.14, 'R':156.19,
                      'S':87.08,  'T':101.11, 'V':99.14,  'W':186.21, 'Y':163.18}
        self.my_filter_resolution = float(self.rtop.my_filter_resolution_counter.get())

        val = self.rtop.modifymutants_radio.getcurselection()
        if val == 'Yes':
            modify_mutants_val = 1
        else:
            modify_mutants_val = 0
        args = [float(string.split(self.rtop.mw_string.get())[0]),     # molecular weight
                int(self.rtop.mutantcheck_int.get()),                  # check primary sequence 
                int(self.rtop.nonnatural_int.get()),                   # allow nonnatural amino acids
                int(self.rtop.indel_int.get()),                        # allow insertions and deletions
                int(self.rtop.mutant_range_counter.get()),             # number of mutations tolerated
                int(self.rtop.modscheck_int.get()),                    # check post-translational modifications
                int(self.rtop.nonspecific_int.get()),                  # allow nonspecific reactions
                int(self.rtop.mods_range_counter.get()),               # number of post-translational modifications tolerated
                int(self.rtop.fragweightlimit_counter.get()),          # start with fragments within this limit of the selected peak
                int(self.rtop.ignorechangesbelow_counter.get()),       # any modifications or mutations that adjust the molecular weight by less than this are ignored
                modify_mutants_val,                                    # if selected ('YES'), apply the chemical modifications to mutants and indels
                int(self.rtop.crosslinks_int.get()),                   # allow crosslinks
                int(self.rtop.uncommon_int.get())]                     # include mods dict 

        frags = self.experiment.get_all_possible_fragment_objects()
        # make a hash of all possible fragment sequences, indexed by molecular weight
        if len(calculation_fragments.keys()) == 0:
            self.starting_fragments = {}
            sequences = self.experiment.get_protein_sequences()
            for frag in frags:
                # determine if it's the cterminal fragment for c-terminus-specific modifications
                is_cterm = 0
                if len(sequences[frag.get_chain()]) == frag.get_cterm_index():
                    is_cterm = 1
                if frag.get_weight() > args[0] - args[8] and frag.get_weight() < args[0] + args[8]:
                    print 'consider mw %s fragment %s'%(frag.get_weight(), frag.get_sequence())
                    try:
                        self.starting_fragments[int(round(frag.get_weight(), 0))]
                    except KeyError:
                        self.starting_fragments[int(round(frag.get_weight(), 0))] = [[frag.get_chain(), frag.get_sequence(), 0, 0, '', frag.get_nterm_index(), frag.get_cterm_index(), is_cterm, frag.get_weight()]]
                    else:
                        self.starting_fragments[int(round(frag.get_weight(), 0))].append([frag.get_chain(), frag.get_sequence(), 0, 0, '', frag.get_nterm_index(), frag.get_cterm_index(), is_cterm, frag.get_weight()])

        # here is the initial list. append new fragments to it where index 2 and 3 retain how many mutations/indels
        # and chemical adducts (respectively) have been added. This block generates all possibilities.
        if len(calculation_fragments.keys()) == 0:
            rationalization_fragments = copy.deepcopy(self.starting_fragments)
        else:
            rationalization_fragments = copy.deepcopy(calculation_fragments)
            
        if args[1]:
            print 'collecting mutations'
            if args[3]:
                print 'do deletions first so that mutations and insertions are not deleted'
                for i in range(args[4]):
                    count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'deletion', args[4], args[9])
            if args[2]:
                print 'add nonnatural mutations second because the X is ignored by the mutation producing code'
                for i in range(args[4]):
                    count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'nonnatural_mutation', args[4], args[9])
            if 1:
                print 'add normal mutations'
                for i in range(args[4]):
                    count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'mutation', args[4], args[9])
            if args[3]:
                print 'do insertions last so they are not mutated'
                for i in range(args[4]):
                    count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'insertion', args[4],args[9])
        mutant_holdout = {}
        if not args[10]:
            mutant_holdout = copy.deepcopy(rationalization_fragments)
            if len(calculation_fragments.keys()) == 0:
                rationalization_fragments = copy.deepcopy(self.starting_fragments)
            else:
                rationalization_fragments = copy.deepcopy(calculation_fragments)

        if args[5]:
            print 'collecting post-translational modifications'
            for i in range(args[7]):
                count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'post-translational modification', args[7],args[9], args[12])
            if args[6]:
                print 'and nonspecific post-translational modifications'
                for i in range(args[7]):
                    count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'nonspecific reactions', args[7],args[9])
            if args[11]:
                print 'and crosslinks'
                for i in range(args[7]):
                    count, rationalization_fragments = self._collect_rationalization_fragset(rationalization_fragments, 'crosslinks', args[7],args[9])

        if len(mutant_holdout.keys()) > 0:              # if mutants are not being modified
            # merge mutant and modification dictionaries
            for mkey in mutant_holdout.keys():
                for rkey in rationalization_fragments.keys():
                    if mkey == rkey:
                        for mutant_fragment in mutant_holdout[mkey]:
                            rationalization_fragments[rkey].append(mutant_fragment)
                        break
                else:                            
                    rationalization_fragments[mkey] = mutant_holdout[mkey]

        # now that all possiblities have been enumerated, filter out those with an appropriate molecular weight
        filtered_mutants = self._filter_rationalization_fragset(rationalization_fragments)
        total_count = 0
        for key in filtered_mutants.keys():
            total_count += len(filtered_mutants[key])
        print '%s after filtering with %s molecular weights'%(total_count, len(filtered_mutants))
        self.update_rationalization_output(filtered_mutants)

    def _collect_rationalization_fragset(self, input_fragset, type, how_many, ignore_changes_below, include_uncommon=0):
        aa_weights = {'A':71.09,  'C':103.15, 'D':115.09, 'E':129.12, 'F':147.18,
                      'G':57.05,  'H':137.14, 'I':113.16, 'K':128.17, 'L':113.16,
                      'M':131.19, 'N':114.11, 'P':97.12,  'Q':128.14, 'R':156.19,
                      'S':87.08,  'T':101.11, 'V':99.14,  'W':186.21, 'Y':163.18}

        # first load the common modifications and build the mods dictionary appropriately
        file = open('./Tools/MassSpectrometry/reactions_dict.pkl')
        common_mods = pickle.load(file)
        file.close()
        mods = {'A':{}, 'C':{}, 'D':{}, 'E':{}, 'F':{}, 'G':{}, 'H':{},
                'I':{}, 'K':{}, 'L':{}, 'M':{}, 'N':{}, 'P':{}, 'Q':{},
                'R':{}, 'S':{}, 'T':{}, 'V':{}, 'W':{}, 'Y':{}}
        for key in common_mods.keys():
            for target_AA in common_mods[key]['target_AA']:
                if target_AA not in mods.keys():
                    mods[target_AA] = {key:float(common_mods[key]['added_weight'])}
                else:
                    mods[target_AA][key] = float(common_mods[key]['added_weight'])
                
        file = open('./Tools/MassSpectrometry/mods.pkl', 'rb')
        uncommon_mods = pickle.load(file)
        file.close()
        for token in ['nn', 'crosslinks', 'X', 'cterm', 'nterm', 'nsp']:
            mods[token] = copy.deepcopy(uncommon_mods[token])
            
        if include_uncommon:
            # open the modifications database
            for key in uncommon_mods.keys():
                if key not in mods.keys():
                    mods[key] = copy.deepcopy(uncommon_mods[key])
                for reaction_key in uncommon_mods[key].keys():
                    if reaction_key not in mods[key].keys():
                        mods[key][reaction_key] = copy.deepcopy(uncommon_mods[key][reaction_key])

        return_set = {}
        # have the original fragment
        if type == 'deletion':
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[2] < how_many:
                        for i in range(len(sf[1])):
                            if aa_weights[sf[1][i]] > ignore_changes_below:
                                try:
                                    return_set[int(round(sf_key-aa_weights[sf[1][i]]))]
                                except KeyError:
                                    return_set[int(round(sf_key-aa_weights[sf[1][i]]))] = [[sf[0], sf[1][:i] + sf[1][i+1:], sf[2]+1, sf[3], sf[4]+'%s deletion at position %s (-%s)'%(sf[1][i], i, aa_weights[sf[1][i]]),sf[5],sf[6],sf[7], sf[8]]]
                                else:
                                    return_set[int(round(sf_key-aa_weights[sf[1][i]]))].append([sf[0], sf[1][:i] + sf[1][i+1:], sf[2]+1, sf[3], sf[4]+'%s deletion at position %s (-%s)'%(sf[1][i], i, aa_weights[sf[1][i]]),sf[5],sf[6],sf[7], sf[8]])

        elif type == 'insertion':
            # have the original fragment and all possible single or multiple deletions
            #
            # put the insertion in the second position, so as to not obstruct terminus analysis
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[2] < how_many:
                        for aa in aa_weights.keys():
                            if aa_weights[aa] > ignore_changes_below:
                                try:
                                    return_set[int(round(sf_key+aa_weights[aa]))]
                                except KeyError:
                                    return_set[int(round(sf_key+aa_weights[aa]))] = [[sf[0], sf[1][0] + aa + sf[1][1:], sf[2]+1, sf[3], sf[4]+'%s insertion (+%s)'%(aa, aa_weights[aa]),sf[5],sf[6],sf[7], sf[8]]]
                                else:
                                    return_set[int(round(sf_key+aa_weights[aa]))].append([sf[0], sf[1][0] + aa + sf[1][1:], sf[2]+1, sf[3], sf[4] + '%s insertion (+%s)'%(aa, aa_weights[aa]),sf[5],sf[6],sf[7], sf[8]])

        elif type == 'mutation':
            # have the original fragment and all possible single or multiple insertions and deletions
            # generate all possible mutants
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    sequence = sf[1]
                    if sf[2] < how_many:
                        for i in range(len(sequence)):              # parse this fragment's sequence
                            for aa in aa_weights.keys():
                                if abs(aa_weights[aa]-aa_weights[sequence[i]]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key-aa_weights[sequence[i]]+aa_weights[aa]))]
                                    except KeyError:
                                        return_set[int(round(sf_key-aa_weights[sequence[i]]+aa_weights[aa]))] = [[sf[0], sequence[:i] + aa + sequence[i+1:], sf[2]+1, sf[3], sf[4] + '%s to %s mutation at position %s (%s)'%(sequence[i], aa, i, aa_weights[aa]-aa_weights[sequence[i]]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key-aa_weights[sequence[i]]+aa_weights[aa]))].append([sf[0], sequence[:i] + aa + sequence[i+1:], sf[2]+1, sf[3], sf[4] + '%s to %s mutation at position %s (%s)'%(sequence[i], aa, i, aa_weights[aa]-aa_weights[sequence[i]]),sf[5],sf[6],sf[7], sf[8]])

        elif type == 'nonnatural_mutation':                                
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[2] < how_many:
                        for i in range(len(sf[1])):
                            for nn in mods['nn'].keys():
                                if abs(mods['nn'][nn]-aa_weights[sf[1][i]]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key-aa_weights[sf[1][i]]+mods['nn'][nn]))]
                                    except KeyError:
                                        return_set[int(round(sf_key-aa_weights[sf[1][i]]+mods['nn'][nn]))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2]+1, sf[3], sf[4] + '%s to %s mutation at position %s (%s)'%(sf[1][i], nn, i, mods['nn'][nn]-aa_weights[sf[1][i]]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key-aa_weights[sf[1][i]]+mods['nn'][nn]))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2]+1, sf[3], sf[4] + '%s to %s mutation at position %s (%s)'%(sf[1][i], nn, i, mods['nn'][nn]-aa_weights[sf[1][i]]),sf[5],sf[6],sf[7], sf[8]])

        elif type == 'post-translational modification':
            j = 0
            k = 0
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[3] < how_many:
                        for i in range(len(sf[1])):        # for each aa in the sequence
                            modset = mods[sf[1][i]]
                            for modkey in modset.keys():
                                if abs(modset[modkey]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key+modset[modkey], 0))]
                                    except KeyError:
                                        return_set[int(round(sf_key+modset[modkey], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2], sf[3]+1, sf[4] + '%s at position %s (%s)'%(modkey, i, modset[modkey]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key+modset[modkey], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2], sf[3]+1, sf[4] + '%s at position %s (%s)'%(modkey, i, modset[modkey]),sf[5],sf[6],sf[7], sf[8]])
                                    j += 1
                                    if j%1000000 == 999999:
                                        k += 1
                                        print '%sM possibilities'%(k)
                        if sf[5] == 0:    # if is the n-term
                            modset = mods['nterm']
                            for modkey in modset.keys():
                                if abs(modset[modkey]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key+modset[modkey], 0))]
                                    except KeyError:
                                        return_set[int(round(sf_key+modset[modkey], 0))] = [[sf[0], 'X' + sf[1][1:], sf[2], sf[3]+1, sf[4] + '%s at n-terminus 0 (%s)'%(modkey, modset[modkey]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key+modset[modkey], 0))].append([sf[0], 'X' + sf[1][1:], sf[2], sf[3]+1, sf[4] + '%s at n-terminus 0 (%s)'%(modkey, modset[modkey]),sf[5],sf[6],sf[7], sf[8]])
                        if sf[7] == 1:    # if includes the n-terms
                            modset = mods['cterm']
                            for modkey in modset.keys():
                                if abs(modset[modkey]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key+modset[modkey], 0))]
                                    except KeyError:
                                        return_set[int(round(sf_key+modset[modkey], 0))] = [[sf[0], sf[1][:-1] + 'X', sf[2], sf[3]+1, sf[4] + '%s at n-terminus 0 (%s)'%(modkey, modset[modkey]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key+modset[modkey], 0))].append([sf[0], sf[1][:-1] + 'X', sf[2], sf[3]+1, sf[4] + '%s at n-terminus (%s)'%(modkey, modset[modkey]),sf[5],sf[6],sf[7], sf[8]])

        elif type == 'nonspecific reactions':
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[3] < how_many:
                        for i in range(len(sf[1])):        # for each aa in the sequence
                            modset = mods['nsp']
                            for modkey in modset.keys():
                                if abs(mods['nsp'][modkey]-aa_weights[sf[1][i]]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key+mods['nsp'][modkey], 0))]
                                    except KeyError:
                                        return_set[int(round(sf_key+mods['nsp'][modkey], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2], sf[3]+1, sf[4] + 'nonspecific %s at position %s (%s)'%(modkey, i, mods['nsp'][modkey]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key+mods['nsp'][modkey], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2], sf[3]+1, sf[4] + 'nonspecific %s at position %s (%s)'%(modkey, i, mods['nsp'][modkey]),sf[5],sf[6],sf[7], sf[8]])
        elif type == 'X':
            # this block, and the whole 'X' option here is 
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[3] < how_many:
                        for i in range(len(sf[1])):        # for each aa in the sequence
                            modset = mods['X']
                            for modkey in modset.keys():
                                if abs(mods['X'][modkey]-aa_weights[sf[1][i]]) > ignore_changes_below:
                                    try:
                                        return_set[int(round(sf_key+mods['X'][modkey], 0))]
                                    except KeyError:
                                        return_set[int(round(sf_key+mods['X'][modkey], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2], sf[3]+1, sf[4] + 'nonspecific %s at position %s (%s)'%(modkey, i, mods['X'][modkey]),sf[5],sf[6],sf[7], sf[8]]]
                                    else:
                                        return_set[int(round(sf_key+mods['X'][modkey], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:], sf[2], sf[3]+1, sf[4] + 'nonspecific %s at position %s (%s)'%(modkey, i, mods['X'][modkey]),sf[5],sf[6],sf[7], sf[8]])
        elif type == 'crosslinks':
            # look for every combination of the crosslinking partners
            # conceptually, the code works like post translational modifications
            for sf_key in input_fragset.keys():
                for sf in input_fragset[sf_key]:
                    if sf[3] < how_many:
                        k = len(sf[1])
                        for i in range(len(sf[1])-1):        # for each aa in the sequence
                            for j in range(i+1,len(sf[1])):
                                coupling1 = '%s%s'%(sf[1][i], sf[1][j])
                                coupling2 = '%s%s'%(sf[1][j], sf[1][i])
                                coupling3 = 'XX'
                                coupling4 = 'XX'
                                if i==0:
                                    coupling3 == '0%s'%(sf[1][j])
                                    coupling4 == '%s0'%(sf[1][j])
                                elif j==0:
                                    coupling3 == 'O%s'%(sf[1][j])
                                    coupling4 == '%sO'%(sf[1][j])
                                if coupling1 in mods['crosslinks'].keys():
                                    if mods['crosslinks'][coupling1] > ignore_changes_below:
                                        modsdic = mods['crosslinks'][coupling1]
                                        for rxn in modsdic.keys():
                                            try:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))]
                                            except KeyError:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]]]
                                            else:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + ' crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]])
                                elif coupling2 in mods['crosslinks'].keys():
                                    if mods['crosslinks'][coupling2] > ignore_changes_below:
                                        modsdic = mods['crosslinks'][coupling2]
                                        for rxn in modsdic.keys():
                                            try:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))]
                                            except KeyError:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]]]
                                            else:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]])
                                elif coupling3 in mods['crosslinks'].keys():
                                    if mods['crosslinks'][coupling3] > ignore_changes_below:
                                        modsdic = mods['crosslinks'][coupling3]
                                        for rxn in modsdic.keys():
                                            try:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))]
                                            except KeyError:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]]]
                                            else:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]])
                                elif coupling4 in mods['crosslinks'].keys():
                                    if mods['crosslinks'][coupling4] > ignore_changes_below:
                                        modsdic = mods['crosslinks'][coupling4]
                                        for rxn in modsdic.keys():
                                            try:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))]
                                            except KeyError:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))] = [[sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]]]
                                            else:
                                                return_set[int(round(sf_key+modsdic[rxn], 0))].append([sf[0], sf[1][:i] + 'X' + sf[1][i+1:j] + 'X' + sf[1][j+1:], sf[2], sf[3]+1, sf[4] + 'crosslink %s between positions %s and %s (%s)'%(rxn, i, j, modsdic[rxn]),sf[5],sf[6],sf[7], sf[8]])
        # enter the return set into the input set.
        redundancy_count = 0
        redundancy_filter_on = 0
        for key in return_set.keys():
            try:
                input_fragset[key]
            except KeyError:
                input_fragset[key] = return_set[key]  # if its not there, append the whole list
            else:
                if redundancy_filter_on:
                    for return_fragment in return_set[key]:  # if it is there, append the elements of the return list to the existing input list
                        for input_fragment in input_fragset[key]:
                            if input_fragment[5] == return_fragment[5] and input_fragment[6] == return_fragment[6] and input_fragment[0] == return_fragment[0] and input_fragment[4] == return_fragment[4]:
                                redundancy_count += 1
                                break
                        else:
                            input_fragset[key].append(return_fragment)
                else:
                    input_fragset[key] += return_set[key]

        total_count = 0
        for key in input_fragset.keys():
            total_count += len(input_fragset[key])

        print '%s possibilities cumulative'%(total_count)            

        return total_count, input_fragset

    def _filter_rationalization_fragset(self, input_set):
        skeys = input_set.keys()
        skeys.sort()
        actual = float(string.split(self.rtop.mw_string.get())[0])

        # filter and remove redundancies.
        return_set = {}
        redundancy_count = 0
        sizefiltercount = 0
        for key in skeys:
            if abs(key-actual) < self.my_filter_resolution*actual:
                return_set[key] = []
                for frag1 in input_set[key]:
                    if len(return_set[key]) == 0:
                        return_set[key].append(frag1)
                    else:
                        for frag2 in return_set[key]:
                            if frag1[5]==frag2[5] and frag1[6]==frag2[6] and frag1[0]==frag2[0] and frag1[4]==frag2[4]:
                                redundancy_count += 1
                                break
                        else:
                            return_set[key].append(frag1)
            else:
                sizefiltercount += len(input_set[key])

        total_count = 0
        for key in return_set.keys():
            total_count += len(return_set[key])

        skeys = return_set.keys()
        skeys.sort()
        if len(skeys) == 0:
            print 'no possible modified fragments match the peak molecular weight'
            return {}
        print '%s filtered by molecular weight spanning %s to %s'%(sizefiltercount, skeys[0], skeys[-1])
        print '%s remaining in %s molecular weight bins'%(total_count, len(return_set.keys()))
        print '%s redundancies'%(redundancy_count)


        aa_weights = {'A':71.09,  'C':103.15, 'D':115.09, 'E':129.12, 'F':147.18,
                      'G':57.05,  'H':137.14, 'I':113.16, 'K':128.17, 'L':113.16,
                      'M':131.19, 'N':114.11, 'P':97.12,  'Q':128.14, 'R':156.19,
                      'S':87.08,  'T':101.11, 'V':99.14,  'W':186.21, 'Y':163.18, 'X':0}


        fout = open('./rationalization_out.txt', 'w')
        for key in skeys:
            fout.write('\n\n     %s daltons (+/- %s)\n'%(key, key * self.my_filter_resolution))
            for frag in return_set[key]:
                fragwt = 18.0
                for aa in frag[1]:
                    fragwt += aa_weights[aa]
                fout.write('%3s %s %s %s %s %s %s\n'%(frag[5], frag[6], frag[2], frag[3], frag[8], frag[1], frag[4]))
        fout.close()
        return return_set

    def _recalculate_with_fragment_filter(self):
        self.examination_fragments = {}
        cursel = self.rtop.fragselection_radio.getcurselection()
        # capture the selected fragments as input to run_rationalize_peak
        calculation_fragments = {}
        for sel in cursel:
            tokens = string.split(sel)
            chain1, nterm1, cterm1 = tokens[1], tokens[3], tokens[5]
            print chain1, nterm1, cterm1
            for key in self.starting_fragments.keys():
                for st_frag in self.starting_fragments[key]:
                    chain2, nterm2, cterm2 = st_frag[0], st_frag[5], st_frag[6]
                    if int(nterm1) == int(nterm2) and int(cterm1) == int(cterm2) and int(chain1) == int(chain2):
                        wt = key
                        try:
                            calculation_fragments[wt]
                        except KeyError:
                            calculation_fragments[wt] = [st_frag]
                        else:
                            calculation_fragments[wt].append(st_frag)
                        break
        self._run_rationalize_peak(calculation_fragments)

    def toggle_fragments_for_filter(self):
        list = []
        for i in range(self.rtop.fragselection_radio.numbuttons()):
            self.rtop.fragselection_radio.invoke(i)
        
    def update_rationalization_output(self, input_modified_fragments):
        # first copy over the currently selected indices for regeneration
        if len(self.examination_fragments.keys()) == 0:
            self.examination_fragments = input_modified_fragments

        # first make a list of indices that are covered by examination fragments
        # make a table of all fragments that are affected by the indices covered by these
        sequences = self.experiment.get_protein_sequences()
        indices = []
        i = 0
        for seq in sequences:
            indices.append([])
            for aa in seq:
                indices[-1].append([])
            for exkey in self.examination_fragments.keys():
                for exfrag in self.examination_fragments[exkey]:
                    if exfrag[0] == i:
                        for j in range(exfrag[5], exfrag[6]):
                            indices[i][j].append(exfrag)
            i += 1

        stored_fragment_selection_indices = []
        try:
            self.rtop.fragselection_radio
        except AttributeError:
            pass
        else:
            cur_selection = self.rtop.fragselection_radio.getcurselection()
            numbuttons    = self.rtop.fragselection_radio.numbuttons()
            for i in range(numbuttons):
                stored_fragment_selection_indices.append(0)
            for i in range(len(cur_selection)):
                stored_fragment_selection_indices[self.rtop.fragselection_radio.index(cur_selection[i])] = 1
        
        # first destroy any old copy of the output frame contents
        self.rtop.fragment_options_frame.destroy()
        self.rtop.output_frame.destroy()
        # and recreate frames
        self.rtop.fragment_options_frame = Frame(self.rtop.outframe, width=600, height=200)
        self.rtop.output_frame = Frame(self.rtop.outframe, width=600, height=400)
        # and subframes
        self.rtop.fragment_buttons_frame = Frame(self.rtop.fragment_options_frame)
        self.rtop.fragment_radio_frame = Pmw.ScrolledFrame(self.rtop.fragment_options_frame, usehullsize=1, hull_width=600, hull_height=200)
        # and widgets
        self.rtop.apply_fragment_filter_button  = Button(self.rtop.fragment_buttons_frame, text='Recalculate', command=self._recalculate_with_fragment_filter)
        
        # first make fragments available for selection/deselection and filtering
        tokens = {}
        tokens_counter = {}
        for key in self.examination_fragments.keys():
            for frag in self.examination_fragments[key]:
                if frag[0] >= 0 and frag[0] < 10:
                    part1 = '0%s'%(frag[0])
                elif frag[0] >= 10:
                    part1 = '%s'%(frag[0])
                if frag[5] >= 0 and frag[5] < 10:
                    part2 = '000%s'%(frag[5])
                elif frag[5] >= 10 and frag[5] < 100:
                    part2 = '00%s'%(frag[5])
                elif frag[5] >= 100 and frag[5] < 1000:
                    part2 = '0%s'%(frag[5])
                elif frag[5] >= 1000:
                    part2 = '%s'%(frag[5])
                if frag[6] >= 0 and frag[6] < 10:
                    part3 = '000%s'%(frag[6])
                elif frag[6] >= 10 and frag[6] < 100:
                    part3 = '00%s'%(frag[6])
                elif frag[6] >= 100 and frag[6] < 1000:
                    part3 = '0%s'%(frag[6])
                elif frag[6] >= 1000:
                    part3 = '%s'%(frag[6])
                    
                index = '%s %s %s'%(part1, part2, part3)
                if index not in tokens.keys():
                    tokens[index] = [frag[0], frag[1], frag[5], frag[6], frag[7], frag[8]]
                    tokens_counter[index] = 1
                else:
                    tokens_counter[index] += 1

        stokens = tokens.keys()
        stokens.sort()

        tokenlist = []
        for token in stokens:
            stuff = string.split(token, " ")
            sequence = tokens[token][1]
            tokenlist.append('chain %s positions %s to %s %6.2famu %s possibilites - %s'%(int(stuff[0]), int(stuff[1]), int(stuff[2]), tokens[token][5], tokens_counter[token], sequence))
            
        # and create the radiobuttons
        self.rtop.fragselection_radio = Pmw.RadioSelect(self.rtop.fragment_radio_frame.interior(),
                                                        buttontype   = 'checkbutton',
                                                        orient       = 'vertical', 
                                                        labelpos     = 'n',
                                                        pady         = 1,
                                                        command      = None,
                                                        label_text   = 'Deselect uninteresting fragments:',
                                                        selectmode   = 'multiple')
        
        if len(stored_fragment_selection_indices) != len(tokenlist):
            for text in tokenlist:
                self.rtop.fragselection_radio.add(text)
                self.rtop.fragselection_radio.invoke(text)
        else:
            j = 0
            for text in tokenlist:
                self.rtop.fragselection_radio.add(text)
                if stored_fragment_selection_indices[j]:
                    self.rtop.fragselection_radio.invoke(text)
                j += 1

        self.rtop.toggle_fragments_button   = Button(self.rtop.fragment_buttons_frame, text='Toggle selection', command=self.toggle_fragments_for_filter)

        self.rtop.fragselection_radio.pack(side = TOP, expand = YES, fill=BOTH, pady=10)
        self.rtop.apply_fragment_filter_button.pack(side=LEFT, anchor='nw', expand=NO, fill=NONE)
        self.rtop.toggle_fragments_button.pack(side=LEFT, anchor='nw', expand=NO, fill=NONE)
        self.rtop.fragment_buttons_frame.pack(side=TOP, expand=NO, fill=X)
        self.rtop.fragment_radio_frame.pack(side=TOP, expand=YES, fill=X)
        self.rtop.fragment_options_frame.pack(side=TOP, expand=NO, fill=X)

        self.rtop.accuracy_notebook = Pmw.NoteBook(self.rtop.output_frame)
        self.rtop.accuracy_notebook.pack(fill=BOTH, expand=1)

        optimal_weight = float(string.split(self.rtop.mw_string.get())[0])

        # find the range of accuracies
        closest = 1000.0
        furthest = 0.0
        for xfkey in self.examination_fragments.keys():
            if abs(optimal_weight-xfkey) > furthest:
                furthest = int(math.ceil(abs(optimal_weight-xfkey)))
            if abs(optimal_weight-xfkey) < closest:
                closest = int(math.floor(abs(optimal_weight-xfkey)))
        increment = int(math.ceil((furthest-closest)/4.0))

        if increment > 3.0:
            bottomrange = [0.0, closest+increment]
            secondrange = [bottomrange[1],  bottomrange[1]+increment]
            thirdrange =  [secondrange[1], secondrange[1]+increment]
            fourthrange = [thirdrange[1], furthest+1]
            range_keys = ['< %s amu'%(bottomrange[1]),
                          '%s-%s amu'%(secondrange[0],secondrange[1]),
                          '%s-%s amu'%(thirdrange[0], thirdrange[1]),
                          '%s-%s amu'%(fourthrange[0],fourthrange[1])]
            ranges = {range_keys[0]:bottomrange,
                      range_keys[1]:secondrange,
                      range_keys[2]:thirdrange,
                      range_keys[3]:fourthrange}
        elif increment > 2.0:
            increment = int(math.ceil((furthest-closest)/3.0))
            bottomrange = [0.0, closest+increment]
            secondrange = [bottomrange[1],  bottomrange[1]+increment]
            thirdrange =  [secondrange[1], secondrange[1]+increment]
            range_keys = ['< %s amu'%(bottomrange[1]),
                          '%s-%s amu'%(secondrange[0],secondrange[1]),
                          '%s-%s amu'%(thirdrange[0], thirdrange[1])]
            ranges = {range_keys[0]:bottomrange,
                      range_keys[1]:secondrange,
                      range_keys[2]:thirdrange}
        elif increment > 1.0:
            increment = int(math.ceil((furthest-closest)/2.0))
            bottomrange = [0.0, closest+increment]
            secondrange = [bottomrange[1],  bottomrange[1]+increment]
            range_keys = ['< %s amu'%(bottomrange[1]),
                          '%s-%s amu'%(secondrange[0],secondrange[1])]
            ranges = {range_keys[0]:bottomrange,
                      range_keys[1]:secondrange}
        else:
            increment = int(math.ceil((furthest)))
            bottomrange = [0.0, closest+increment]
            range_keys = ['< %s amu'%(bottomrange[1])]
            ranges = {range_keys[0]:bottomrange}

        # use 
        trees  = {}
        apply_buttons = {}
        i = -1
        for range_key in range_keys:
            i += 1
            page = self.rtop.accuracy_notebook.add(range_key)
            self.rtop.accuracy_notebook.tab(range_key).focus_set()
            trees[range_key] = [None, None]
            apply_buttons[range_key] = [None, None]
            treepanes = Pmw.PanedWidget(page,
                                        hull_borderwidth=1,
                                        orient='horizontal',
                                        hull_relief='sunken',
                                        hull_width=250)
            mutation_pane = treepanes.add('mutation pane', min=.2, max=.8, size=0.5)
            modification_pane = treepanes.add('modification pane', min=.2, max=.8, size=0.5)
            # mutation tree
            mut_scframe = Pmw.ScrolledFrame(mutation_pane, horizflex='elastic', vertflex='elastic', vscrollmode='static', usehullsize=1, hull_height=800)
            mod_scframe = Pmw.ScrolledFrame(modification_pane, horizflex='elastic', vertflex='elastic', vscrollmode='static', usehullsize=1, hull_height=800)
            mut_button_box = Pmw.ButtonBox(mutation_pane)
            mod_button_box = Pmw.ButtonBox(modification_pane)
            apply_buttons[range_key][0] = mut_button_box.add('Apply')
            apply_buttons[range_key][0].config(state='disabled')
            apply_buttons[range_key][1] = mod_button_box.add('Apply')
            apply_buttons[range_key][1].config(state='disabled')
            
            def _update_mut_apply_button(widget, args, range_key):
                text = args
                tokens = string.split(args, '} {')
                if len(tokens) == 1 and len(args) > 0:
                    apply_buttons[range_key][0].config(state='active')
                else:
                    apply_buttons[range_key][0].config(state='disabled')
                    
            def _update_mod_apply_button(widget, args, range_key):
                text = args
                tokens = string.split(args, '} {')
                if len(tokens) == 1 and len(args) > 0:
                    apply_buttons[range_key][1].config(state='active')
                else:
                    apply_buttons[range_key][1].config(state='disabled')
                    
            c_lambda = lambda widget, args, x=range_key:_update_mut_apply_button(widget, args, x)
            trees[range_key][0] = bwidget.Tree(mut_scframe.interior(), dragenabled=0, height=800, selectcommand=c_lambda)
            c_lambda = lambda widget, args, x=range_key:_update_mod_apply_button(widget, args, x)
            trees[range_key][1] = bwidget.Tree(mod_scframe.interior(), dragenabled=0, height=800, selectcommand=c_lambda)

            def _apply_mutation(key):
                print trees[key][0]

            def _apply_update_weight_tags(key):
                weight_tags = {}
                # find the node
                x = trees[key][1].selection_get()
                tokens = string.split(x[0])
                # find the sequence
                sequences = self.experiment.get_protein_sequences()
                sequence = sequences[int(tokens[0])]
                # the modification
                modification = ""
                modification_tokens = tokens[3:-8]
                for mod_token in modification_tokens:
                    modification += mod_token
                # and the weight
                weight = float(tokens[-5][1:-1])
                # and send it back to the experiment
                print 'weight tags %s'%(int(tokens[1])+int(tokens[-6]))
                print weight, modification, sequence
                print self.parent.parent.parent.parent.parent
                self.parent.parent.parent.parent.parent.update_modification_weights({sequence:[[int(tokens[1])+int(tokens[-6]), weight]]})
                # and rerun the analysis, closing this rationalization dialog

                # recalculate fragments
                self.experiment._calculate_single_reaction_fragments_dictionary()
                self.experiment._calculate_unreactive_fragments_dictionary()
                self.experiment._calculate_all_possible_fragments_dictionary()
        
                # recognize peaks
                self.parent.parent.parent.set_peak_rec_possible_peaks(self.experiment.get_single_reaction_fragment_weights(0))
                self.parent.parent.parent.set_all_fragment_weights(self.experiment.get_all_possible_fragments(0))
                self.parent.parent.parent.calculate_peaks()
                self.parent.parent.parent.draw_PlotPanels()
                self.parent.parent.parent.draw_peaks()
                

            c_lambda = lambda a=range_key: _apply_mutation(a)
            apply_buttons[range_key][0].config(command=c_lambda)
            c_lambda = lambda a=range_key: _apply_update_weight_tags(a)
            apply_buttons[range_key][1].config(command=c_lambda)
            # collect altered fragments in this accuracy range
            # make a dictionary of fragments, indexed by comment, for this accuracy range
            muts_dict = {}
            mods_dict = {}
            for weight in self.examination_fragments.keys():
                # chain, sequence, # mutations, # modifications, comment, nterm, cterm, iscterm, original weight
                if abs(weight-optimal_weight) >= ranges[range_key][0] and abs(weight-optimal_weight) < ranges[range_key][1]:
                    for altered_frag in self.examination_fragments[weight]:
                        if altered_frag[2] == 1 and altered_frag[3] == 0:
                            tokens = string.split(altered_frag[4])
                            if tokens[1] == 'deletion' or tokens[1] == 'insertion':
                                key = '%s %s'%(tokens[0], tokens[1])
                            else:
                                key = '%s %s %s'%(tokens[0], tokens[1], tokens[2])
                            try:
                                muts_dict[key]
                            except KeyError:
                                muts_dict[key] = [altered_frag]
                            else:
                                muts_dict[key].append(altered_frag)
                        elif altered_frag[2] == 0 and altered_frag[3] == 1:
                            tokens = string.split(altered_frag[4])
                            key = ""
                            for j in range(0,len(tokens)-4):
                                key += '%s '%(tokens[j])
                            try:
                                mods_dict[key]
                            except KeyError:
                                mods_dict[key] = [altered_frag]
                            else:
                                mods_dict[key].append(altered_frag)
            print 'range %s - %s in muts, %s in mods'%(range_key, len(muts_dict.keys()), len(mods_dict.keys()))
            tree_store1 = {}
            tree_store2 = {}
            tree_store3 = {}
            
            # query the host experimental PlotWindow for any peaks that are not rationalized
            unrationalized_peak_indices = self.parent.parent.parent.get_unrationalized_peaks()
            rationalized_peaks = self.parent.parent.parent.get_reactions()

            for muts_key in muts_dict.keys():
                tree_store1[muts_key] = [0, 'root', '(%s)  %s'%(len(muts_dict[muts_key]), muts_key)]
                i = 1
                for mutated_fragment in muts_dict[muts_key]:
                    new_justification_count = 0
                    new_unjustification_count = 0
                    justified_peaks = []
                    unjustified_peaks = []
                    tokens = string.split(mutated_fragment[4])
                    if tokens[-2] == 'insertion':
                        continue
                    new_justification_count = 0
                    for unrationalized_peak_index in unrationalized_peak_indices:
                        position = int(tokens[-2]) + int(mutated_fragment[5])
                        wt_alteration = float(tokens[-1][1:-1])
                        for other_fragment_object in self.experiment.get_all_possible_fragment_objects():
                            nterm, cterm = other_fragment_object.get_nterm_index(), other_fragment_object.get_cterm_index()
                            if position >= nterm and position <= cterm:
                                # consider the alteration to have affected the fragment
                                wt = other_fragment_object.get_weight()
                                if abs((wt + wt_alteration) - self.parent.x[unrationalized_peak_index]) < (wt + wt_alteration) * self.my_filter_resolution:
                                    new_justification_count += 1
                                    justified_peaks.append([self.parent.x[unrationalized_peak_index], other_fragment_object])
                    # next see how many identified peaks would no longer be justified
                    # take the reaction profiles and see if it contains the position
                    for rationalized_peak in rationalized_peaks:
                        position = int(tokens[-2]) + int(mutated_fragment[5])
                        possible_fragments = rationalized_peak.get_possible_fragments()
                        if len(possible_fragments) != 1:
                            continue
                        nterm, cterm = possible_fragments[0].get_nterm_index(), possible_fragments[0].get_cterm_index()
                        if position >= nterm and position <= cterm:
                            new_unjustification_count += 1
                            unjustified_peaks.append([possible_fragments[0].get_weight(), possible_fragments[0]])

                    info = '%s %s %s %s (gain %s, lose %s)'%(mutated_fragment[0], mutated_fragment[5], mutated_fragment[6], mutated_fragment[4], new_justification_count, new_unjustification_count)
                    tree_store2[info] = [i, muts_key, info, new_justification_count, new_unjustification_count]
                    i += 1
                    for justified_peak in justified_peaks:
                        new_info = '%s justified by %s %s %s'%(justified_peak[0], justified_peak[1].get_nterm_index(), justified_peak[1].get_cterm_index(), justified_peak[1].get_chain())
                        tree_store3[new_info] = [i, info, new_info]
                        i += 1
                    for unjustified_peak in unjustified_peaks:
                        new_info = '%s is unjustified by %s %s %s'%(unjustified_peak[0], unjustified_peak[1].get_nterm_index(), unjustified_peak[1].get_cterm_index(), unjustified_peak[1].get_chain())
                        tree_store3[new_info] = [i, info, new_info]
                        i += 1
            # now build the tree with the stored nodes, adding the number of justified peaks to the top level node info
            for node1_key in tree_store1.keys():
                node1 = tree_store1[node1_key]
                # find the highest justification count
                max_just = 0
                for node2_key in tree_store2.keys():
                    node2 = tree_store2[node2_key]
                    if node2[1] == node1_key:
                        if node2[3] > max_just:
                            max_just = node2[3]
                max_unjust = 0
                for node2_key in tree_store2.keys():
                    node2 = tree_store2[node2_key]
                    if node2[1] == node1_key:
                        if node2[4] > max_unjust:
                            max_unjust = node2[4]
                node1[2] += ' (gain %s, lose %s)'%(max_just, max_unjust)
                trees[range_key][0].insert(0, 'root', node1_key, text=node1[2])
            for node2_key in tree_store2.keys():
                node2 = tree_store2[node2_key]
                trees[range_key][0].insert(node2[0], node2[1], node2[2], text=node2[2])
            for node3_key in tree_store3.keys():
                node3 = tree_store3[node3_key]
                trees[range_key][0].insert(node3[0], node3[1], node3[2], text=node3[2])

            tree_store1 = {}
            tree_store2 = {}
            tree_store3 = {}
            for mods_key in mods_dict.keys():
                tree_store1[mods_key] = [0, 'root', '(%s)  %s'%(len(mods_dict[mods_key]), mods_key)]
                i = 1
                for modified_fragment in mods_dict[mods_key]:
                    new_justification_count = 0
                    justified_peaks = []
                    new_unjustification_count = 0
                    unjustified_peaks = []
                    tokens = string.split(modified_fragment[4])
                    if tokens[-2] == 'insertion':
                        continue
                    for unrationalized_peak_index in unrationalized_peak_indices:
                        if tokens[-2] == 'n-terminus':
                            position = int(modified_fragment[5])
                        else:
                            position = int(tokens[-2]) + int(modified_fragment[5])
                        wt_alteration = float(tokens[-1][1:-1])
                        for other_fragment_object in self.experiment.get_all_possible_fragment_objects():
                            nterm, cterm = other_fragment_object.get_nterm_index(), other_fragment_object.get_cterm_index()
                            if position >= nterm and position <= cterm:
                                # consider the alteration to have affected the fragment
                                wt = other_fragment_object.get_weight()
                                if abs((wt + wt_alteration) - self.parent.x[unrationalized_peak_index]) < (wt + wt_alteration) * self.experiment.get_filter_resolution():
                                    new_justification_count += 1
                                    justified_peaks.append([self.parent.x[unrationalized_peak_index], other_fragment_object])
                    # next see how many identified peaks would no longer be justified
                    # take the reaction profiles and see if it contains the position
                    for rationalized_peak in rationalized_peaks:
                        position = int(tokens[-2]) + int(modified_fragment[5])
                        possible_fragments = rationalized_peak.get_possible_fragments()
                        if len(possible_fragments) != 1:
                            continue
                        nterm, cterm = possible_fragments[0].get_nterm_index(), possible_fragments[0].get_cterm_index()
                        if position >= nterm and position <= cterm:
                            new_unjustification_count += 1
                            unjustified_peaks.append([possible_fragments[0].get_weight(), possible_fragments[0]])

                    info = '%s %s %s %s (gain %s, lose %s)'%(modified_fragment[0], modified_fragment[5], modified_fragment[6], modified_fragment[4], new_justification_count, new_unjustification_count)
                    tree_store2[info] = [i, mods_key, info, new_justification_count, new_unjustification_count]
                    i += 1
                    for justified_peak in justified_peaks:
                        new_info = '%s is justified by %s %s %s %s'%(justified_peak[0], justified_peak[1].get_nterm_index(), justified_peak[1].get_cterm_index(), justified_peak[1].get_chain(), modified_fragment[4])
                        tree_store3[new_info] = [i, info, new_info]
                        i += 1
                    for unjustified_peak in unjustified_peaks:
                        new_info = '%s is unjustified by %s %s %s %s'%(unjustified_peak[0], unjustified_peak[1].get_nterm_index(), unjustified_peak[1].get_cterm_index(), unjustified_peak[1].get_chain(), modified_fragment[4])
                        tree_store3[new_info] = [i, info, new_info]
                        i += 1

            # now build the tree with the stored nodes, adding the number of justified peaks to the top level node info
            for node1_key in tree_store1.keys():
                node1 = tree_store1[node1_key]
                # find the highest justification count
                max_just = 0
                for node2_key in tree_store2.keys():
                    node2 = tree_store2[node2_key]
                    if node2[1] == node1_key:
                        if node2[3] > max_just:
                            max_just = node2[3]
                # find the highest unjustification count
                max_unjust = 0
                for node2_key in tree_store2.keys():
                    node2 = tree_store2[node2_key]
                    if node2[1] == node1_key:
                        if node2[4] > max_unjust:
                            max_unjust = node2[4]
                node1[2] += ' (gain %s, lose %s)'%(max_just, max_unjust)
                trees[range_key][1].insert(0, 'root', node1_key, text=node1[2], selectable=0)
            for node2_key in tree_store2.keys():
                node2 = tree_store2[node2_key]
                trees[range_key][1].insert(node2[0], node2[1], node2[2], text=node2[2], selectable=1)
            for node3_key in tree_store3.keys():
                node3 = tree_store3[node3_key]
                trees[range_key][1].insert(node3[0], node3[1], node3[2], text=node3[2], selectable=0)

            trees[range_key][0].pack(side=TOP, expand=YES, fill=BOTH)
            trees[range_key][1].pack(side=TOP, expand=YES, fill=BOTH)
            mut_button_box.pack(side=TOP, expand=NO, fill=NONE)
            mod_button_box.pack(side=TOP, expand=NO, fill=NONE)
            mut_scframe.pack(side=TOP, expand=YES, fill=BOTH)
            mod_scframe.pack(side=TOP, expand=YES, fill=BOTH)
            treepanes.pack(side=TOP, expand=YES, fill=BOTH)
        self.rtop.output_frame.pack(side=TOP, expand=YES, fill=BOTH)















