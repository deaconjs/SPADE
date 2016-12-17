import Tkinter
import os
import pickle
import sys
sys.path.append('Dependencies')
import Pmw
sys.path.append('Tools/Widgets')
import FunctionalPmwCounter
import RAVE
import copy
sys.path.append('Applications/Modification')
import FragmentViewer
reload(FragmentViewer)


class Viewer:
    def __init__(self, parent, sequences, modification_tags):
        self.parent = parent
        self.sequences = sequences
        self.modification_tags = modification_tags
        self.expcond_top = Tkinter.Toplevel(self.parent)
        self.expcond_top.title('Experimental Condition Suggestions')
        self.expcond_top.wm_transient(self.parent)
        self.expcond_top.top_width=800
        self.expcond_top.top_height=500
        geometry_string = "%dx%d%+d%+d"%(self.expcond_top.top_width,self.expcond_top.top_height,100,100) # width,height,x-offset,y-offset
        self.expcond_top.geometry(geometry_string)


        # frame1 is on the left
        self.frame1 = Tkinter.Frame(self.expcond_top, relief='raised',bd=2)
        self.frame1a = Tkinter.Frame(self.frame1, relief='groove',bd=2)
        self.frame1b = Tkinter.Frame(self.frame1, relief='groove',bd=2)
        self.frame1c = Tkinter.Frame(self.frame1, relief='groove',bd=2)
        self.frame2 = Tkinter.Frame(self.expcond_top, relief='raised',bd=2)

        # available reagents
        self.reactions_dict = {}
        filename = './Tools/MassSpectrometry/reactions_dict.pkl'
        if os.path.exists(filename):
            reactions_file = open(filename, 'rb')
            self.reactions_dict = pickle.load(reactions_file)
            reactions_file.close()
        else:
            print 'no reactions file %s found'%(filename)

        # put the available reagents into a scrolled frame
        self.reagents_sframe = Pmw.ScrolledFrame(self.frame1a, usehullsize=1, hull_width=200, hull_height=200)
        rsframe_interior = self.reagents_sframe.interior()

        # frame1a has available reagents (radiobutton)
        self.reagents_radio = Pmw.RadioSelect(rsframe_interior,
                                              buttontype = 'radiobutton',
                                              orient = 'vertical',
                                              labelpos = 'nw',
                                              label_text = 'Available\n Reagents',
                                              selectmode='single')
        self.reagents_radio.pack(side='top', expand = 0, fill='none')
        self.reagents_sframe.pack(side='top', expand=1, fill='both')

        keys = self.reactions_dict.keys()
        keys.sort()
        keys.reverse()
        for key in keys:
            self.reagents_radio.add(key)
        self.reagents_radio.invoke(keys[0])


        # available proteases
        self.proteases_dict = {}
        filename = './Tools/MassSpectrometry/protease_dict.pkl'
        if os.path.exists(filename):
            protease_file = open(filename, 'rb')
            self.proteases_dict = pickle.load(protease_file)
            protease_file.close()
        else:
            print 'no protease file %s found'%(filename)

        self.protease_sframe = Pmw.ScrolledFrame(self.frame1b, usehullsize=1, hull_width=200, hull_height=200)
        psframe_interior = self.protease_sframe.interior()

        self.protease_radio = Pmw.RadioSelect(psframe_interior,
                                              buttontype = 'checkbutton',
                                              orient = 'vertical',
                                              labelpos = 'nw',
                                              label_text = 'Available\n Proteases',
                                              selectmode='multiple')
        self.protease_radio.pack(side='top', expand = 0, fill='none')
        self.protease_sframe.pack(side='top', expand = 1, fill='both')
        keys = self.proteases_dict.keys()
        keys.sort()
        keys.reverse()
        for key in keys:
            self.protease_radio.add(key)
        self.protease_radio.invoke(keys[0])
        
        # frame1c has the missed cut sites counter
        self.missed_cut_sites_counter = FunctionalPmwCounter.functional_pmw_counter(self.frame1c,
                                                    'n',
                                                    'Missed cut sites tolerated:',
                                                    2,
                                                    {'validator':'integer', 'min':0, 'max':99},
                                                    'integer',
                                                    1)

        self.missed_cut_sites_counter.pack(side='top', expand=1, fill='both')

        # frame2 gets the report and buttonbox
        self.report_textbox = Pmw.ScrolledText(self.frame2,usehullsize = 1,hull_width = 400,hull_height = 350)
        self.report_textbox.pack(side='top', fill='both', expand=1)

        self.report_text = self.report_textbox.component('text')        

        self.buttonBox = Pmw.ButtonBox(self.frame2)
        self.buttonBox.add('Calculate', command = self.calculate_report)
        self.buttonBox.add('Cancel', command = self.cancel_calculation)
        self.buttonBox.alignbuttons()

        self.buttonBox.pack(side='top', fill='x', expand=0)

        self.frame1a.pack(side='top', expand=1, fill='both')
        self.frame1b.pack(side='top', expand=1, fill='both')
        self.frame1c.pack(side='top', expand=1, fill='both')
        self.frame1.pack(side='left', expand=0, fill='y')
        self.frame2.pack(side='left', expand=1, fill='both')
        self.old_stdout = None


    def cancel_calculation(self):
        if self.old_stdout:
            sys.stdout = self.old_stdout        # reset
        self.expcond_top.destroy()

    def write(self, text):
        """ write text to the outbox of the codebox. This functions is present
            to override stdout. Use python's normal print() command instead.
        """
        self.report_text.configure(state='normal')
        self.report_text.insert('end', text)
        self.report_text.configure(state='disabled')
        self.report_text.see('end')
        self.report_text.update_idletasks()

    def view_fragments(self, fragments, reagent_specificity=None, protease_specificity=None):
        self.fragview_top = Tkinter.Toplevel(self.expcond_top)
        self.fragview_top.title('Fragment Viewer')
        self.fragview_top.wm_transient(self.expcond_top)
        self.fragview_top.top_width=700
        self.fragview_top.top_height=500
        geometry_string = "%dx%d%+d%+d" %(self.fragview_top.top_width,self.fragview_top.top_height,200,200) # width,height,x-offset,y-offset
        self.fragview_top.geometry(geometry_string)

        # sort the fragments first by nterm then by cterm eg 1_5 comes before 2_5
        frag_index_dict = {}
        for frag in fragments:
            key = '%s_%s'%(frag.get_nterm_index(), frag.get_cterm_index())
            frag_index_dict[key] = frag
            
        keys = frag_index_dict.keys()
        keys.sort()


        frag_list = []
        for key in keys:
            #frag = frag_index_dict[key]
            #try:
            #    frag_list[frag.get_chain()]
            #except IndexError:
            #    while len(frag_list) < frag.get_chain()+1:
            #        frag_list.append([])
            #frag_list[int(frag.get_chain())].append(frag)#[frag.get_nterm_index(), frag.get_cterm_index(), frag.get_sequence(), frag.get_chain(), frag.get_color()])
            frag_list.append(frag_index_dict[key])
            
        chain_names = []
        for i in range(len(frag_list)):
            chain_names.append('%s'%(i))

        reagent_target_lists = None
        if reagent_specificity:
            reagent_target_lists = []
            for sequence in self.sequences:
                reagent_target_lists.append([])
                i = 0
                for aa in sequence:
                    if aa in reagent_specificity:
                        reagent_target_lists[-1].append(i)
                    i += 1

        protease_target_lists = None
        if protease_specificity:
            protease_target_lists = []
            for sequence in self.sequences:
                protease_target_lists.append([])
                i = 0
                for aa in sequence:
                    if aa in protease_specificity:
                        protease_target_lists[-1].append(i)
                    i += 1
            
        # build a menu for controlling the view of reagent specificity vs protease specificity vs both

        self.fragview_top.fview = FragmentViewer.Viewer(self.fragview_top, self.fragview_top, self.fragview_top.top_width, self.fragview_top.top_height, self.sequences, frag_list, reagent_target_lists, protease_target_lists)

        self.fragview_top.menuBar = Pmw.MenuBar(self.fragview_top, hull_relief = 'raised',hull_borderwidth = 2)
        self.fragview_top.menuBar.addmenu('vertical lines', 'vertical lines')
        self.fragview_top.menuBar.addmenuitem('vertical lines', 'command', label='reagent specificity', command=self.fragview_top.fview.draw_vertical_lines_1)
        self.fragview_top.menuBar.addmenuitem('vertical lines', 'command', label='protease specificity', command=self.fragview_top.fview.draw_vertical_lines_2)
        self.fragview_top.menuBar.addmenuitem('vertical lines', 'command', label='both', command=self.fragview_top.fview.draw_both_vertical_lines)
        self.fragview_top.menuBar.addmenuitem('vertical lines', 'command', label='neither', command=self.fragview_top.fview.remove_vertical_lines)
        self.fragview_top.menuBar.pack(fill='x', side='top', expand=0)

        self.fragview_top.fview.pack(side='top', fill='both', expand=1)
        Tkinter.Label(self.fragview_top, text='fragments without overlap', fg='black').pack(side='top', expand=0, fill='none')
        Tkinter.Label(self.fragview_top, text='unmodified fragment overlap with an unmodified fragment', fg='red').pack(side='top', expand=0, fill='none')
        Tkinter.Label(self.fragview_top, text='modified fragment overlap with an unmodified fragment', fg='blue').pack(side='top', expand=0, fill='none')
        Tkinter.Label(self.fragview_top, text='unmodified fragment overlap with a modified fragment', fg='purple').pack(side='top', expand=0, fill='none')
        Tkinter.Label(self.fragview_top, text='modified fragment overlap with a modified fragment', fg='green').pack(side='top', expand=0, fill='none')

    def calculate_report(self):
        # capture stdout, saving the old for resetting later
        self.old_stdout = sys.stdout
        sys.stdout = self

        self.report_textbox.clear()

        reagent = self.reagents_radio.getcurselection()

        proteases = []
        p = self.protease_radio.getcurselection()
        for prot in p:
            proteases.append(prot)

        missed_cuts = self.missed_cut_sites_counter.get()

        print 'reagent:        %s'%(reagent)
        print 'proteases:    ',
        for protease in proteases:
            print '%s, '%(protease),
        print '\nmissed cuts:  %s'%(missed_cuts)

        # set the experimental parameters
        exp = RAVE.MassSpecExperiment('exp', self.modification_tags)
        exp.set_modifying_reagent(reagent)
        exp.set_missed_cut_site_tolerance(missed_cuts)
        exp.set_digestion_type('exhaustive')
        exp.set_protein_sequences(self.sequences)

        print 'mass range from %s to %s'%(exp.spectrum_start, exp.spectrum_end)

        rxn_count = 0
        reaction_specificity = self.reactions_dict[reagent]['target_AA']
        reaction_weight      = float(self.reactions_dict[reagent]['added_weight'])

        # rxn_indices is a list of dictionaries, one for each sequence
        # where keys are indices of reactive residues, and 1 or 0 is stored to
        # represent whether the index has been blocked
        rxn_indices = []
        for seq in self.sequences:
            rxn_indices.append({})              # one for each sequence
            i = 0
            for res in seq:
                if res in reaction_specificity:
                    rxn_count += 1              # get the total number of reactive residues
                    rxn_indices[-1][i] = 0      # initialize a container to monitor them
                i += 1

        print '%s %s reactive residues available at indices:'%(rxn_count, reaction_specificity)
        i = 0
        for indexset in rxn_indices:
            keys = indexset.keys()
            keys.sort()
            print 'chain %s - '%(i),
            for ind in keys:
                print ind,
            i += 1
            print ''

        # start the calculation
        viewer_buttons = []
        for protease in proteases:
            print '\n\n%s:'%(protease)
            protease_specificity = self.proteases_dict[protease]
            exp.set_protease(protease)
            
            # calculate the single-reaction, n-missed cut site fragments
            protease_reagent_overlap = 0
            for r_spec in reaction_specificity:
                if r_spec in protease_specificity:
                    protease_reagent_overlap = 1

            exp._calculate_single_reaction_fragments_dictionary()
            
            # store 1 in rxn_indices if a fragment covers the index
            for frag in exp.fragment_objects:
                c = frag.get_cterm_index()
                n = frag.get_nterm_index()
                chain = frag.get_chain()
                keys = rxn_indices[chain].keys()
                keys.sort()
                for key in keys:
                    if protease_reagent_overlap:
                        if key >n and key < c:
                            rxn_indices[chain][key] = 1
                    else:
                        if key >n and key <= c:
                            rxn_indices[chain][key] = 1

            # make available_indices and unavailable_indices lists
            cnt = 0
            cnt2 = 0
            available_indices = []
            unavailable_indices = []
            for list in rxn_indices:
                available_indices.append([])
                unavailable_indices.append([])
                keys = list.keys()
                keys.sort()
                for key in keys:
                    cnt += list[key]
                    if list[key] and key not in available_indices[-1]:
                        available_indices[-1].append(key)
                        cnt2 += 1
                    elif not list[key] and key not in unavailable_indices[-1]:
                        unavailable_indices[-1].append(key)
            print '%s single-reactive-site fragments cover %s reaction indices:'%(len(exp.fragment_objects), cnt)
            
            i = 0
            for indexset in available_indices:
                print 'chain %s - '%(i),
                for ind in indexset:
                    print ind,
                print ''
                i += 1
            print ''
            
            # reset the reaction index-presence/absence monitor
            j = 0
            for seq in self.sequences:
                i = 0
                for res in seq:
                    if res in reaction_specificity:
                        rxn_indices[j][i] = 0
                    i += 1
                j += 1                    

            exp._calculate_all_possible_fragments_dictionary()
            all_frags_keys = exp.all_fragment_dict.keys()
            
            # check which frags (and their modifications) are obscured by other unmodified fragments
            # obscured_frags ends up being a list the size of the Experiment's all_fragments list with
            #  0 if not obscured by other fragments
            print 'testing unmodified and unmodified versions for overlap with unmodified peaks'
            obscured_frags = []
            for frag in exp.fragment_objects:
                frag_weight = frag.get_weight()
                break1 = 0
                break2 = 0
                for key in all_frags_keys:
                    if abs(frag_weight - key) < exp.get_filter_resolution()*key and abs(frag_weight-key) != 0.0:
                        print 'unmodified reactive fragment %s overlaps with fragment %s'%(frag_weight, key)
                        break1 = 1
                    if abs(frag_weight + reaction_weight - key) < exp.get_filter_resolution()*key and abs(frag_weight-key) != 0.0:
                        print 'modified reactive fragment %s overlaps with fragment %s'%(frag_weight+reaction_weight, key)
                        break2 = 1
                    if break2:
                        frag.set_color('blue')
                        obscured_frags.append(1)
                    elif break1:
                        frag.set_color('red')
                        obscured_frags.append(1)
                    if break1 or break2:
                        break
                else:
                    obscured_frags.append(0)        # not obscured

            unobscured_frag_count = 0
            for obs in obscured_frags:
                if not obs:
                    unobscured_frag_count += 1

            # count the remaining number of reaction sites on the unobscured fragments
            # rxn_indices gets 1 if the reactive site is not obscured
            i = 0
            for frag in exp.fragment_objects:
                if not obscured_frags[i]:       # only check the unobscured ones
                    c = frag.get_cterm_index()
                    n = frag.get_nterm_index()
                    list = rxn_indices[frag.get_chain()]
                    keys = list.keys()      # these keys are the list of indices for reactive rez
                    keys.sort()
                    for key in keys:
                        if protease_reagent_overlap:
                            if key >n and key < c:
                                list[key] = 1
                        else:
                            if key >n and key <= c:
                                list[key] = 1
                i += 1

            unobscured_coverage_count = 0
            for list in rxn_indices:
                keys = list.keys()
                for key in keys:
                    unobscured_coverage_count += list[key]

            print '%s reactive indices from %s fragments stand out from the %s unmodified fragments\n'%(unobscured_coverage_count, unobscured_frag_count, len(all_frags_keys))

            # reset the reaction index-presence/absence monitor
            j = 0
            for seq in self.sequences:
                i = 0
                for res in seq:
                    if res in reaction_specificity:
                        rxn_indices[j][i] = 0
                    i += 1
                j += 1

            print 'testing unmodified and modified versions for overlap with modified peaks'
            # test the frags against the modified versions of all available fragments
            i = 0
            mod_obscured_frags = []
            test_against = 0
            for frag in exp.fragment_objects:
                if not obscured_frags[i]:
                    test_against = 0
                    frag_weight = frag.get_weight()
                    breakout = 0
                    break1 = 0
                    break2 = 0
                    for key in all_frags_keys:
                        sec_frag = exp.all_fragment_dict[key]
                        sec_seq = sec_frag['sequence']
                        for frag_object in exp.all_fragment_objects:
                            if sec_seq == frag_object.get_sequence():
                                sec_frag_object = frag_object
                        rxn_count = 0
                        for res in sec_seq:
                            if res in reaction_specificity:
                                rxn_count += 1
                        if rxn_count >= 1:
                            for j in range(rxn_count):
                                added_weight = j * reaction_weight
                                sec_mod_weight = sec_frag_object.get_weight() + added_weight
                                test_against += 1

                                if abs(frag_weight - sec_mod_weight) < 7.0 and frag_weight != sec_mod_weight:
                                    print 'unmodified reactive fragment %s overlaps with modified fragment %s'%(frag_weight, sec_mod_weight)
                                    break1 = 1
                                if abs(frag_weight + reaction_weight - sec_mod_weight) < 7.0 and frag_weight+reaction_weight != sec_mod_weight:
                                    print 'modified reactive fragment %s overlaps with modified fragment %s'%(frag_weight+reaction_weight, sec_mod_weight)
                                    break2 = 1
                                if break2:
                                    frag.set_color('green')
                                    mod_obscured_frags.append(1)
                                elif break1:
                                    frag.set_color('purple')
                                    mod_obscured_frags.append(1)
                                if break1 or break2:
                                    breakout = 1
                                    break
                                
                        if breakout:
                            break
                    else:
                        mod_obscured_frags.append(0)
                else:
                    mod_obscured_frags.append(1)
                i += 1

            print 'tested for overlap with %s modified fragments, from %s total frags'%(test_against, len(all_frags_keys))

            # count the remaining number of reaction sites on the probable fragments
            # considering modifications of all_fragments
            i = 0
            for frag in exp.fragment_objects:
                if not obscured_frags[i] and not mod_obscured_frags[i]:
                    c = frag.get_cterm_index()
                    n = frag.get_nterm_index()
                    list = rxn_indices[frag.get_chain()]
                    keys = list.keys()      # these keys are the list of indices for reactive rez
                    keys.sort()
                    for key in keys:
                        if protease_reagent_overlap:
                            if key >n and key < c:
                                list[key] = 1
                        else:
                            if key >n and key <= c:
                                list[key] = 1
                i += 1

            cnt = 0
            available_indices = []
            j = 0
            for seq in self.sequences:
                available_indices.append([])
                list = rxn_indices[j]
                keys = list.keys()
                keys.sort()
                for key in keys:
                    cnt += list[key]
                    if list[key]:
                        available_indices[-1].append(key)
                j += 1

            print '%s of the reactive sites should be represented without overlap:'%(cnt)
            j = 0
            for seq in self.sequences:
                print 'chain %s '%(j),
                for ind in available_indices[j]:
                    print ind,
                j += 1
                print

            c_lambda = lambda x=copy.deepcopy(exp.fragment_objects), y = copy.deepcopy(reaction_specificity), z = copy.deepcopy(protease_specificity): self.view_fragments(x, y, z)
            viewer_button = Tkinter.Button(self.report_text, text='View', command=c_lambda)
            viewer_button.pack(side='top', expand=0, fill='none')
            viewer_buttons.append(viewer_button)
            self.report_text.window_create('end', window=viewer_button)
            self.report_text.insert('end', '\n\n')

            # calculate for follow-up with other proteases
            useful_second_indices = []
            for protease2 in proteases:
                if protease2 == protease:
                    continue
                print '\n\n%s -> %s'%(protease, protease2)
                second_protease_specificity = self.proteases_dict[protease2]
                frags = exp.fragment_objects

                # digest the fragments, collecting all within the molecular weight range that have
                # only one modification on them in new_fragset and all of the fragments in new_fragset_all
                new_fragset = []
                new_fragset_sequences = []
                new_fragset_all = []
                new_fragset_all_sequences = []
                for frag in frags:
                    seq = frag.get_sequence()
                    growing_seq = []
                    nterm       = frag.get_nterm_index()
                    last_nterm  = nterm
                    total_length = 0
                    length = 0
                    for aa in seq:
                        growing_seq.append(aa)
                        total_length += 1
                        length += 1
                        if aa in second_protease_specificity or total_length == len(seq):
                            # make sure only one modification
                            rxns = 0
                            for aa2 in growing_seq:
                                if aa2 in reaction_specificity:
                                    rxns += 1
                            # store and reset
                            new_frag = RAVE.Fragment(copy.copy(growing_seq), -1, last_nterm, last_nterm+length, None, frag.get_chain())
                            last_nterm = last_nterm+length
                            length = 0
                            growing_seq = []
                            if new_frag.get_weight() > 800.0 and new_frag.get_weight() < 10000.0:
                                if new_frag.get_sequence() not in new_fragset_all_sequences:
                                    new_fragset_all.append(new_frag)
                                    new_fragset_all_sequences.append(new_frag.get_sequence())
                                    if rxns == 1:
                                        if new_frag.get_sequence() not in new_fragset_sequences:
                                            new_fragset.append(new_frag)
                                            new_fragset_sequences.append(new_frag.get_sequence())

                print '%s single reaction fragments, among %s total'%(len(new_fragset), len(new_fragset_all))

                # collect the reaction indices and the fragments that cover them
                indexed_fragments = {}
                for j in range(len(self.sequences)):
                    seq = self.sequences[j]
                    rxn_index_keys = rxn_indices[j].keys()
                    rxn_index_keys.sort()
                    for key in rxn_index_keys:
                        for new_frag in new_fragset:
                            nt = new_frag.get_nterm_index()
                            ct = new_frag.get_cterm_index()
                            if key > nt and key < ct:
                                if key not in indexed_fragments.keys():
                                    indexed_fragments[key] = [new_frag]
                                else:
                                    if new_frag not in indexed_fragments[key]:
                                        indexed_fragments[key].append(new_frag)

                idx_keys = indexed_fragments.keys()
                idx_keys.sort()
                frag_cnt = 0
                print 'reactive indices included:'
                for idx in idx_keys:
                    print idx,
                    frag_cnt += len(indexed_fragments[idx])
                print 'covered by %s fragments'%(frag_cnt)
                
                # make sure none of the other fragments overlap with the recovery ones or their modifications
                frags_to_view = []
                index_keys = indexed_fragments.keys()
                index_keys.sort()
                for index in index_keys:
                    frags_to_recover = indexed_fragments[index]
                    remaining_chances = len(frags_to_recover)
                    remaining_holder = []
                    for j in range(remaining_chances):
                        remaining_holder.append(0)
                        
                    for test_frag in frags_to_recover:
                        frag_weight = test_frag.get_weight()
                        j = 0
                        for new_fragment in new_fragset_all:
                            sec_weight = new_frag.get_weight()
                            if new_fragment.get_sequence() != test_frag.get_sequence():
                                if abs(frag_weight-sec_weight) < 7.0:
                                    print 'unmodified reactive fragment %s overlaps with unmodified fragment %s'%(frag_weight, sec_weight)
                                    test_frag.set_color('blue')
                                    remaining_chances -= 1
                                    remaining_holder[j] = 1
                                    break
                                elif abs(frag_weight+reaction_weight-sec_weight) < 7.0:
                                    print 'modified reactive fragment %s overlaps with unmodified fragment %s'%(frag_weight+reaction_weight, sec_weight)
                                    test_frag.set_color('red')
                                    remaining_chances -= 1
                                    remaining_holder[j] = 1
                                    break
                            j += 1
                    final_count = 0
                    for j in range(len(remaining_holder)):
                        if not remaining_holder[j]:
                            final_count += 1
                    #print '%s fragments containing blocked index %s are visible through unmodified fragments'%(final_count, index)

                    # now check modified versions
                    for new_fragment in new_fragset_all:
                        new_seq = new_fragment.get_sequence()
                        rxncnt = 0
                        for aa in new_seq:
                            if aa in reaction_specificity:
                                rxncnt += 1
                        if rxncnt == 0:
                            continue
                        for j in range(1,rxncnt):
                            added_weight = j * reaction_weight
                            mod_weight = new_fragment.get_weight() + added_weight
                            k = 0
                            for test_frag in frags_to_recover:
                                if test_frag.get_sequence() != new_fragment.get_sequence():
                                    if not remaining_holder[k]:
                                        if abs(mod_weight-test_frag.get_weight()) < 7.0:
                                            print 'unmodified reactive fragment %s overlaps with modified fragment %s'%(test_frag.get_weight(), mod_weight)
                                            test_frag.set_color('green')
                                            remaining_chances -= 1
                                            remaining_holder[k] = 1
                                            break
                                        elif abs(mod_weight-(test_frag.get_weight()+reaction_weight)) < 7.0:
                                            print 'modified reactive fragment %s overlaps with modified fragment %s'%(test_frag.get_weight()+reaction_weight, mod_weight)
                                            test_frag.set_color('purple')
                                            remaining_chances -= 1
                                            remaining_holder[k] = 1
                                            break
                                k += 1
                    for test_frag in frags_to_recover:
                        frags_to_view.append(test_frag)
                    final_count = 0
                    for j in range(len(remaining_holder)):
                        if not remaining_holder[j]:
                            final_count += 1
                    
                    #print '%s fragments containing blocked index %s are visible'%(final_count, index)
                    #print final_count, index, available_indices
                    if final_count:
                        useful_second_indices.append(index)
                    if final_count and index not in available_indices:
                        print 'recovered index %s in %s fragments'%(index, final_count)
                print '%s frags to view, covering indices %s'%(len(frags_to_view), useful_second_indices)
                
                c_lambda = lambda x=copy.deepcopy(frags_to_view), y = copy.deepcopy(reaction_specificity), z = copy.deepcopy(second_protease_specificity): self.view_fragments(x, y, z)
                viewer_button = Tkinter.Button(self.report_text, text='View', command=c_lambda)
                viewer_button.pack(side='top', expand=0, fill='none')
                viewer_buttons.append(viewer_button)
                self.report_text.window_create('end', window=viewer_button)
                self.report_text.insert('end', '\n\n')


                        
            # reset the reaction index-presence/absence monitor            
            for seq in self.sequences:
                i = 0
                for res in seq:
                    if res in reaction_specificity:
                        rxn_indices[-1][i] = 0
                    i += 1
            
        # reset stdout
        sys.stdout = self.old_stdout


