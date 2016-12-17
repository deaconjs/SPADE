import Tkinter
import os
import sys
import Pmw
import copy
import pickle
import string


class Viewer:
    def __init__(self, parent, sequences):
        self.parent_frame = parent.parent_frame
        self.sequences = sequences
        self.parent = parent
        self.reactions_top = Tkinter.Toplevel(self.parent_frame)
        self.reactions_top.title('Modification Manager')
        self.reactions_top.wm_transient(self.parent_frame)
        p_geo = self.parent_frame.geometry()
        args = string.split(p_geo, 'x')
        arg1 = args[0]
        args = string.split(args[1], '+')
        geometry_string = "%dx%d%+d%+d"%(500, 450, int(args[1])+25,int(args[2])+25) # width,height,x-offset,y-offset
        self.reactions_top.geometry(geometry_string)
        self.colormap = { 'color1': "#aaf",
                          'color2': "#afa",
                          'color3': "#aff",
                          'color4': "#faa",
                          'color5': "#faf",
                          'color6': "#ffa",
                          'color7': "#88b",
                          'color8': "#8b8",
                          'color9': "#8bb",
                          'color10':"#b88",
                          'color11':"#b8b",
                          'color12':"#bb8",
                          'color13':"#aab",
                          'color14':"#aba",
                          'color15':"#abb",
                          'color16':"#baa",
                          'color17':"#bab",
                          'color18':"#bba"}
        self.colorkeys = self.colormap.keys()
        self.colorkeys.sort()
        # available reagents
        self.reactions_dict = {}
        filename = './Tools/MassSpectrometry/reactions_dict.pkl'
        if os.path.exists(filename):
            reactions_file = open(filename, 'rb')
            self.reactions_dict = pickle.load(reactions_file)
            reactions_file.close()
        else:
            print 'no reactions file %s found'%(filename)
        keys = self.reactions_dict.keys()

        self.modifying_reagents = keys
        
        self.temp_reactant_dict = copy.deepcopy(self.reactions_dict)
        #frames
        self.tiptop_frame = Tkinter.Frame(self.reactions_top, borderwidth=2, relief='groove')
        self.top_frame = Tkinter.Frame(self.reactions_top)
        self.left_frame = Tkinter.Frame(self.top_frame,
                                        borderwidth = 2,
                                        relief = 'groove')
        self.right_frame = Tkinter.Frame(self.top_frame,
                                         borderwidth = 2,
                                         relief = 'groove')
        self.bottom_frame = Tkinter.Frame(self.reactions_top)

        # create a scrolledlistbox to hold the reactions
        rxn_keys = self.reactions_dict.keys()
        self.reactant_titles_dropdown = Pmw.ComboBox(self.left_frame,
                                                     label_text = 'Common Reactions',
                                                     labelpos = 'nw',
                                                     selectioncommand = self._view_reactant,
                                                     scrolledlist_items = rxn_keys)
        if len(rxn_keys) > 0:
            self.reactant_titles_dropdown.selectitem(rxn_keys[0])        
        self.reactant_titles_dropdown.pack(side='top', fill='x', expand=1)

        self.reactants_button_box = Pmw.ButtonBox(self.left_frame)
        self.reactants_save_button = self.reactants_button_box.add('Add', command = self._add_reactant)
        self.reactant_manager_cancel_button = self.reactants_button_box.add('Remove', command = self._remove_reactant)
        self.reactants_button_box.add('Import', command=self._import_reactants)
        self.reactants_button_box.pack(side='top', fill='x', expand=1)

        if len(rxn_keys) > 0:
            self.reactant_name_message = rxn_keys[0]
        else:
            self.reactant_name_message = ""
        self.reactant_name_text    = Tkinter.StringVar(self.right_frame)
        self.reactant_name_text.set(self.reactant_name_message)

        self.reactant_name_label   = Tkinter.Label(self.right_frame, textvariable=self.reactant_name_text, bg=self.colormap[self.colorkeys[0]])

        self.reactant_manager_button_box = Pmw.ButtonBox(self.bottom_frame)
        # support the import of reactants from the pickled mods dictionary
        file = open('./Tools/MassSpectrometry/mods.pkl', 'rb')
        mods = pickle.load(file)
        file.close()
        modskeys = mods.keys()
        modskeys.sort()
        all_keys = []
        for modskey in modskeys:
            newkeys = mods[modskey].keys()
            newkeys.sort()
            for key in newkeys:
                all_keys.append(key)
        self.which_mod = Pmw.ComboBoxDialog(self.reactions_top,
                                            title = 'Select a reaction to import',
                                            buttons = ('OK', 'Cancel'),
                                            defaultbutton = 'OK',
                                            combobox_labelpos = 'n',
                                            scrolledlist_items = (all_keys))
        self.which_mod.withdraw()
        self.reactant_manager_save_button = self.reactant_manager_button_box.add('Save', command = self._save_reactants)
        self.reactant_manager_save_button.config(state='disabled')
        
        self.reactant_manager_cancel_button = self.reactant_manager_button_box.add('Discard Changes', command = self._discard_reactant_manager_changes)
        self.reactant_manager_cancel_button = self.reactant_manager_button_box.add('Close', command = self._close_reactant_manager)

        self.reactant_added_weight_entry = Pmw.EntryField(self.right_frame,
                                                       labelpos = 'w',
                                                       label_text = 'Weight added to residue',
                                                       validate = {'validator':'real'},
                                                       modifiedcommand = self._update_temp_reactant_dict_weight)
        if len(rxn_keys) > 0:
            self.reactant_added_weight_entry.setvalue(self.reactions_dict[rxn_keys[0]]['added_weight'])
        
        self.reactant_target_AA_entry = Pmw.EntryField(self.right_frame,
                                                       labelpos = 'w',
                                                       label_text = 'Target Residue 1-letter codes',
                                                       validate = {'validator' : self._validate_one_letter_codes},
                                                       modifiedcommand = self._update_temp_reactant_dict_target_AA)
        if len(rxn_keys) > 0:
            self.reactant_target_AA_entry.setvalue(self.reactions_dict[rxn_keys[0]]['target_AA'])

        # now build a text box to hold the sequences for modification
        textfont = ('courier', 10)
        self.sequence_box = Pmw.ScrolledText(self.tiptop_frame, text_font=textfont, text_height=20)
        for sequence in self.sequences:
            self.sequence_box.appendtext('%s\n\n'%(sequence))
        self.add_reactant_popup = Tkinter.Menu(self.reactions_top, tearoff=0)
        self.sequence_box.component('text').bind("<Button-3>", self.popup_menu)

        self.sequence_box.pack(side='top', expand=1, fill='x')
        self.reactant_name_label.pack(side='top', expand=0, fill='x')
        self.reactant_added_weight_entry.pack(side='top', expand=1, fill='x')
        self.reactant_target_AA_entry.pack(side='top', expand=1, fill='x')
        self.tiptop_frame.pack(side='top', expand=1, fill='x')
        self.left_frame.pack(side='left', expand=1, fill='both')
        self.right_frame.pack(side='left', expand=1, fill='both')
        self.top_frame.pack(side='top', expand=1, fill='both')
        self.reactant_manager_button_box.pack(fill='x', expand=1)
        self.bottom_frame.pack(side='top', expand=1, fill='x')
        self._save_reactants()


    def popup_menu(self, event):
        # figure out which residue
        x = self.sequence_box.get("current", "end")
        name_table = {'A':'ala', 'C':'cys', 'D':'asp', 'E':'glu',
                      'F':'phe', 'G':'gly', 'H':'his', 'I':'ile',
                      'K':'lys', 'L':'leu', 'M':'met', 'N':'asn',
                      'P':'pro', 'Q':'gln', 'R':'arg', 'S':'ser',
                      'T':'thr', 'V':'val', 'W':'trp', 'Y':'tyr',
                      '0':'nterm', 'Z':'cterm'}
        keys = self.reactions_dict.keys()
        # nterm and cterm are treated as special amino acid types 0 (zero) and Z. In effect,
        # their reactive groups are considered a special amino acid type. The reactions_dict
        # has an extra entry for these reactions, so they shouldn't block the reactions of
        # sidechains
        chain = len(self.sequences) - len(string.split(x))
        position = len(self.sequences[chain]) - len(string.split(x)[0])
        for key in keys:
            self.add_reactant_popup.delete(0)
        for key in keys:
            name = '%s %s'%(x[0], key)
            if (x[0] in self.reactions_dict[key]['target_AA']) or ((position == 0) and ('0' in self.reactions_dict[key]['target_AA'])) or ((position == len(self.sequences[chain])) and ('Z' in self.reactions_dict[key]['target_AA'])):
                if (position == 0) and ('0' in self.reactions_dict[key]['target_AA']):
                    selection_key = '0'
                elif ((position == len(self.sequences[chain])) and ('Z' in self.reactions_dict[key]['target_AA'])):
                    selection_key = 'Z'
                else:
                    selection_key = x[0]
                submenu = Tkinter.Menu(self.add_reactant_popup, tearoff=0)
                added_weight = self.reactions_dict[key]['added_weight']
                target_AA    = self.reactions_dict[key]['target_AA']
                c_lambda = lambda c=chain, p=position, w=added_weight, k=key: self._react(c, p, w, k)
                submenu.add_command(label="%s this %s"%(key, name_table[x[0]]), command=c_lambda)
                d_lambda = lambda p=target_AA, w=added_weight, k=selection_key: self._react_all(p, w, k)
                submenu.add_command(label="%s all %s"%(key, self.reactions_dict[key]['target_AA']), command=d_lambda)
                self.add_reactant_popup.add_cascade(label=name, menu=submenu)
                #self.add_reactant_popup.add_separator()

        self.add_reactant_popup.post(event.x_root, event.y_root)

    def _react(self, chain, position, weight, reactant):
        # first figure out what color to make the background
        rxn_keys = self.reactions_dict.keys()
        index = 0
        for i in range(len(rxn_keys)):
            if reactant == rxn_keys[i]:
                index = i
                break
        tag_color = self.colormap[self.colorkeys[index%len(self.colorkeys)]]

        # create a list to add to the Experiment object's list
        # iterate through the chains and residues, collecting tags
        tags = {}
        for i in range(len(self.sequences)):
            sequence = self.sequences[i]
            tags[sequence] = []
            if chain == i:
                tags[sequence].append([int(position), float(weight)])
                print 'adding tag at %s'%('%s.%s'%(i, position))
                if i == 0:
                    i = 1
                else:
                    i = 2*i+1
                self.sequence_box.tag_add(reactant, '%d.%d'%(i, position))
                
        self.sequence_box.tag_config(reactant, background=tag_color)
        self.parent.update_modification_weights(tags) # update ProjectManager

    def _react_all(self, types, weight, reactant):
        # first figure out what color to make the background
        rxn_keys = self.reactions_dict.keys()
        index = 0
        for i in range(len(rxn_keys)):
            if reactant == rxn_keys[i]:
                index = i
                break
        tag_color = self.colormap[self.colorkeys[index]]

        # create a list to add to the Experiment object's list
        # iterate through the chains and residues, collecting tags
        tags = {}
        for i in range(len(self.sequences)):
            sequence = self.sequences[i]
            tags[sequence] = []
            for j in range(len(sequence)):
                rez = sequence[j]
                c = 1
                if rez in types:
                    tags[sequence].append([j, weight])
                    if i == 0:
                        c = 1
                    else:
                        c = 2*i+1
                    print 'adding tag at %s'%('%s.%s'%(c,j))
                    self.sequence_box.tag_add(reactant, '%d.%d'%(c, j))
        self.sequence_box.tag_config(reactant, background=tag_color)
        self.parent.update_modification_weights(tags)
        
    def _add_reactant(self):
        self.reactant_name_dialog = Pmw.PromptDialog(self.reactions_top,
                                       title = 'New Reactant',
                                       label_text = 'Reactant Title',
                                       entryfield_labelpos = 'n',
                                       defaultbutton = 0,
                                       buttons = ('OK', 'Cancel'),
                                       command = self._add_reactant_callback)
        
    def _add_reactant_callback(self, button_call):
        new_key = self.reactant_name_dialog.get()
        if len(new_key) > 0 and button_call == 'OK':
            self.temp_reactant_dict[new_key] = {'target_AA':"", 'added_weight':0}
            rxn_keys = self.temp_reactant_dict.keys()
            rxn_keys.sort()
            self.reactant_titles_dropdown.setlist(rxn_keys)
            self.reactant_titles_dropdown.selectitem(new_key)
            self.reactant_added_weight_entry.setvalue(self.temp_reactant_dict[new_key]['added_weight'])
            self.reactant_target_AA_entry.setvalue(self.temp_reactant_dict[new_key]['target_AA'])
        self.reactant_name_dialog.destroy()
        self.reactant_manager_save_button.config(state='normal')

    def _remove_reactant(self):
        new_key = self.reactant_name_text.get()
        # dialog popup to confirm removal
        if len(new_key) > 0:
            self.removal_confirmation_dialog = Pmw.MessageDialog(self.reactions_top,
                                                        title = 'Confirm removal',
                                                        message_text = 'deleting reactant %s'%(new_key),
                                                        command = self._remove_reactant_callback,
                                                        buttons = ['Confirm', 'Cancel'])
        
    def _remove_reactant_callback(self, button_call):
        new_key = self.reactant_name_text.get()
        if len(new_key) > 0:
            if button_call == 'Confirm':
                del self.temp_reactant_dict[new_key]
        self.removal_confirmation_dialog.destroy()
        rxn_keys = self.temp_reactant_dict.keys()
        rxn_keys.sort()
        self.reactant_titles_dropdown.setlist(rxn_keys)
        self.reactant_titles_dropdown.selectitem(rxn_keys[0])
        self.reactant_added_weight_entry.setvalue(self.temp_reactant_dict[rxn_keys[0]]['added_weight'])
        self.reactant_target_AA_entry.setvalue(self.temp_reactant_dict[rxn_keys[0]]['target_AA'])
        self.reactant_manager_save_button.config(state='normal')

    def _discard_reactant_manager_changes(self):
        """ discards all changes since the last save
        """
        self.temp_reactant_dict = copy.deepcopy(self.reactions_dict)
        self.reactant_added_weight_entry.setvalue(self.reactions_dict[self.reactant_name_message]['added_weight'])
        self.reactant_target_AA_entry.setvalue(self.reactions_dict[self.reactant_name_message]['target_AA'])
        self.reactant_manager_save_button.config(state='disabled')
    
    def _close_reactant_manager(self):
        self.reactions_top.destroy()

    def _view_reactant(self, reactant):
        try:
            self.temp_reactant_dict[reactant]
        except KeyError:
            print 'no such reactant'
            return
        else:
            rxn_keys = self.reactions_dict.keys()
            index = 0
            for i in range(len(rxn_keys)):
                if reactant == rxn_keys[i]:
                    index = i
                    break
            self.reactant_name_message = reactant # this is a Label... change bkg color
            self.reactant_name_label.config(bg=self.colormap[self.colorkeys[index]])
            self.reactant_name_text.set(self.reactant_name_message)
            self.reactant_target_AA_entry.setvalue(self.temp_reactant_dict[reactant]['target_AA'])
            self.reactant_added_weight_entry.setvalue(self.temp_reactant_dict[reactant]['added_weight'])

    def _validate_one_letter_codes(self, letters):
        for i in range(len(letters)):
            # includes '0' for nterm, 'Z' for cterm
            if letters[i] not in ['0', 'Z', 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']:
                return Pmw.PARTIAL
        else:
            return Pmw.OK
        
    def _update_temp_reactant_dict_weight(self):
        reactant_title = self.reactant_titles_dropdown.get()
        t2 = self.temp_reactant_dict[reactant_title]['added_weight']
        v2 = self.reactant_added_weight_entry.getvalue()
        if v2 != t2:
            self.temp_reactant_dict[reactant_title]['added_weight'] =  self.reactant_added_weight_entry.getvalue()
            self.reactant_manager_save_button.config(state='normal')
    
    def _update_temp_reactant_dict_target_AA(self):
        reactant_title = self.reactant_titles_dropdown.get ()
        t1 = self.temp_reactant_dict[reactant_title]['target_AA']
        v1 = self.reactant_target_AA_entry.getvalue()
        if v1 != t1:
            self.temp_reactant_dict[reactant_title]['target_AA'] = self.reactant_target_AA_entry.getvalue()
            self.reactant_manager_save_button.config(state='normal')

    def _save_reactants(self):
        self.reactions_dict = copy.deepcopy(self.temp_reactant_dict)
        reactions_file = open('./Tools/MassSpectrometry/reactions_dict.pkl', 'wb')
        pickle.dump(self.reactions_dict, reactions_file)
        self.reactant_manager_save_button.config(state='disabled')

    def _import_reactants(self):
        result = self.which_mod.activate()
        if result == 'OK':
            new_key = self.which_mod.get()
        file = open('./Tools/MassSpectrometry/mods.pkl', 'rb')
        mods = pickle.load(file)
        file.close()
        modskeys = mods.keys()
        modskeys.sort()
        specificities = ""
        modweight = 0
        for mod_key in modskeys:
            if len(mod_key) == 1:
                if new_key in mods[mod_key].keys():
                    specificities += mod_key
                    modweight = mods[mod_key][new_key]

        if len(new_key) > 0:
            self.temp_reactant_dict[new_key] = {'target_AA':specificities, 'added_weight':modweight}
            rxn_keys = self.temp_reactant_dict.keys()
            rxn_keys.sort()
            self.reactant_titles_dropdown.setlist(rxn_keys)
            self.reactant_titles_dropdown.selectitem(new_key)
            self.reactant_added_weight_entry.setvalue(self.temp_reactant_dict[new_key]['added_weight'])
            self.reactant_target_AA_entry.setvalue(self.temp_reactant_dict[new_key]['target_AA'])
        self.reactant_manager_save_button.config(state='normal')


