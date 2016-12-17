import Tkinter
from Tkinter import END
import string
import sys
sys.path.append('Dependencies')
import Pmw

class SystemSelectionDialog(Tkinter.Frame):
    def __init__(self, viewer, top, system, return_selected=0, return_list=[], see_shortcuts=0, exit_function=None):
        """ can be used as either a query tool that returns a single molecule,
            residue or atom, or as a general selection tool for changing the
            selection state of the same objects. return_selected can be 0 for general
            selection, or 'atom', 'residue', or 'molecule', for the query tool.
            see_shortcuts != 0 shows a variety of shortcuts to selection groups
            exit_function specifies a function to run when the dialog box is closed.
            This should be recreated to work without a user interface            
        """
        # shortcuts to add:
        # molecules    : chain A-N, nucleotide chains, proteins, ligands, water
        # amino acids  : A-Y, acidic, basic, neutral, aromatic, other (polar)
        # atom         : alpha C, beta C, backbone, sidechain, type, name

        # secondary    : beta, alpha, nonsecondary
        # data driven  : .features or .features combo threshold
        # distance_from: popup w/ named sets and single choice from below

        # set options  : replace, union, intersection, not
        # named sets   : save set as name, offer named sets
        # transparency control for selected (remembers old opacity with an OpacityVisitor)

        # molecules shortcut

        Tkinter.Frame.__init__(self, top)
        self.system = system
        self.top = top
        if exit_function != None:
            top.protocol("WM_DELETE_WINDOW",exit_function)
        self.viewer = viewer
        self.return_selected = return_selected
        self.return_list = return_list
        frame1  = Tkinter.Frame(top)
        frame2  = Tkinter.Frame(top)
        frame3  = Tkinter.Frame(top)

        self.mlabel_text = Tkinter.StringVar()
        self.mlabel = Tkinter.Label(frame1, textvariable=self.mlabel_text)
        self.mlabel_text.set('Molecules')      
        self.mlabel.pack()
        self.rlabel_text = Tkinter.StringVar()
        self.rlabel = Tkinter.Label(frame2, textvariable=self.rlabel_text)
        self.rlabel_text.set('Residues')
        self.rlabel.pack()
        self.alabel_text = Tkinter.StringVar()
        self.alabel = Tkinter.Label(frame3, textvariable=self.alabel_text)
        self.alabel_text.set('Atoms')
        self.alabel.pack()

        mframe1 = Tkinter.Frame(frame1)
        rframe1 = Tkinter.Frame(frame2)
        aframe1 = Tkinter.Frame(frame3)
        self.moleculeWorkScroll = Tkinter.Scrollbar(mframe1, orient='vertical')
        self.residueWorkScroll  = Tkinter.Scrollbar(rframe1, orient='vertical')
        self.atomWorkScroll     = Tkinter.Scrollbar(aframe1, orient='vertical')
        self.moleculeWorkBox    = Tkinter.Listbox(mframe1, selectmode='multiple', height=4, width=12, yscrollcommand=self.moleculeWorkScroll.set, exportselection=0)
        self.residueWorkBox     = Tkinter.Listbox(rframe1, selectmode='multiple', height=4, width=12, yscrollcommand=self.residueWorkScroll.set, exportselection=0)
        self.atomWorkBox        = Tkinter.Listbox(aframe1, selectmode='multiple', height=4, width=12, yscrollcommand=self.atomWorkScroll.set, exportselection=0)
        self.moleculeWorkScroll.config(command=self.moleculeWorkBox.yview)
        self.residueWorkScroll.config(command=self.residueWorkBox.yview)
        self.atomWorkScroll.config(command=self.atomWorkBox.yview)

        self.moleculeButtons = Pmw.ButtonBox(frame1, orient='vertical')
        self.residueButtons  = Pmw.ButtonBox(frame2, orient='vertical')
        self.atomButtons     = Pmw.ButtonBox(frame3, orient='vertical')
        
        for pchain in self.system.ProteinList:
            self.moleculeWorkBox.insert('end', 'Pro %s'%(pchain.chain_name))
        for nchain in self.system.NucleotideChainList:
            self.moleculeWorkBox.insert('end', 'Nuc %s'%(nchain.chain_name))
        for lig in self.system.LigandList:
            self.moleculeWorkBox.insert('end', 'Lig %s'%(lig.key))

        self.moleculeButtons.add('Select',   command=self._select_molecules)
        self.moleculeButtons.add('Deselect', command=self._deselect_molecules)
        self.moleculeButtons.add('Expand',   command=self._expand_molecules)
       
        self.residueButtons.add('Select',   command=self._select_residues)
        self.residueButtons.add('Deselect', command=self._deselect_residues)
        self.residueButtons.add('Expand',   command=self._expand_residues)
       
        self.atomButtons.add('Select',   command=self._select_atoms)
        self.atomButtons.add('Deselect', command=self._deselect_atoms)
        # turn off unnecessary select/deselect buttons
        if return_selected in ['atom', 'residue']:
            self.moleculeButtons.button('Select').config(state=DISABLED)
        if return_selected in ['atom', 'molecule']:
            self.residueButtons.button('Select').config(state=DISABLED)
        if return_selected in ['residue', 'molecule']:
            self.atomButtons.button('Select').config(state=DISABLED)
        if return_selected:
            self.moleculeButtons.button('Deselect').config(state=DISABLED)
            self.residueButtons.button('Deselect').config(state=DISABLED)
            self.atomButtons.button('Deselect').config(state=DISABLED)
        
        mframe2 = Tkinter.Frame(frame1)
        rframe2 = Tkinter.Frame(frame2)
        aframe2 = Tkinter.Frame(frame3)
        self.moleculeSelectionScroll = Tkinter.Scrollbar(mframe2, orient='vertical')
        self.residueSelectionScroll  = Tkinter.Scrollbar(rframe2, orient='vertical')
        self.atomSelectionScroll     = Tkinter.Scrollbar(aframe2, orient='vertical')
        self.moleculeSelectionText   = Tkinter.Text(mframe2, height=4, width=12, yscrollcommand=self.moleculeSelectionScroll.set)
        self.residueSelectionText    = Tkinter.Text(rframe2, height=4, width=12, yscrollcommand=self.residueSelectionScroll.set)
        self.atomSelectionText       = Tkinter.Text(aframe2, height=4, width=12, yscrollcommand=self.atomSelectionScroll.set)
        self.moleculeSelectionScroll.config(command=self.moleculeSelectionText.yview)
        self.residueSelectionScroll.config(command=self.residueSelectionText.yview)
        self.atomSelectionScroll.config(command=self.atomSelectionText.yview)

        for pchain in self.system.ProteinList:
            self.moleculeSelectionText.insert('end', 'Pro %s\n'%(pchain.chain_name))
        for nchain in self.system.NucleotideChainList:
            self.moleculeSelectionText.insert('end', 'Nuc %s\n'%(nchain.chain_name))
        for lig in self.system.LigandList:
            self.moleculeSelectionText.insert('end', 'Lig %s\n'%(lig.key))

        self.moleculeWorkScroll.pack(side='right', expand='yes', fill='y')
        self.residueWorkScroll.pack(side='right', expand='yes', fill='y')
        self.atomWorkScroll.pack(side='right', expand='yes', fill='y')
        self.moleculeWorkBox.pack(side='right', expand='yes', fill='both')
        self.residueWorkBox.pack(side='right', expand='yes', fill='both')
        self.atomWorkBox.pack(side='right', expand='yes', fill='both')
        
        self.moleculeSelectionScroll.pack(side='right', expand='yes', fill='y')
        self.residueSelectionScroll.pack(side='right', expand='yes', fill='y')
        self.atomSelectionScroll.pack(side='right', expand='yes', fill='y')
        self.moleculeSelectionText.pack(side='right', expand='yes', fill='both')
        self.residueSelectionText.pack(side='right', expand='yes', fill='both')
        self.atomSelectionText.pack(side='right', expand='yes', fill='both')
        
        mframe1.pack(side='top', expand='yes', fill='both')
        rframe1.pack(side='top', expand='yes', fill='both')
        aframe1.pack(side='top', expand='yes', fill='both')

        self.moleculeButtons.pack(side='top', expand='no', fill='x')
        self.residueButtons.pack(side='top', expand='no', fill='x')
        self.atomButtons.pack(side='top', expand='no', fill='x')

        mframe2.pack(side='top', expand='yes', fill='both')
        rframe2.pack(side='top', expand='yes', fill='both')
        aframe2.pack(side='top', expand='yes', fill='both')

        frame1.pack(side='left', expand='yes', fill='both')
        frame2.pack(side='left', expand='yes', fill='both')
        frame3.pack(side='left', expand='yes', fill='both')
        # if used as a query tool, don't need to start out w/ initial selection
        if not return_selected:
            self._update_selection_state()        
        
    def _update_selection_state(self):
        self.moleculeSelectionText.tag_delete('partly selected')
        self.moleculeSelectionText.tag_delete('fully selected')
        # pass through each of the lower boxes and mark those items that are selected
        str = self.moleculeSelectionText.get('1.0', END+'-1c')
        lines = string.split(str, '\n')
        index = 0
        for mol in self.system.MoleculeList:
            start = "%d.0"%(index+1)
            end = "%d.%d"%(index+1, len(lines[index]))
            selected_count = 0
            for atom in mol.atoms:
                if atom.selected == 1:
                    selected_count = selected_count + 1
            if selected_count == len(mol.atoms):
                self.moleculeSelectionText.tag_add('fully selected', start, end)
            elif selected_count == 0:
                pass
            elif selected_count < len(mol.atoms):
                self.moleculeSelectionText.tag_add('partly selected', start, end)
            index = index + 1
        tk_rgb = "#%02x%02x%02x" % (21, 95, 189)
        self.moleculeSelectionText.tag_config('fully selected', background=tk_rgb)
        self.moleculeSelectionText.tag_config('fully selected', foreground='white')
        tk_rgb = "#%02x%02x%02x" % (63, 240, 255)
        self.moleculeSelectionText.tag_config('partly selected', background=tk_rgb)
        self.moleculeSelectionText.tag_config('partly selected', foreground='black')

        self.residueSelectionText.tag_delete('partly selected')
        self.residueSelectionText.tag_delete('fully selected')
        str = self.residueSelectionText.get('1.0', END+'-1c')
        lines = string.split(str, '\n')
        res_label = self.rlabel_text.get()
        tokens = string.split(res_label)
        if len(tokens) == 2:
            tokens.append('')
       
        if len(tokens) > 1:         # if a polymer has been expanded
            index = 0
            print tokens
            for res in self.system.PolymerDict[tokens[2]].residues:
                start = "%d.0"%(index+1)
                end = "%d.%d"%(index+1, len(lines[index]))
                selected_count = 0
                for atom in res.atoms:
                    if atom.selected == 1:
                        selected_count = selected_count + 1
                if selected_count == len(res.atoms):
                    self.residueSelectionText.tag_add('fully selected', start, end)
                elif selected_count == 0:
                    pass
                else:       # only some are selected
                    self.residueSelectionText.tag_add('partly selected', start, end)
                index = index + 1
            tk_rgb = "#%02x%02x%02x" % (21, 95, 189)
            self.residueSelectionText.tag_config('fully selected', background=tk_rgb)
            self.residueSelectionText.tag_config('fully selected', foreground='white')
            tk_rgb = "#%02x%02x%02x" % (10, 47, 94)
            self.residueSelectionText.tag_config('partly selected', background=tk_rgb)
            self.residueSelectionText.tag_config('partly selected', foreground='white')
        
    def _select_molecules(self):
        selected_molecules = self._get_selected_molecules()
        if not self.return_selected:
            for mol in selected_molecules:
                mol.select()
            self.viewer.update_view()
            self._update_selection_state()
        else:
            if len(selected_molecules) == 1:
                self.top.destroy()
                self.return_list.append(selected_molecules[0])
                return self.return_list
            else:
                print 'Select just one molecule'
        
    def _deselect_molecules(self):
        for mol in self._get_selected_molecules():
            mol.deselect()
        self.viewer.update_view()
        self._update_selection_state()
        
    def _expand_molecules(self):
        mols = self._get_selected_molecules()
        # for now, just do the first one
        mol = mols[0]
        if mol.__module__ in ['MolecularComponents.classProtein', 'MolecularComponents.classNucleotideChain']:
            self.residueWorkBox.delete(0,END)
            self.residueSelectionText.delete(1.0,END)
            self.atomWorkBox.delete(0, END)
            self.atomSelectionText.delete(1.0, END)
            for res in mol.residues:
                self.residueWorkBox.insert(END, '%3s %s'%(res.res_type, res.res_number))
                self.residueSelectionText.insert(END, '%3s %s\n'%(res.res_type, res.res_number))
            if mol.__module__ == 'MolecularComponents.classProtein':
                self.alabel_text.set('Atoms')
                self.rlabel_text.set('Residues: pchain %s'%(mol.chain_name))
            elif mol.__module__ == 'MolecularComponents.classNucleotideChain':
                self.alabel_text.set('Atoms')
                self.rlabel_text.set('Residues: nchain %s'%(mol.chain_name))           
        elif mol.__module__ == 'MolecularComponents.classLigand':
            self.residueWorkBox.delete(0,END)
            self.residueSelectionText.delete(1.0,END)
            self.atomWorkBox.delete(0, END)
            self.atomSelectionText.delete(1.0,END)
            for atom in mol.atoms:
                self.atomWorkBox.insert(END, '%3s %s'%(atom.atom_type, atom.atom_number))
                self.atomSelectionText.insert(END, '%3s %s\n'%(atom.atom_type, atom.atom_number))
            self.rlabel_text.set('Residues')
            self.alabel_text.set('Atoms: ligand %s %s'%(mol.res_number, mol.chain_name))
        self._update_selection_state()
        
    def _get_selected_molecules(self):
        sel_list = self.moleculeWorkBox.curselection()
        items = []
        index = 0
        for item in [self.system.ProteinList, self.system.NucleotideChainList, self.system.LigandList]:
            for mol in item:
                if '%s'%(index) in sel_list:
                    items.append(mol)
                index = index + 1
        return items
    
    def _select_residues(self):
        selected_residues = self._get_selected_residues()
        if not self.return_selected:
            for res in selected_residues:
                res.select()
            self.viewer.update_view()
            self._update_selection_state()
        else:
            if len(selected_residues) == 1:
                self.top.destroy()
                self.return_list.append(selected_residues[0])
                return self.return_list
            else:
                print 'Select just one residue'
        
    def _deselect_residues(self):
        for res in self._get_selected_residues():
            res.deselect()
        self.viewer.update_view()
        self._update_selection_state()
        
    def _expand_residues(self):
        mols = self._get_selected_residues()
        # for now, just do the first one
        mol = mols[0]
        self.atomWorkBox.delete(0, END)
        self.atomSelectionText.delete(1.0, END)
        for atom in mol.atoms:
            self.atomWorkBox.insert(END, '%3s %s'%(atom.atom_number, atom.atom_type))
            self.atomSelectionText.insert(END, '%3s %s\n'%(atom.atom_number, atom.atom_type))
        if mol.__module__ == 'MolecularComponents.classAminoAcid':
            self.alabel_text.set('Atoms: aminoacid %s %s'%(mol.res_number, mol.chain_name))
        elif mol.__module__ == 'MolecularComponents.classNucleotide':
            self.alabel_text.set('Atoms: nucleotide %s %s'%(mol.res_number, mol.chain_name))
        self._update_selection_state()

    def _get_selected_residues(self):
        sel_list = self.residueWorkBox.curselection()
        chain_name = ''
        chain = None
        if len(string.split(self.rlabel_text.get())) == 3:
            chain_name = string.split(self.rlabel_text.get())[2]
        for item in [self.system.ProteinList, self.system.NucleotideChainList]:
            for mol in item:
                if string.strip(mol.chain_name) == string.strip(chain_name):
                    chain = mol
                    break
            if chain != None:
                break
        items = []
        index = 0
        for res in chain.residues:
            if '%s'%(index) in sel_list:
                items.append(res)
            index = index + 1
        return items

    def _select_atoms(self):
        atoms = self._get_selected_atoms()
        if not self.return_selected:
            for atom in atoms:
                atom.select()
            self.viewer.update_view()
            self._update_selection_state()
        else:
            if len(atoms) == 1:
                self.top.destroy()
                self.return_list.append(atoms[0])
                return self.return_list
            else:
                print 'Select just one atom'

    def _deselect_atoms(self):
        atoms = self._get_selected_atoms()
        for atom in atoms:
            atom.deselect()
        self.viewer.update_view()
        self._update_selection_state()

    def _get_selected_atoms(self):
        sel_list = self.atomWorkBox.curselection()
        chain_name = ''
        chain = None
        tokens = string.split(self.alabel_text.get())
        print tokens
        if len(tokens) == 4:
            chain_name = tokens[3]
        res_number = tokens[2]

        items = []
        if tokens[1] == 'aminoacid':
            index = 0
            for atom in self.system.ProteinDict[chain_name].residue_dict[string.atoi(res_number)].atoms:
                if '%s'%(index) in sel_list:
                    items.append(atom)
                index = index + 1
        elif tokens[1] == 'nucleotide':
            index = 0
            for atom in self.system.NucleotideChainDict[chain_name].residue_dict[string.atoi(res_number)].atoms:
                if '%s'%(index) in sel_list:
                    items.append(atom)
                index = index + 1
        elif tokens[1] == 'ligand':
            index = 0
            for atom in self.system.LigandDict['%s%s'%(res_number, chain_name)].atoms:
                if '%s'%(index) in sel_list:
                    items.append(atom)
                index = index + 1
        return items
