import copy

class SystemMemory:
    def __init__(self, system):
        # constructs self.data by calling _get_data on every component including the system
        self.system = system
        self.data = {'system':{}, 'prots':{}, 'nuccs':{}, 'ligs':{}, 'wats':{}}
        self.data['system'] = {'hbonds':{}}
        self.data['system'] = copy.deepcopy(self._get_data(system))
        for pchain in system.ProteinList:
            self.data['prots'][pchain.key] = {'atoms':{},'residues':{},'molecule':'None'}
            for atom in pchain.atoms:
                self.data['prots'][pchain.key]['atoms'][atom.atom_number] = copy.deepcopy(self._get_data(atom))
            for res in pchain.residues:
                self.data['prots'][pchain.key]['residues'][res.res_number] = copy.deepcopy(self._get_data(res))
            self.data['prots'][pchain.key]['molecule'] = copy.deepcopy(self._get_data(pchain))
        for nchain in system.NucleotideChainList:
            self.data['nuccs'][nchain.key] = {'atoms':{},'residues':{},'molecule':'None'}
            for atom in nchain.atoms:
                self.data['nuccs'][nchain.key]['atoms'][atom.atom_number] = copy.deepcopy(self._get_data(atom))
            for res in nchain.residues:
                self.data['nuccs'][nchain.key]['residues'][res.res_number] = copy.deepcopy(self._get_data(res))
            self.data['nuccs'][nchain.key]['molecule'] = copy.deepcopy(self._get_data(nchain))
        for lig in system.LigandList:
            self.data['ligs'][lig.key] = {'atoms':{},'molecule':'None'}
            for atom in lig.atoms:
                self.data['ligs'][lig.key]['atoms'][atom.atom_number] = copy.deepcopy(self._get_data(atom))
            self.data['ligs'][lig.key]['molecule'] = copy.deepcopy(self._get_data(lig))
        for wat in system.WaterList:
            self.data['wats'][wat.key] = {'atoms':{},'molecule':'None'}
            for atom in wat.atoms:
                self.data['wats'][wat.key]['atoms'][atom.atom_number] = copy.deepcopy(self._get_data(atom))
            self.data['wats'][wat.key]['molecule'] = copy.deepcopy(self._get_data(wat))
    def restore_system(self):
        self._set_data(self.system, copy.deepcopy(self.data['system']))
        for pchain in self.system.ProteinList:
            for atom in pchain.atoms:
                self._set_data(atom, copy.deepcopy(self.data['prots'][pchain.key]['atoms'][atom.atom_number]))
            for res in pchain.residues:
                self._set_data(res, copy.deepcopy(self.data['prots'][pchain.key]['residues'][res.res_number]))
            self._set_data(pchain, copy.deepcopy(self.data['prots'][pchain.key]['molecule']))
        for nchain in self.system.NucleotideChainList:
            for atom in nchain.atoms:
                self._set_data(atom, copy.deepcopy(self.data['nuccs'][nchain.key]['atoms'][atom.atom_number]))
            for res in nchain.residues:
                self._set_data(res, copy.deepcopy(self.data['nuccs'][nchain.key]['residues'][res.res_number]))
            self._set_data(nchain, copy.deepcopy(self.data['nuccs'][nchain.key]['molecule']))
        for lig in self.system.LigandList:
            for atom in lig.atoms:
                self._set_data(atom, copy.deepcopy(self.data['ligs'][lig.key]['atoms'][atom.atom_number]))
            self._set_data(lig, copy.deepcopy(self.data['ligs'][lig.key]['molecule']))
        for wat in self.system.WaterList:
            for atom in wat.atoms:
                self._set_data(atom, copy.deepcopy(self.data['wats'][wat.key]['atoms'][atom.atom_number]))
            self._set_data(wat, copy.deepcopy(self.data['wats'][wat.key]['molecule']))
    def _get_data(self, item):
        pass
    def _set_data(self, item, data):
        pass
        

class SelectionMemory(SystemMemory):
    def __init__(self, system):
        SystemMemory.__init__(self, system)
    def _get_data(self, item):
        return item.selected
    def _set_data(self, item, data):
        item.selected = data
    def is_molecule_selection_current(self, mol):
        key_dict = {'System':'system',
                    'MolecularComponents.classProtein':'prots',
                    'MolecularComponents.classNucleotideChain':'nuccs',
                    'MolecularComponents.classLigand':'ligs',
                    'MolecularComponents.classWater':'wats'}
        for atom in mol.atoms:
            if atom.selected != self.data[key_dict[mol.__module__]][mol.key]['atoms'][atom.atom_number]:
                return 0
        else:
            return 1

class ObjectSelectionVisitor(SystemMemory):
    def __init__(self, system):
        SystemMemory.__init__(self, system)
        self.restore_system()
    def _get_data(self, item):
        return 1
    def _set_data(self, item, data):
        item.selected = data

class AtomSelectionVisitor(SystemMemory):
    def __init__(self, system):
        SystemMemory.__init__(self, system)
        self.restore_system()
    def _get_data(self, item):
        if item.__module__ == 'MolecularComponents.classAtom':
            item.selected = 1
        else:
            item.selected = 0
    def _set_data(self, item, data):
        item.selected = data
    
class QuickRenderVisitor(SystemMemory):
    def __init__(self, system):
        SystemMemory.__init__(self, system)
        self.restore_system()
    def _get_data(self, item):
        return 1
    def _set_data(self, item, data):
        if item.__module__ == 'System':
            item.vtk_arg_list['hbonds']['selected'] = 1
            item.vtk_arg_list['hbonds']['representation'] = 'line'
        if item.__module__ in ['MolecularComponents.classProtein', 'MolecularComponents.classNucleotideChain']:
            item.vtk_arg_list['trace']['selected'] = 1
            item.vtk_arg_list['trace']['representation'] = 'line'
            item.vtk_arg_list['volume']['selected'] = 0
            item.vtk_arg_list['atoms']['selected'] = 0
        if item.__module__ == "MolecularComponents.classLigand":
            item.vtk_arg_list['volume']['selected'] = 0
            item.vtk_arg_list['atoms']['selected'] = 1
            item.vtk_arg_list['atoms']['representation'] = 'wireframe'
        if item.__module__ == "MolecularComponents.ClassWater":
            item.vtk_arg_list['volume']['selected'] = 0
            item.vtk_arg_list['atoms']['selected'] = 0
        if item.__module__ == "MolecularComponents.classAtom":
            if item.parent.__module__ in ['MolecularComponents.classAminoAcid', 'MolecularComponents.classNucleotide']:
                item.vtk_arg_list['atoms']['color'] = item.parent.parent.default_chain_color
                item.vtk_arg_list['atoms']['opacity'] = 1.0
            else:
                item.vtk_arg_list['atoms']['color'] = [0.8,0.8,0.8]
                item.vtk_arg_list['atoms']['opacity'] = 1.0
        if item.__module__ in ["MolecularComponents.classAminoAcid", "MolecularComponents.classNucleotide"]:
            item.vtk_arg_list['trace']['color'] = item.parent.default_chain_color

class GraphicsMemory(SystemMemory):
    def __init__(self, system):
        SystemMemory.__init__(self, system)
        self.key_dict = {'System':'system',
                         'MolecularComponents.classProtein':'prots',
                         'MolecularComponents.classNucleotideChain':'nuccs',
                         'MolecularComponents.classLigand':'ligs',
                         'MolecularComponents.classWater':'wats'}
    def _get_data(self, item):
        return item.vtk_arg_list
    def _set_data(self, item, data):
        item.vtk_arg_list = data
    def get_volume_args(self, mol):
        # specialized for Molecule surfaces
        return self.data[self.key_dict[mol.__module__]][mol.key]['molecule']
    def get_representation(self, obj, type):
        if type == 'atoms':
            self.get_atoms_representation(obj)
        elif type == 'trace':
            self.get_trace_representation(obj)
        elif type == 'volumes':
            self.get_volumes_representation(obj)
        elif type == 'hbonds':
            self.get_hbonds_representation(obj)
    def get_hbonds_representation(self, obj):
        return self.data[self.key_dict[obj.__module__]]['hbonds']['representation']
    def get_atoms_representation(self, obj):
        return self.data[self.key_dict[obj.__module__]][obj.key]['molecule']['atoms']['representation']
    def get_trace_representation(self, obj):
        return self.data[self.key_dict[obj.__module__]][obj.key]['molecule']['trace']['representation']
    def get_volumes_representation(self, obj):
        return self.data[self.key_dict[obj.__module__]][obj.key]['molecule']['volumes']['representation']

