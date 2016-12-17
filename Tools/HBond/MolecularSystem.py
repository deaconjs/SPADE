# python imports
import sys
import os
import string
import re
import parms
import pickle
import copy
# internal imports
sys.path.append(os.getcwd())
from MolecularComponents.classAtom            import Atom
from MolecularComponents.classMolecule        import Molecule
from MolecularComponents.classWater           import Water
from MolecularComponents.classLigand          import Ligand
from MolecularComponents.classAminoAcid       import AminoAcid
from MolecularComponents.classNucleotide      import Nucleotide
from MolecularComponents.classPolymer         import Polymer
from MolecularComponents.classProtein         import Protein
from MolecularComponents.classNucleotideChain import NucleotideChain

from MolecularComponents.MathFunc             import *

verbose = 0

class System:
    def __init__(self, parent):
        """ initialize a system object from a pdb file """
        self.parent = parent
        self.__module__ = "System"
        
    def load(self, filename, type):
        if type == 'bps':
            self.load_system(filename)
        elif type == 'pdb':
            self.load_pdb(filename)
        
    def load_pdb(self, filename):
        print "Opening pdb file %s\n"%(os.path.split(filename)[1])
        """ parse a pdb file and load it as a system object."""
        self.filename = filename
        self.selected = 1
        self.ProteinList         = []
        self.NucleotideChainList = []
        self.LigandList          = []
        self.WaterList           = []
        self.HeaderLines         = []
        self.HBonds              = []
        self.header             = ''
        # open the file
        if len(filename) == 0:
            return
        PDBfile = open(filename, 'r')
        # separate atoms, HETATMS, and other
        ATOMlines   = []
        HETATMlines = []
        ANISlines   = []        # anisotropy cards are ignored
        otherlines  = []
        while 1:
            line = PDBfile.readline()
            # should be able to tell from the first 4 characters of each line
            if line[0:4] == 'ATOM':
                ATOMlines.append(line)
            elif line[0:4] == 'HETA':
                HETATMlines.append(line)
            elif line[0:4] == 'ANIS':
                ANISlines.append(line)
            elif line[0:4] in ['TER ','MAST','CONE','END ']:
                pass
            else:
                self.HeaderLines.append(line)
            if not line: break
        for line in self.HeaderLines:
            if line[0:4] == 'HEAD':
                self.header = string.strip(line[6:])
                
        if verbose:
            print "%d atoms" % len(ATOMlines)
            print "%d hetero atoms" % len(HETATMlines)
            print "%d header lines" % len(otherlines)

        # next make a list of lines for the amino acids and nucleotides
        AATypes = parms.get('AATypes')
        NUCTypes = parms.get('NUCTypes')
        WATTypes = parms.get('WATTypes')
        AAlines = []
        NUClines = []
        WATlines = []
        OTHlines = []
        for line in ATOMlines:
            if line[17:20] in AATypes:
                AAlines.append(line)
            elif line[17:20] in NUCTypes:
                NUClines.append(line)
            elif line[17:20] in WATTypes:
                WATlines.append(line)
            else:
                OTHlines.append(line)
        for line in HETATMlines:
            if line[17:20] in WATTypes:
                WATlines.append(line)
            else:
                OTHlines.append(line)

        # now separate the chains
        # first collect the chain names for AAlines
        if len(AAlines) > 0:
            proteinchainlist = []
            for line in AAlines:
                if line[21:22] in proteinchainlist:
                    pass
                else:
                    proteinchainlist.append(line[21:22])
            # and create a list for each
            separatedlists = []
            for chainname in proteinchainlist:
                newlist = []
                for line in AAlines:
                    if chainname == line[21:22]:
                        newlist.append(line)
                separatedlists.append(newlist)
            # now create all of the Protein objects
            for slist in separatedlists:
                newprotein = Protein(slist, self)
                self.ProteinList.append(newprotein)

        # next make a list of all nucleotide lines
        # and create all of the NucleotideChain objects
        if len(NUClines) > 0:
            NUCchainlist = []
            for line in NUClines:
                if line[21:22]in NUCchainlist:
                    pass
                else:
                    NUCchainlist.append(line[21:22])
            if verbose:
                if len(NUCchainlist) == 1:
                    print "1 nucleotide chain loading"
                else:
                    print "%d nucleotide chains loading" % len(NUCchainlist)
            # and create a list for each
            separatedlists = []
            for chainname in NUCchainlist:
                newlist = []
                for line in NUClines:
                    if chainname == line[21:22]:
                        newlist.append(line)
                separatedlists.append(newlist)
            # now create all of the Nucleotide chain objects
            for slist in separatedlists:
                newnucchain = NucleotideChain(slist, self)
                self.NucleotideChainList.append(newnucchain)

        # now make a list of all other molecule types
        # and create all of the ligand Molecule objects
        # pass through the list, collecting by residue number
        if len(OTHlines) > 0:
            last_rez_num = None
            current_molecule = []
            for line in OTHlines:
                current_rez_num = string.atoi(line[23:26])
                if current_rez_num == last_rez_num:
                    current_molecule.append(line)
                else:
                    if len(current_molecule)>0:
                        new_molecule = Ligand(current_molecule, self)
                        self.LigandList.append(new_molecule)
                    current_molecule = []
                    current_molecule.append(line)
                last_rez_num = current_rez_num
            new_molecule = Ligand(current_molecule, self)
            self.LigandList.append(new_molecule)

            mol_names = []
            for mol in self.LigandList:
                mol_names.append(mol.res_type)
            if verbose:
                if len(self.LigandList) == 1:
                    print "1 ligand molecule %s" % (mol_names)
                else:
                    print "%d ligand molecules %s" % (len(self.LigandList), mol_names)

        # once ligands have been read in, then process the water list
        # so that they can be numbered sequentially
        if len(WATlines) > 0:
            last_rez_num = None
            current_molecule = []
            for line in WATlines:
                current_rez_num = line[23:26]
                if current_rez_num == last_rez_num:
                    current_molecule.append(line)
                else:
                    if len(current_molecule)>0:
                        new_molecule = Water(current_molecule,self)
                        self.WaterList.append(new_molecule)
                    current_molecule = []
                    current_molecule.append(line)
                last_rez_num = current_rez_num
            if len(current_molecule) > 0:
                new_molecule = Water(current_molecule,self)
                self.WaterList.append(new_molecule)
            if verbose:
                if len(self.WaterList) == 1:
                     print "1 water molecule"
                else:
                    print "%d water molecules" % (len(self.WaterList))
        PDBfile.close()
        self._supplementary_initialization()

    def _supplementary_initialization(self):
        """ Completes initialization by rebuilding the object from stored data """
        self.MoleculeList        = []
        self.SmallMoleculeList   = []
        self.PolymerList         = []
        self.PolymerDict         = {}
        self.ProteinDict         = {}
        self.NucleotideChainDict = {}
        self.LigandDict          = {}
        self.WaterDict           = {}
        # now create dictionaries out of the lists
        for pchain in self.ProteinList:
            self.ProteinDict[pchain.key] = pchain
        for nchain in self.NucleotideChainList:
            self.NucleotideChainDict[nchain.key] = nchain
        for lig in self.LigandList:
            self.LigandDict[lig.key] = lig
        for wat in self.WaterList:
            self.WaterDict[wat.key] = wat

        # these can be removed to add water functionality (much of which is not yet written,
        # but which will follow the same patterns as LigandList... maybe the WaterList should
        # be added to ligandList...
        #self.WaterList = []
        #self.WaterDict = {}

        # count the atoms first, for a printout and to set 
        items = [self.ProteinList, self.NucleotideChainList, self.LigandList, self.WaterList]
        count = 0
        for item in items:
            for mol in item:
                count = count + len(mol.atoms)
        self.atom_count = count
        # set up the MoleculeList
        items = [self.ProteinList, self.NucleotideChainList, self.LigandList]
        count = 0
        for item in items:
            for mol in item:
                self.MoleculeList.append(mol)
        # use MoleculeList to create atom_dict features for all molecules
        for mol in self.MoleculeList:
            mol.atom_dict = {}
            for atom in mol.atoms:
                mol.atom_dict[atom.atom_number] = atom
        # set up the SmallMoleculeList
        items = [self.LigandList, self.WaterList]
        count = 0
        for item in items:
            for mol in item:
                self.SmallMoleculeList.append(mol)
        # set up PolymerList and PolymerDict
        items = [self.ProteinList, self.NucleotideChainList]
        count = 0
        for item in items:
            for mol in item:
                self.PolymerList.append(mol)
                self.PolymerDict[mol.chain_name] = mol
        # create b_factor features
        self.normalize_b_factors()
        print 'loaded %s atoms in %d polymers, %d ligands, %d waters'%(self.atom_count, len(self.ProteinList)+len(self.NucleotideChainList), len(self.LigandList), len(self.WaterList))
                
    def select(self):
        """ select all of the objects in this molecular system """
        for pchain in self.ProteinList:
            pchain.select()
        for nchain in self.NucleotideChainList:
            nchain.select()
        for lig in self.LigandList:
            lig.select()
        for wat in self.WaterList:
            wat.select()
            
    def deselect(self):
        """ deselect all of the objects in this molecular system """
        for pchain in self.ProteinList:
            pchain.deselect()
        for nchain in self.NucleotideChainList:
            nchain.deselect()
        for lig in self.LigandList:
            lig.deselect()
        for wat in self.WaterList:
            wat.deselect()

               
    def save_pdb(self, filename="tmp.pdb"):
        """ write out the structure to the given file name, in pdb-style """
        if verbose:
            print "currently only printing the protein atom lines"
        pdb_file = open(filename, 'w')
        pdb_file.writelines(self.HeaderLines)
        for pchain in self.ProteinList:
            pdb_file.writelines(pchain.get_pdb_lines())
        for nchain in self.NucleotideChainList:
            pdb_file.writelines(nchain.get_pdb_lines())
        for lig in self.LigandList:
            pdb_file.writelines(lig.get_pdb_lines())
        for wat in self.WaterList:
            pdb_file.writelines(wat.get_pdb_lines())
        pdb_file.close()
                
    def save_system(self, filename="tmp.bps"):
        save_dict = {}
        save_dict['filename'] = self.filename
        save_dict['selected'] = self.selected
        save_dict['ProteinList'] = self.ProteinList
        save_dict['NucleotideChainList'] = self.NucleotideChainList
        save_dict['LigandList'] = self.LigandList
        save_dict['WaterList'] = self.WaterList
        save_dict['HeaderLines'] = self.HeaderLines
        save_dict['header'] = self.header

        # save the file    
        if len(filename) == 0:
            return
        sys_file = open(filename, 'w')
        pickle.dump(save_dict, sys_file, 0)
        
    def load_system(self, filename):
        print 'opening %s'%(os.path.split(filename)[1])
        new_file = open(filename, 'r')
        save_dict = pickle.load(new_file)
        self.filename = save_dict['filename']
        self.selected = save_dict['selected']
        self.ProteinList = save_dict['ProteinList']
        self.NucleotideChainList = save_dict['NucleotideChainList']
        self.LigandList = save_dict['LigandList']
        self.WaterList = save_dict['WaterList']
        self.HeaderLines = save_dict['HeaderLines']
        self.header = save_dict['header']
        self._supplementary_initialization()
        
    def get_filename_by_extension(self, ext, chain_name = 'None'):
        """ use the pdb's filename to produce a filename of type .ext  It figures out whether to add a .
        chain_name is a unique identifier to tag on to the filename, i.e. molecule.key"""
        af = self.filename
        index = af.rfind('.')
        slash = af.rfind(os.sep)
        if( index == -1 or index < slash ):
            index = len( af )
        if chain_name == 'None':
            if ext[0] == ".":
                df = af[0:index] + ext
            else:
                df = af[0:index] + "." + ext
        else:
            if ext[0] == ".":
                df = af[0:index] + chain_name + ext
            else:
                df = af[0:index] + chain_name + "." + ext
            
        return df

    def get_bounds(self):
        x = 0.0
        y = 0.0
        z = 0.0
        count = 0
        for items in [self.ProteinList, self.NucleotideChainList, self.LigandList]:
            for mol in items:
                for atom in mol.atoms:
                    x = x + atom.x
                    y = y + atom.y
                    z = z + atom.z
                    count = count + 1
                    
        center_x = x / float(count)
        center_y = y / float(count)
        center_z = z / float(count)
        
        min_x = max_x = center_x
        min_y = max_y = center_y
        min_z = max_z = center_z
        
        for items in [self.ProteinList, self.NucleotideChainList, self.LigandList]:
            for mol in items:
                for atom in mol.atoms:
                    if atom.x >= max_x:
                        max_x = atom.x
                    if atom.x < min_x:
                        min_x = atom.x
                    if atom.y >= max_y:
                        max_y = atom.y
                    if atom.y < min_y:
                        min_y = atom.y
                    if atom.z >= max_z:
                        max_z = atom.z
                    if atom.z < min_z:
                        min_z = atom.z
        return [min_x,max_x,min_y,max_y,min_z,max_z]
    
    def normalize_b_factors(self):
        max_b = 0.0
        min_b = 1000.0
        if self.MoleculeList[0].atoms[0].b_factor != None:
            for mol in self.MoleculeList:
                for atom in mol.atoms:
                    if atom.b_factor < min_b:
                        min_b = atom.b_factor
                    if atom.b_factor > max_b:
                        max_b = atom.b_factor
            print 'b_factor range %s %s'%(min_b, max_b)
            if max_b != 0.0:
                for mol in self.MoleculeList:
                    for atom in mol.atoms:
                        if atom.b_factor < (0.5*max_b):
                            atom.features['b_factor'] = 0.0
                        else:
                            atom.features['b_factor'] =(atom.b_factor-(0.5*max_b)) / (max_b-(0.5*max_b))
    def renumber_atoms(self):
        atom_number = 1
        for pchain in self.ProteinList:
            for atom in pchain.atoms:
                atom.atom_number = atom_number
                atom_number += 1
        for nchain in self.NucleotideChainList:
            for atom in nchain.atoms:
                atom.atom_number = atom_number
                atom_number += 1
        for lig in self.LigandList:
            for atom in lig.atoms:
                atom.atom_number = atom_number
                atom_number += 1
        for wat in self.WaterList:
            for atom in wat.atoms:
                atom.atom_number = atom_number
                atom_number += 1

    def create_hydrogen(self, chain_name, res_number, atom_name, coordinates_list, redo=False):
        atom_place = atom_name[1:]
        hyd_present = 1
        target_res = self.PolymerDict[chain_name].residues_dict[str(res_number)] 
        for atom in target_res.atoms:
            if atom.atom_type[0] == 'H':
                if atom_place in atom.atom_type:
                    hyd_present += 1
        atom_number = 0
        atom_type   = 'H' + atom_place + '%s'%(hyd_present)
        res_type    = target_res.res_type
        x,y,z = coordinates_list[0], coordinates_list[1], coordinates_list[2]
        new_atom = Atom(self.PolymerDict[chain_name].residues_dict[str(res_number)], "ATOM %6s%5s %3s %1s%4s %11.3f %7.3f %7.3f%26s\n"%(atom_number,atom_type,res_type,chain_name,res_number,x,y,z, " "))
        target_res.atoms.append(new_atom)
        target_res.atoms_dict[new_atom.atom_type] = new_atom
        self.PolymerDict[chain_name].atoms.append(new_atom)
        return new_atom

    def create_hbonds (self,availAcceptors=None,availDonors=None,strict=True):
        """
         Create all of the hydrogen bonds. Right now it only handles protein atoms.
         Sets member 'HBonds':
             {'accAtom', 'donorAtom', 'hydAtom', 'HBondInfo'}
         Parameters: strict=>False... ignore angle requirements
        """
        # Get all of the available acceptors and donors
        avAccs = []
        avDonors = []
        if (availAcceptors != None):
            avAccs = availAcceptors
        if (availDonors != None):
            avDonors = availDonors

        if (availAcceptors == None or availDonors == None):
            for protein in self.ProteinList:
                protein.protonate ("protonate.dat")
                if (availAcceptors == None):
                    avAccs[len(avAccs):]=protein.get_avail_acceptors ("acceptors.dat")
                if (availDonors == None):
                    avDonors[len(avDonors):]=protein.get_avail_donors ("donors.dat")

        # Initialize all of the hydrogen bond arrays to empty
        self.HBonds=[]
        for acceptor in avAccs:
            acceptor['accAtom'].Acc_HBonds=[]
        for donor in avDonors:
            donor['donorAtom'].Donor_HBonds=[]
    
        # Enumerate all of the acceptor/donor pairs that could hydrogen bond
        for acceptor in avAccs:
            for donor in avDonors:
                # Cannot hbond to yourself
                if (donor['donorAtom'] == acceptor['accAtom']):
                    continue
                for hydAtom in donor['hydAtoms']:
                    HBondInfo = self.measure_hbond (acceptor,donor,hydAtom)
                    if (self.is_valid_hbond (HBondInfo,strict)):
                        hbond = {'accAtom': acceptor['accAtom'], 'donorAtom': donor['donorAtom'],
                                 'hydAtom': hydAtom, 'HBondInfo': HBondInfo}
                        acceptor['accAtom'].Acc_HBonds.append(hbond)
                        donor['donorAtom'].Donor_HBonds.append(hbond)
                        self.HBonds.append (hbond)

    def measure_hbond (self,acc,donor,hydAtom):
        """
         Measures the possible hbond between the acceptor and donor.
        """
        donorAtom=donor['donorAtom']
        accAtom=acc['accAtom']
        AAAtom=acc['AAAtom']
        Dv=r_[donorAtom.x,donorAtom.y,donorAtom.z]
        Hv=r_[hydAtom.x,hydAtom.y,hydAtom.z]
        Av=r_[accAtom.x,accAtom.y,accAtom.z]
        # There is no AAAtom for water molecules
        if (AAAtom == None and accAtom.res_type in ("HOH")):
            AAv=None
        else:
            AAv=r_[AAAtom.x,AAAtom.y,AAAtom.z]
        HBondInfo={}
        HBondInfo['dist_D_A']= distance (Dv,Av)
        HBondInfo['dist_H_A'] = distance (Hv,Av)
        HBondInfo['angle_D_H_A'] = angle (Dv,Hv,Av)
        if (AAv != None):
            # print Dv[0],Dv[1],Dv[2],Av[0],Av[1],Av[2]
            HBondInfo['angle_D_A_AA'] = angle (Dv,Av,AAv)
            HBondInfo['angle_H_A_AA'] = angle (Hv,Av,AAv)
        else:
            HBondInfo['angle_D_A_AA'] = 'Water'
            HBondInfo['angle_H_A_AA'] = 'Water'
        return HBondInfo

    def is_valid_hbond (self,HBondInfo,strict=True):
        """
            If strict == False then the angle criteria is ignored
        """
        self.maxDist_D_A = 3.9
        self.maxDist_H_A = 2.5
        self.minDist_D_A = 2.9
        self.minDist_H_A = 1.4
        self.minAngle_D_H_A = 90.0
        self.minAngle_D_A_AA = 90.0
        self.minAngle_H_A_AA = 90.0 
        
        if (HBondInfo['dist_D_A'] > self.maxDist_D_A or HBondInfo['dist_D_A'] < self.minDist_D_A):
            return False
        if (HBondInfo['dist_H_A'] > self.maxDist_H_A or HBondInfo['dist_H_A'] < self.minDist_H_A):
            return False
        
        if (HBondInfo['angle_D_H_A'] < self.minAngle_D_H_A and strict == True):
            return False 
        if (HBondInfo['angle_D_A_AA'] < self.minAngle_D_A_AA and strict == True):
            return False
        if (HBondInfo['angle_H_A_AA'] < self.minAngle_H_A_AA and strict == True):
            return False
        return True

    def print_hbonds (self,output=None):
        f = sys.stdout
        if (output != None):
            f = open (output,"w")
      
        f.write("%6s|%8s|%8s|%10s|%8s|%8s|%11s|%12s|%12s\n"\
                %('Acc AA','Acc Atom','Donor AA','Donor Atom','Dist D-A',
                  'Dist H-A','Angle D-H-A','Angle D-A-AA','Angle H-A-AA'))
        for HBond in self.HBonds:
            accAtom=HBond['accAtom']
            donorAtom=HBond['donorAtom']
            HBondInfo=HBond['HBondInfo']
            f.write("%6s|%8s|%8s|%10s|%8.2f|%8.2f|%11.2f"\
                    %(accAtom.res_type,accAtom.atom_type,donorAtom.res_type,donorAtom.atom_type,
                      HBondInfo['dist_D_A'],HBondInfo['dist_H_A'],HBondInfo['angle_D_H_A']))
            if (HBondInfo['angle_D_A_AA'] == 'Water'):
                f.write("|%12s|%12s\n"\
                        %(HBondInfo['angle_D_A_AA'],HBondInfo['angle_H_A_AA']))
            else:
                f.write("|%12.2f|%12.2f\n"\
                        %(HBondInfo['angle_D_A_AA'],HBondInfo['angle_H_A_AA']))

class ModificationSystem(System):
    def __init__(self, system):
        self.parent = system.parent
        self.filename = system.filename
        self.selected = system.selected
        self.system = system

    def fill_proteolysis_fragments(self, plys_type):
        """ passes through each of the protein chains, making a collection of fragments for the given proteolysis type """
        self.frag_weight_list   = []
        self.frag_start_list    = []
        self.frag_end_list      = []
        self.frag_info_by_chain = []
        for pchain in self.system.ProteinList:
            pchain.fill_proteolysis_fragments(plys_type)
            self.frag_info_by_chain.append(pchain.fragment_list)
            
    def get_proteolysis_fragments_within(self, queryweight, within_weight):
        """ returns only those fragments within the given range of molecular weights """
        qw = string.atof(queryweight)
        # count how many fragments are within the given tolerance
        within_count = 0
        for fragments_by_chain in self.frag_info_by_chain:
            for frag_info in fragments_by_chain:
                if abs(frag_info[2]-qw) <= within_weight:
                    within_count = within_count + 1
        # then use get_nearest_proteolysis_fragments to get the sorted results
        print 'getting the top %d'%(within_count)
        return self.get_nearest_proteolysis_fragments(queryweight, within_count)
        
    # return a list of size return_count of the closest fragments
    def get_nearest_proteolysis_fragments(self, queryweight, return_count=1):
        """ returns the specified number of fragments closest to the query weight given """
        if return_count == 0:        # an empty list
            return []
        qw             = string.atof(queryweight)
        if verbose:
            print "looking for %s"%(queryweight)
        frag_index     = 0
        nearest_index  = 0
        furthest_saved_dist = 1000000000.0
        saved_fragments   = []
        nearest_a      = 0
        nearest_b      = 0
        nearest_c      = 0
        return_list    = []
        chain_index = 0
        for fragments_by_chain in self.frag_info_by_chain:
            for frag_info in fragments_by_chain:
                frag_info.append(chain_index)
                if abs(frag_info[2]-qw) < furthest_saved_dist:
                    if len(saved_fragments) == 0:
                        if verbose:
                            print "inserting %s 1"%(frag_info)
                        saved_fragments.append(frag_info)
                    else:                   # sorted insert
                        saved_index = 0
                        for fragment in saved_fragments:
                            # if the one being inserted is closer than this saved one
                            if abs(frag_info[0]-qw) < abs(fragment[2]-qw):
                                # insert before the more distant one
                                if verbose:
                                    print "inserting %s 2"%(frag_info)
                                saved_fragments[saved_index:saved_index] = [frag_info]
                                break
                            saved_index = saved_index + 1
                        else:       # still want to insert it, at the end
                            if verbose:
                                print "inserting %s 3"%(frag_info)
                            saved_fragments.append(frag_info)
                        if len(saved_fragments) == return_count + 1: # if just inserted one
                            # if just inserted, update furthest_saved and delete the last one
                            furthest_saved_dist = abs(saved_fragments[return_count-1][2]-qw)
                            if verbose:
                                print "removing %s"%(saved_fragments[len(saved_fragments)-1])
                            saved_fragments.pop()
            chain_index = chain_index + 1
        if verbose:
            print "returning %d fragments:"%(len(saved_fragments))
        return saved_fragments
                       
    def color_blue(self):
        """ just color all of the residues blue """
        for pchain in self.system.ProteinList:
            for rez in pchain.residues:
                rez.vtk_arg_list['trace']['color'] = [0.1,0.1,1.0]

