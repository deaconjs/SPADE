# Python imports
import pickle
import os
import sys
import parms
import copy
import string
import math
# dependency imports
from Tkinter import *
import tkFileDialog
import Pmw
# module imports
sys.path.append('./Applications/AlignmentEditor')
import TreeSystem
sys.path.append('./Applications/Dynamics')
import DynamicsViewer
sys.path.append('./Tools/Selection')
import SystemSelectionDialog

# holds a set of GUI Cards
class CardFrame(Frame):
    def __init__(self, parent, my_system, master_window, standard_cards=1):
        """ initialize the information (card) holder of the system window """
        self.parent = parent
        self.system = my_system
        self.standard_cards = standard_cards
        self.master_window = master_window      # any transient windows should be transient to this
        Frame.__init__(self, parent, height=3)
        self.proteinCardList = []
        self.nucleotideChainCardList = []
        self.CFSCardlist = []
        if self.system != None:
            self.reloadSystem(self.system, 0)

    def removeCards(self):
        try:
            self.systemCard
        except (AttributeError):
            return                              # no system is present
        else:
            self.systemCard.pack_forget()
            self.systemCard.destroy()
            del self.systemCard

        try:
            self.ligandCardList
        except (AttributeError):
            pass
        else:
            self.ligandCardList.pack_forget()
            self.ligandCardList.destroy()
            del self.ligandCardList
        
        for card in self.proteinCardList:
            card.pack_forget()
            card.destroy()
            del card
        for card in self.nucleotideChainCardList:
            card.pack_forget()
            card.destroy()
            del card
        for card in self.CFSCardlist:
            card.pack_forget()
            card.destroy()
            del card
        
    def reloadSystem(self, my_system, reload=1):    # reload should be 0 if initializing
        """ destroy any old cards (if reload) and open GUIcards for a new system"""
        if (reload):
            self.removeCards()
        self.system = my_system
        if self.standard_cards:
            self.systemCard = SystemCard(self, self.system)
            self.systemCard.pack(expand=YES, fill=X)
            self.proteinCardList = []
            for p in self.system.ProteinList:
                self.proteinCardList.append(ProteinCard(self, p))
                self.proteinCardList[len(self.proteinCardList)-1].pack(expand=YES, fill=X)
            self.nucleotideChainCardList = []
            for n in self.system.NucleotideChainList:
                self.nucleotideChainCardList.append(NucleotideChainCard(self, n))
                self.nucleotideChainCardList[len(self.nucleotideChainCardList)-1].pack(expand=YES, fill=X)
            if len(self.system.LigandList) > 0:
                self.ligandCardList = LigandCard(self, self.system.LigandList)
                self.ligandCardList.pack(expand=YES, fill=X)
        self.CFSCardlist = []

class Card(Frame):
    def __init__(self, parent, select_command):
        Frame.__init__(self, parent, borderwidth=0, height=1)
        self.parent = parent
        if select_command != None:
            # make a button for selection
            self.select_button = Button(self, bitmap='gray12', width=14, height=14)
            self.select_button.config(command=select_command)
            self.select_button.pack(side=LEFT, anchor=W, expand=NO)
    def tearoff(self):
        m = self.item.__module__
        pchain_buffer = []       # holds pchains,nchains,ligands,waters
        if m=="System":
            for pchain_index in range(len(self.item.ProteinList)):
                if self.item.ProteinList[pchain_index].selected > 0:
                    pchain_buffer.append(self.item.ProteinList[pchain_index])
                    self.item.ProteinList = self.item.ProteinList[:pchain_index] + self.item.ProteinList[pchain_index+1:]
    def remove_tags(self):
        for tagn in self.block_holder.tag_names():
            self.block_holder.tag_delete(tagn)
    def select(self):
        # deselect the card's item
        if self.select_button.cget('bitmap') == 'gray12':
            self.select_button.config(bitmap='gray75')
            m = self.item.__module__
            if m=="MolecularComponents.classProtein" or m=="MolecularComponents.classNucleotideChain" or m=="System":
                if self.item.selected == 1:
                    self.item.deselect()
            else:
                for i in self.item:
                    if i.selected == 1:
                        i.deselect()
        # select the whole chain
        elif self.select_button.cget('bitmap') == 'gray75':
            self.select_button.config(bitmap='gray12')
            m = self.item.__module__
            if m == "MolecularComponents.classProtein" or m=="MolecularComponents.classNucleotideChain" or m=="System":
                self.item.select()
            else:
                for i in self.item:
                    i.select()
        self.parent.parent.parent.viewer.update_view()
        
class SystemCard(Card):
    def __init__(self, parent, system):
        Card.__init__(self, parent, self.select)
        self.item = system
        self.parent = parent
        self.info_button   = Menubutton(self, bitmap='info', direction=RIGHT,
                                        relief=RAISED,width=17, height=17)
        self.info_button.pack(side=LEFT, anchor=W, expand=NO)
        self.info_button.menu = Menu(self.info_button, tearoff=0)
        self.info_button["menu"] = self.info_button.menu
        self.info_button.menu.add_command(label="save system", command=self.save_system)
        self.info_button.menu.add_command(label="show header", command=self.view_header)
        # make a button for selection
        self.selection_button = Button(self, text='SELECTION', height=1)
        self.selection_button.config(command=self.launch_selection_dialog)
        self.selection_button.pack(side=LEFT,expand=NO)
        # make a button for tearing off
        self.tear_off_button = Button(self, bitmap='gray12', width=14, height=14)
        self.tear_off_button.config(command=self.tearoff)
        self.tear_off_button.pack(side=RIGHT, anchor=E, expand=NO)
        # create a label object as the title
        name_string = 'System'
        self.title = Label(self, height=1, width=6, text=name_string)
        self.title.pack(side=LEFT, anchor=W, expand=NO)
        
    def view_MD(self):
        reload(DynamicsViewer)
        self.MD_viewer = DynamicsViewer.MDViewer(self.item)
    
    def launch_selection_dialog(self):
        reload(SystemSelectionDialog)
        self.selection_top = Toplevel()
        geometry_string = "%dx%d%+d%+d" %(440,300,1,220) # width,height,x-offset,y-offset
        self.selection_top.geometry(geometry_string)
        self.selection_top.wm_transient(self.parent.master_window) # transient to the System window
        self.selection_top.title('Select Molecules, Residues, and Atoms')
        self.selection_dialog = SystemSelectionDialog.SystemSelectionDialog(self.parent.master_window.viewer, self.selection_top, self.item)
        self.selection_dialog.pack(expand='yes', fill='both')

    def view_header(self):
        win_title = 'Header for %s'%(self.item.filename)
        dialog = Pmw.TextDialog(self.parent, scrolledtext_labelpos = 'n', title = win_title, defaultbutton = 0, master=self.parent.master_window)
        hdr = ""
        for line in self.item.HeaderLines:
          hdr = hdr + string.strip(line) + '\n'
        dialog.insert('end', hdr)
    
    def save_system(self):
        system_save_name   = tkFileDialog.asksaveasfilename(title = 'Save System', defaultextension='.sps', filetypes=[("SPADE Pickled System", "*.sps"),("all files", "*")])
        self.item.save_system(system_save_name)

class backCard(Frame):
  def __init__(self, parent, frame_to_remove):
    Frame.__init__(self, parent, borderwidth=0, height=1)
    self.frame_to_remove = frame_to_remove
    self.parent = parent
    # make a button to recreate the system cards
    self.back_button = Button(self, text='System View', width=11, height=1)
    self.back_button.config(command=self.recreate_system_cards)
    self.back_button.pack(side=LEFT, anchor=W, expand=NO)
    
  def recreate_system_cards(self):
    self.frame_to_remove.pack_forget()
    self.frame_to_remove.destroy()
    del self.frame_to_remove
    self.pack_forget()
    self.destroy()
    self.parent.reloadSystem(self.parent.system, 0)
      
class ProteinCard(Card):
  def __init__(self, parent, chain):
    Card.__init__(self, parent, self.select)
    self.chain = chain
    # make a button for showing additional information
    self.info_button   = Menubutton(self, bitmap='info', direction=RIGHT,
                                    relief=RAISED,width=17, height=17)
    self.info_button.pack(side=LEFT, anchor=W, expand=NO)
    self.info_button.menu = Menu(self.info_button, tearoff=0)
    self.info_button["menu"] = self.info_button.menu
    self.info_button.menu.add_command(label="default colors", command=self.default_colors)
    self.info_button.menu.add_command(label="color by type", command=self.color_by_type)

    # make a button for tearing off
    self.tear_off_button = Button(self, bitmap='gray12', width=14, height=14)
    self.tear_off_button.config(command=self.tearoff)
    self.tear_off_button.pack(side=RIGHT, anchor=E, expand=NO)
    self.item = chain
    # create a label object as the title
    name_string = 'Prot ' + self.item.chain_name
    col3     = self.item.get_chain_color()
    bg_color = "#%02x%02x%02x"%(int(col3[0]*256),int(col3[1]*256),int(col3[2]*256))
    self.title = Label(self, height=1, width=6, text=name_string, bg=bg_color)
    self.title.pack(side=LEFT, anchor=W, expand=NO)
    # add a scrollbar
    self.xsbar = Scrollbar(self, orient='horizontal')
    self.xsbar.pack(side=LEFT)
    # create a text object for the sequence
    self.block_holder = Text(self, height=1, wrap=NONE)
    self.block_holder.insert('1.0', self.item.get_sequence())
    self.block_holder.pack(expand=YES, fill=X, anchor=N)
    # associate the scrollbar with the textbox
    self.xsbar.config(command=self.block_holder.xview)
    self.block_holder.config(xscrollcommand=self.xsbar.set)
    self.block_holder.config(state=DISABLED)
    
  def default_colors(self):
    self.remove_tags()
    self.parent.system.color_ribbon_by_chain()
    for rez in self.item.residues:
      for atom in rez.atoms:
        atom.vtk_arg_list[0][1] = atom.get_color('cpk')
    self.item.update_atom_actors()
    for rez in self.item.residues:
      rez.update_actors()
    self.parent.parent.parent.viewer.update()
    
  def color_by_type(self):
    buried = ['A','C','F','I','L','M','V','W']
    uncharged = ['G','N','P','Q','S','T','Y']
    chargedminus = ['D','E']
    chargedplus  = ['H','K','R']
    colors = ['cyan','orange','red','yellow']
    # first color the letters in the card
    cntr = 0
    for c in self.block_holder.get('1.0', END+'-1c'):
      str = "1.%d"%cntr
      if c in buried:
        self.block_holder.tag_add('buried', str)
      elif c in uncharged:
        self.block_holder.tag_add('uncharged', str)
      elif c in chargedminus:
        self.block_holder.tag_add('chargedminus', str)
      elif c in chargedplus:
        self.block_holder.tag_add('chargedplus', str)
      cntr = cntr + 1
    self.block_holder.tag_config('buried', background=colors[0])
    self.block_holder.tag_config('uncharged', background=colors[1])
    self.block_holder.tag_config('chargedminus', background=colors[2])
    self.block_holder.tag_config('chargedplus', background=colors[3])
    # now change the colors of all of the residues
    for rez in self.item.residues:
      c = rez.res_type1
      if c in buried:
        rez.vtk_arg_list[2][2] = [0.5,0.75,1.0]
      elif c in uncharged:
        rez.vtk_arg_list[2][2] = [1.0,0.5,0.0]
      elif c in chargedminus:
        rez.vtk_arg_list[2][2] = [0.9,0.2,0.2]
      elif c in chargedplus:
        rez.vtk_arg_list[2][2] = [1.0,1.0,0.2]
    # and recolor the atoms by type
    for rez in self.item.residues:
      for atom in rez.atoms:
        atom.vtk_arg_list[0][1] = atom.get_color('typ')
      rez.update_actors()
    self.item.update_atom_actors()
    self.parent.parent.parent.viewer.screen.Render()

class NucleotideChainCard(Card):
  def __init__(self, parent, chain):
    Card.__init__(self, parent, self.select)
    # make a button for showing additional information
    self.info_button   = Menubutton(self, bitmap='info', direction=RIGHT,
                                    relief=RAISED,width=17, height=17)
    self.info_button.pack(side=LEFT, anchor=W, expand=NO)
    self.info_button.menu = Menu(self.info_button, tearoff=0)
    self.info_button["menu"] = self.info_button.menu
    self.info_button.menu.add_command(label="default colors", command=self.remove_tags)
    # make a button for tearing off
    self.tear_off_button = Button(self, bitmap='gray12', width=14, height=14)
    self.tear_off_button.config(command=self.tearoff)
    self.tear_off_button.pack(side=RIGHT, anchor=E, expand=NO)
    self.item = chain
    # create a label object as the title
    name_string = 'Nuc ' + self.item.chain_name + ' '
    self.title = Label(self, height=1, width=6, text=name_string)
    self.title.pack(side=LEFT, anchor=W, expand=NO)
    # add a scrollbar
    self.xsbar = Scrollbar(self, orient='horizontal')
    self.xsbar.pack(side=LEFT)
    # create a text object for the sequence
    self.block_holder = Text(self, height=1, wrap=NONE)
    self.block_holder.insert('1.0', self.item.get_sequence())
    self.block_holder.pack(expand=YES, fill=X, anchor=N)
    # associate the scrollbar with the textbox
    self.xsbar.config(command=self.block_holder.xview)
    self.block_holder.config(xscrollcommand=self.xsbar.set)
    self.block_holder.config(state=DISABLED)

class LigandCard(Card):
  def __init__(self, parent, ligandList):
    Card.__init__(self, parent, self.select_all)
    # make a button for showing additional information
    self.info_button   = Menubutton(self, bitmap='info', direction=RIGHT,
                                    relief=RAISED,width=17, height=17)
    self.info_button.pack(side=LEFT, anchor=W, expand=NO)
    self.info_button.menu = Menu(self.info_button, tearoff=0)
    self.info_button["menu"] = self.info_button.menu
    self.info_button.menu.add_command(label="default colors", command=self.remove_tags)
    # make a button for tearing off
    self.tear_off_button = Button(self, bitmap='gray12', width=14, height=14)
    self.tear_off_button.config(command=self.tearoff)
    self.tear_off_button.pack(side=RIGHT, anchor=E, expand=NO)
    self.item = ligandList
    # create a label object as the title
    name_string = 'Ligand'
    self.title = Label(self, height=1, width=6, text=name_string)
    self.title.pack(side=LEFT, anchor=W, expand=NO)
    # add a scrollbar
    self.xsbar = Scrollbar(self, orient='horizontal')
    self.xsbar.pack(side=LEFT)
    # create a text object for the sequence
    self.block_holder = Text(self, height=1, wrap=NONE)
    self.block_holder.pack(expand=NO, fill=X, anchor=N)
    self.ligbuttonlist = []
    lig_index = 0
    for lig in self.item:
      # first enter the three letter name for the molecule
      self.block_holder.insert(END,lig.res_type)
      c_lambda = lambda lig_i=lig_index: self.select(lig_i)
      self.ligbuttonlist.append(Button(self.block_holder, bitmap='gray12', width=8, height=8, command=c_lambda))
      self.ligbuttonlist[len(self.ligbuttonlist)-1].pack(expand=NO,side=TOP,anchor=E)
      self.block_holder.window_create(END, window=self.ligbuttonlist[len(self.ligbuttonlist)-1])
      lig_index = lig_index + 1
    # associate the scrollbar with the textbox
    self.xsbar.config(command=self.block_holder.xview)
    self.block_holder.config(xscrollcommand=self.xsbar.set)
    self.block_holder.config(state=DISABLED)
    
  def select_all(self):
    if self.select_button.cget('bitmap') == 'gray12':
      self.select_button.config(bitmap='gray75')
      for i in range(0,len(self.item)):
        if self.ligbuttonlist[i].cget('bitmap') == 'gray12':
          self.ligbuttonlist[i].config(bitmap='gray75')
          if self.item[i].selected == 1:
            self.item[i].deselect()
    elif self.select_button.cget('bitmap') == 'gray75':
      self.select_button.config(bitmap='gray12')
      for i in range(0,len(self.item)):
        if self.ligbuttonlist[i].cget('bitmap') == 'gray75':
          self.ligbuttonlist[i].config(bitmap='gray12')
          if self.item[i].selected == 0:
            self.item[i].select()
    self.parent.parent.parent.viewer.update_view()
    
  def select(self, lig_num):
    # deselect the card's item
    if self.ligbuttonlist[lig_num].cget('bitmap') == 'gray12':
      self.ligbuttonlist[lig_num].config(bitmap='gray75')
      if self.item[lig_num].selected == 1:
        self.item[lig_num].deselect()
    # select the whole chain
    elif self.ligbuttonlist[lig_num].cget('bitmap') == 'gray75':
      self.ligbuttonlist[lig_num].config(bitmap='gray12')
      if self.item[lig_num].selected == 0:
        self.item[lig_num].select()
    self.parent.parent.parent.viewer.update_view()
    
