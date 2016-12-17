
import os
import sys
import string
from Tkinter import *
from tkFileDialog import *
import profile
import MolecularSystem
from GUICards import *
from MolnirMolecularViewer import MolnirMolecularViewer
import Pmw
from Bio import SCOP
from ScopDomainViewer import *
from DatabaseDomainViewer import *

if __name__ == '__main__':
	molnir = Tk()
	geometry_string = "%dx%d%+d%+d" %(1000,22,1,1) # width,height,x-offset,y-offset
	molnir.geometry(geometry_string)


class MolnirGUI(Frame):
  def __init__(self, parent=None):
    Frame.__init__(self, parent)
    self.config()
    self.master.title('Molnir Metasuite')
    self.pack(expand=YES, fill=BOTH)
    self.makeWidgets()
  def makeWidgets(self):
    # start a menu
    self.menu = MenuSystem(self)
    self.system_windows = []
    self.scopViewer_windows = []
    self.databaseViewer_windows = []
class MenuSystem(Frame):
  def __init__(self, parent):
    Frame.__init__(self, parent)
    self.pack(side=TOP, fill=X, expand=NO, anchor=N)
    self.parent = parent
    # File operations
    file_button=Menubutton(self, text='File', underline=0)
    file_button.pack(side=LEFT, fill=X, expand=NO)
    file=Menu(file_button)
    file.add_command(label='Open a PDB', command=self.openNewSystem, underline=0)
    file_button.config(menu=file)
    c_lambda = lambda i=0:self.reloadSystem(i)
    file.add_command(label='Reopen a PDB', command=c_lambda, underline=0)
    file_button.config(menu=file)

    dbs_button = Menubutton(self, text='Databases', underline=0)
    dbs_button.pack(side=LEFT, fill=X, expand=NO)    
    databases=Menu(dbs_button, tearoff=0)     
#    scop = Menu(dbs_button, tearoff=0)
#    scop.add_command(label='Entire SCOP Structure', command=(lambda strg = "": self.openScopBrowser(strg)), underline=0)
#    scop.add_command(label='All SCOP Leaf Nodes', command=(lambda strg = "": self.openScopBrowser(strg)), underline=0)
#    databases.add_cascade(label='Browse SCOP', menu=scop, underline=0)   
    databases.add_command(label='SCOP Database', command=(lambda strg = "": self.openScopBrowser(strg)), underline=0)
    databases.add_command(label='User Defined Database', command=(lambda strg = "": self.openUserDbBrowser(strg)), underline=0)
    dbs_button.config(menu=databases)
  
  def openScopBrowser(self, strg):
   
    new_window = Toplevel()
    pw = Pmw.PanedWidget(new_window)	
    pw.add('top', min=10)
    pw.add('bottom', min = 10)
    
    geometry_string = "%dx%d%+d%+d" %(700,800,50+(25*len(self.parent.system_windows)),30+(25*len(self.parent.system_windows))) # width,height,x-offset,y-offset
    new_window.geometry(geometry_string)
    new_window.title("SCOP Domain Viewer")
    
    ScopViewer(pw, new_window, self).pack(side=TOP, expand=YES)
    pw.pack(expand=1, fill=BOTH) 
    
    bottomPane = pw.pane('bottom')
    
    bottomPane.system = MolecularSystem.System(self, bottomPane)
    bottomPane.system.load_pdb("")
    bottomPane.system.color_ribbon_by_chain()
  
    bottomPane.viewer = MolnirMolecularViewer(bottomPane, bottomPane.system)
    bottomPane.viewer.loadSystem(bottomPane.system)
    bottomPane.cardframe = CardFrame(bottomPane, bottomPane.system)
    bottomPane.cardframe.pack(expand=NO, fill=X)
    
    self.parent.scopViewer_windows.append(new_window)
   
   
  def openUserDbBrowser(self, strg):
    root_path = "C:\\CourseWork\\CS499\\MyTest\\SCOP_Pdbs\\1-10"  
    root_path = askdirectory(initialdir=root_path, title="select directory", mustexist=1)
    print root_path
    
    new_window = Toplevel()
    pw = Pmw.PanedWidget(new_window)	
    pw.add('top', min=10)
    pw.add('bottom', min = 10)
    
    geometry_string = "%dx%d%+d%+d" %(700,800,50+(25*len(self.parent.system_windows)),30+(25*len(self.parent.system_windows))) # width,height,x-offset,y-offset
    new_window.geometry(geometry_string)
    new_window.title("User Database Viewer")
  
    UserDbViewer(pw, root_path, new_window).pack(side=TOP, expand=YES)
    pw.pack(expand=1, fill=BOTH) 
    bottomPane = pw.pane('bottom')
    
    bottomPane.system = MolecularSystem.System(self, bottomPane)
    bottomPane.system.load_pdb("")
    bottomPane.system.color_ribbon_by_chain()
  
    bottomPane.viewer = MolnirMolecularViewer(bottomPane, bottomPane.system)
    bottomPane.viewer.loadSystem(bottomPane.system)
    bottomPane.cardframe = CardFrame(bottomPane, bottomPane.system)
    bottomPane.cardframe.pack(expand=NO, fill=X)

    
    self.parent.system_windows.append(new_window)
     
  def openNewSystem(self):
    strg = askopenfilename(title = 'open a new PDB', defaultextension='.pdb', filetypes=[("Protein Data Bank", "*.pdb"),("all files", "*")])
    # stores systems
    if len(strg)>0:
      # create a window
      new_window = Toplevel()
      geometry_string = "%dx%d%+d%+d" %(700,600,50+(25*len(self.parent.system_windows)),30+(25*len(self.parent.system_windows))) # width,height,x-offset,y-offset
      new_window.geometry(geometry_string)
      # create a system
      new_system = MolecularSystem.System(self, new_window)
      new_system.load_pdb(strg)
      new_system.color_ribbon_by_chain()
      new_window.title(new_system.header)
      # store it in the window
      new_window.system = new_system
      # open a molecular viewer
      new_viewer = MolnirMolecularViewer(new_window, new_window.system)
      new_window.viewer = new_viewer
      new_cardframe = CardFrame(new_window, new_window.system)
      new_cardframe.pack(expand=NO, fill=X)
      new_window.cardframe = new_cardframe 
      self.parent.system_windows.append(new_window)
  
  def reloadSystem(self, which_window):
    strg = askopenfilename(title = 'open a new PDB', defaultextension='.pdb', filetypes=[("Protein Data Bank", "*.pdb"),("all files", "*")])
    self.loadSystem(strg)
    
  def loadSystem(self, strg, obj=None):
	if len(strg)>0:
	#      self.closeSystem(which_window)
		if obj is None:
			obj = Toplevel()
					
		obj.string = strg 	#identifier
		obj.system = MolecularSystem.System(self, obj)
		obj.system.load_pdb(strg)
		obj.system.color_ribbon_by_chain()
		obj.title(obj.system.header)
		obj.viewer = MolnirMolecularViewer(obj, obj.system)
		obj.viewer.loadSystem(obj.system)
		obj.cardframe = CardFrame(obj, obj.system)
		obj.cardframe.pack(expand=NO, fill=X)
		self.parent.system_windows.append(obj)  
        
  def closeSystem(self, which_window):
    print "len before %d"%(len(self.parent.system_windows))
    self.parent.system_windows[which_window].destroy()
    if len(self.parent.system_windows) > which_window+1:
      self.parent.system_windows = self.parent.system_windows[:which_window] + self.parent.system_windows[which_window+1:]
    else:
      self.parent.system_windows = self.parent.system_windows[:which_window]
    print "len after %d"%(len(self.parent.system_windows))
    
# holds a set of GUI Cards
class CardFrame(Frame):
  def __init__(self, parent, my_system):
    self.parent = parent
    self.system = my_system
    ht = len(self.system.ProteinList + self.system.NucleotideChainList)
    Frame.__init__(self, parent, height=ht)
    if len(my_system.LigandList) > 0:
      ht = ht + 1
    if len(my_system.WaterList) > 0:
      ht = ht + 1
    Frame.__init__(self, parent, height=ht)
    self.proteinCardList = []
    self.nucleotideChainCardList = []
    self.CFSCardlist = []
    self.reloadSystem(self.system, 0)
  def reloadSystem(self, my_system, reload=1):
    if (reload):
      self.systemCard.pack_forget()
      self.systemCard.destroy()
      del self.systemCard
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
    self.system = my_system
    ht = len(self.system.ProteinList + self.system.NucleotideChainList)
    if len(my_system.LigandList) > 0:
      ht = ht + 1
    if len(my_system.WaterList) > 0:
      ht = ht + 1
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
    self.ligandCardList = LigandCard(self, self.system.LigandList)
    self.ligandCardList.pack(expand=YES, fill=X)
    self.CFSCardlist = []
def notdone():
  print('not yet implemented')

if __name__ == '__main__':
  molnir.ui = MolnirGUI()
  molnir.ui.mainloop()
