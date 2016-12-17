from Tkinter import *
import urllib
import re
import parms
import Pmw
import string
import SPADE

import MolecularSystem
from GUICards import *
import MolecularViewer
from os import *

#location of image: plusImage, minusImage
pdb_dir = parms.get('pdb_dir')
gifs_dir = path.abspath('./Tools/DatabaseManager') + path.normpath('/')

class DirNode(Frame):
    """    This class represents each Directory/File item in the treeWinText textbox.
           Creates a frame with:
            1.path to the item (path)
            2.actual item (name) 
            3.Pointer to the actual Node in SCOP Data structure
    """
    def __init__(self, pth, obj, parent=None):
        Frame.__init__(self, parent, bg='white')
        self.pack()
        self.toplevel= parent
        self.thisNode = obj
        self.thisPath = pth
        self.isDisplayed = 1
        self.isExpanded = 0
        self.isALeaf = 0
        self.pdbFiles = []
        self.createButton()
        self.createLabel()
            
    def createButton(self):
        self.plusImage = PhotoImage(file=os.path.join(gifs_dir,"plusnode.gif"))
        self.minusImage = PhotoImage(file=os.path.join(gifs_dir,"minusnode.gif"))
        self.pdbImage = PhotoImage(file=os.path.join(gifs_dir,"pdb.gif"))
        # determine the padding based on location of the item in the tree structure
        regex = re.compile('\\\\')
        padding = len(regex.split(self.thisPath))- len(regex.split(self.toplevel.root_path))
        self.emptyLabel = Label(self, text='', bg='white')
        self.emptyLabel.pack(side=LEFT, padx=((padding * 3) + 1))
        if path.isdir(os.path.join(self.thisPath,self.thisNode)):
            # if its a directory, add a button
            self.button = Button(self, image=self.plusImage, command=(lambda node=self: self.toplevel.selectNode(node)))			
        else:
            # if not, add a button then mark it as a leaf
            self.button = Button(self, image=self.minusImage, command=(lambda node=self: self.toplevel.selectNode(node)))
            self.isALeaf = 1
        # add a different graphic if it's a leaf
        if self.isCurrentNodeALeaf() == 1:
            self.button.config(image=self.pdbImage)

        self.button.pack(side=LEFT, padx=1)
        self.button.bind('<Button-1>', self.singleClick)

    def createLabel(self):
        self.label = Label(self, text='  '+self.thisNode+'  ', bg='white')
        self.label.pack(side=RIGHT)
        self.label.bind('<Double-1>', self.doubleClick)
        self.label.bind('<Button-1>', self.singleClick)

    def singleClick(self, event):
        self.toplevel.toggleSelection(self)
        self.toplevel.viewButton.config(state=DISABLED)
        #Determine if the node is a leaf, set self.isALeaf accordingly
        self.isCurrentNodeALeaf()
        if self.isALeaf == 1:
            self.toplevel.viewButton.config(state=NORMAL)
        self.displayLineage()
    
    def isCurrentNodeALeaf(self):
        """ Conditions for Leaf: (1) has no directories, and (2) has atleast 1 pdb file  """
        if path.isdir(os.path.join(self.thisPath,self.thisNode)):
            filesInSelected = listdir(os.path.join(self.thisPath,self.thisNode))
            self.isALeaf = 1
            for eachFile in filesInSelected:
                if path.isdir(os.path.join(self.thisPath,self.thisNode,eachFile)):
                    self.isALeaf = 0
                    return 0
                extension = string.split(eachFile, '.')[-1]
                if (string.upper(extension) in ['PDB', 'ENT']) and self.pdbFiles.count(eachFile) == 0:
                    self.pdbFiles.append(eachFile)
                #if (eachFile[len(eachFile)-4:].find('.pdb') == 0 or
                #    eachFile[len(eachFile)-4:].find('.ent') == 0) and
                #   self.pdbFiles.count(eachFile) == 0:
                #        self.pdbFiles.append(eachFile)
            if self.isALeaf == 1 and len(filesInSelected) > 0 and len(self.pdbFiles) > 0:
                return 1
        return 0

    def isCurrentNodeValid(self):
        """ if its a node, it should have a directory """
        filesInCurrent = listdir(os.path.join(self.thisPath,self.thisNode))
        for eachFile in filesInCurrent:
            if path.isdir(os.path.join(self.thisPath,self.thisNode,eachFile)):
                return 1
        return 0
                    
    def displayLineage(self):
        """The path to the selected Node and its contents (if it is a leaf) are displayed"""
        self.toplevel.lineageBox.config(state=NORMAL)
        self.toplevel.lineageBox.delete(1.0, END)
        if self.isALeaf == 1:
            # if its a leaf, add all of its parents
            filesInThisNode = listdir(os.path.join(self.thisPath,self.thisNode))
    
            for eachFile in filesInThisNode:
                self.toplevel.lineageBox.insert('1.0', "  "+eachFile+'\n')

            indx = '1.'+str(len("Files assoc. with selected Node: "))				
            self.toplevel.lineageBox.insert('1.0', "Files assoc. with selected Node: \n")
            self.toplevel.lineageBox.tag_add('label', '1.0', indx)
            self.toplevel.lineageBox.tag_config('label', font=('veranda', 8, 'bold'))				

        indx = '1.'+str(5)
        self.toplevel.lineageBox.insert('1.0', "Path: "+self.thisPath+self.thisNode+'\n\n')
        self.toplevel.lineageBox.tag_add('label', '1.0', indx)
        self.toplevel.lineageBox.tag_config('label', font=('veranda', 8, 'bold'))
                    
            
        self.toplevel.lineageBox.config(state=DISABLED)
            
    def doubleClick(self, event):
        self.toplevel.selectNode(self)
                    
    def expand(self):
        self.button.config(image=self.minusImage)
        self.isExpanded = 1
    
    def contract(self):
        self.button.config(image=self.plusImage)
        self.isExpanded = 0
                
                                
class UserDbViewer(Frame):
    """ Arguments:
                    parent  - the paned window object
                    pth     - the user-selected path
                    toplevel- the window in which Viewer is built
    """
    def __init__(self, parent, pth, toplevel):
        Frame.__init__(self, parent.pane('top'), bg='white')
        self.parent = parent
        self.focus()		
        self.spade_level = toplevel #parent.pane('bottom')
        self.root_path = pth
        self.toplevelwindow = toplevel
        self.pdbPath = ""
        self.treeWinText = Pmw.ScrolledText(self, labelpos = 'n', usehullsize = 1, label_text='Local Database Domain Nodes',
                        hull_width = 100, hull_height = 15, text_padx = 1, text_pady = 1, text_wrap='none', text_cursor='arrow')
        self.treeWinText.pack(side=LEFT, expand=YES, fill=BOTH)
        self.pack(expand=YES, fill=BOTH)
        self.viewButton = Button(self, text="View", state=DISABLED, height=1, command=self.displaySelected, width=5)
        self.viewButton.pack(side=LEFT, anchor=E, padx=10)
        self.lineageBox = Text(self, relief=FLAT, cursor = "arrow")
        self.lineageBox.pack(side=TOP, anchor=N)
        self.lineageBox.config(height=12, width=38, state=DISABLED)
        #Build toplevel node structure
        root_list = self.getDirectoryListing(self.root_path)
        self.lines = 1
        self.nodesList = []
        if len(root_list) > 0:
            for item in root_list:
                newDirNode = DirNode(self.root_path, item, self)
                if path.isdir(os.path.join(self.root_path,item)) and (newDirNode.isCurrentNodeValid() == 1 or newDirNode.isCurrentNodeALeaf() == 1):
                    self.nodesList.append(newDirNode)
                    indx = "%0.1f"%(self.lines)
                    self.treeWinText.window_create(indx, window=self.nodesList[-1])
                    indx = "%0.1f"%(self.lines+0.1)
                    self.treeWinText.insert(indx, '\n')
                    self.lines = self.lines + 1
                else:
                    newDirNode.destroy()
        
        self.currentSelection = self.nodesList[0]
        
        self.treeWinText.pack()
        self.treeWinText.configure(text_state='disabled')
        self.pack()		

    def getDirectoryListing(self, dest_path):
        """ Return list of files and directories in the path """
        objs = listdir(dest_path)
        directories = []
        files = []
        
        for item in objs:
            if path.isdir(os.path.join(dest_path,item)):
                directories.append(item)
            else:
                if string.upper(string.split(item, '.')[-1]) in ['PDB', 'ENT']:
                    files.append(item)

        for item in files:
            directories.append(item)
                
        return directories
    
    def selectNode(self, nde):
        """ This function implements the functionality when the +/- buttons are clicked """
        self.toggleSelection(nde)
        self.treeWinText.configure(text_state='normal')
        if nde.isExpanded == 1:
            pos = 0
            while nde.thisNode != self.nodesList[pos].thisNode:
                pos = pos + 1

            self.nodesList[pos].contract()
            
            pos = 0
            while len(self.nodesList) > pos:
                selectedNode = nde.thisPath+nde.thisNode
                iterNode = self.nodesList[pos].thisPath+self.nodesList[pos].thisNode
                if iterNode.find(selectedNode) == 0 and selectedNode != iterNode:
                    print pos
                    self.nodesList.pop(pos)
                    indx = "%0.1f"%(pos+1)
                    indx1 = "%0.1f"%(pos+2)
                    self.treeWinText.delete(indx, indx1)
                    pos = pos - 1
                pos = pos + 1
                
        else:
            pos = 0
            len_nodesList = len(self.nodesList)
            while pos < len_nodesList:
                item = self.nodesList[pos]
                if item.thisNode == nde.thisNode and item.thisPath == nde.thisPath:
                    chld_pos = 1
                    self.nodesList[pos].expand()
                    selectedNode = os.path.join(self.nodesList[pos].thisPath,self.nodesList[pos].thisNode)
                    if path.isdir(selectedNode):
                        subdirs = self.getDirectoryListing(selectedNode)
                        for child in subdirs:
                            newDirNode = DirNode(selectedNode, child, self)
                            if path.isdir(os.path.join(selectedNode,child)) and (newDirNode.isCurrentNodeValid() == 1 or newDirNode.isCurrentNodeALeaf() == 1):
                                self.nodesList.insert(pos+chld_pos, newDirNode)
                                indx = "%0.1f"%(pos+chld_pos+1)
                                self.treeWinText.window_create(indx, window=self.nodesList[pos+chld_pos])
                                indx = "%0.1f"%(pos+chld_pos+1+0.1)
                                self.treeWinText.insert(indx, '\n')
                                chld_pos = chld_pos + 1
                            else:
                                newDirNode.destroy()
                        pos = len_nodesList
                pos = pos + 1
                        
        self.treeWinText.configure(text_state='disabled')

    def toggleSelection(self, nde):
        self.currentSelection.label.config(bg='white', fg='black')
        nde.label.config(bg='blue', fg='white')
        self.currentSelection = nde	
    
    def displaySelected(self):
        """ If the selected node has more than 1 pdb files, the following function prompts the user to select one """
        if len(self.currentSelection.pdbFiles) > 1:
            self.pdbSelectdialog = Pmw.SelectionDialog(self,
                title = 'Select One PDB',
                buttons = ('OK', 'Cancel'),
                defaultbutton = 'OK',
                scrolledlist_labelpos = 'n',
                label_text = 'This Domain folder has '+str(len(self.currentSelection.pdbFiles))+' PDB files, \nPlease Select One',
                scrolledlist_items = self.currentSelection.pdbFiles,
                command = self.displayPDB)
            self.pdbSelectdialog.activate()
        else:
            self.displayPDB(0)
    
    
    #Based on the selection, display the PDB file
    def displayPDB(self, selection):		
        if selection == 0:
            selectedFile = self.currentSelection.pdbFiles[0]
            self.pdbPath = os.path.join(self.currentSelection.thisPath,self.currentSelection.thisNode,selectedFile)
        elif selection == 'OK':
            if len(self.pdbSelectdialog.getcurselection()) == 0:
                dialog = Pmw.MessageDialog(self, title = 'Incorrect Selection',
                    defaultbutton = 0,
                    buttons = ('OK', ),
                    message_text = 'Please select one PDB or press \'CANCEL\'')
                dialog.activate()			
                print "Select a PDB file"
                return
            else:
                selectedFile = self.pdbSelectdialog.getcurselection()[0]
                self.pdbPath = os.path.join(self.currentSelection.thisPath,self.currentSelection.thisNode,selectedFile)
        if selection in ['OK', 0]:
            if len(self.pdbPath) > 0:
                try:
                    self.spade_level.pw
                except AttributeError:
                    pass
                else:
                    self.spade_level.pw.forget()
                self.spade_level.string = self.pdbPath #identifier
                self.spade_level.system = MolecularSystem.System(self.pdbPath)
                self.toplevelwindow.title(self.spade_level.system.header)
                self.spade_level.pw = Pmw.PanedWidget(self.spade_level.panes.pane('bottom'),orient='vertical')
                self.spade_level.viewer_pane  = self.spade_level.pw.add('viewer', min=.1,max=.9,size=.88)
                self.spade_level.info_pane = self.spade_level.pw.add('info', min=.1,max=.9,size=.12)
                self.spade_level.info_pane.parent = self.spade_level
                self.spade_level.pw.setnaturalsize()
                self.spade_level.pw.pack(side=TOP,expand=1,fill='both')
                
                self.spade_level.application_notebook = Pmw.NoteBook(self.spade_level.info_pane, borderwidth=1, tabpos='n')
                self.spade_level.application_notebook.pack(expand=YES, fill=BOTH)
                app_keys = self.spade_level.application_pages.keys()
                self.spade_level.application_pages = {}
                self.spade_level.application_pages['Info'] = self.spade_level.application_notebook.add('Info')
                #self.spade_level.application_notebook.setnaturalsize()
                new_cardframe = CardFrame(self.spade_level.application_pages[self.spade_level.application_pages.keys()[0]], self.spade_level.system, self.spade_level)
                new_cardframe.pack(expand=NO, fill=X)
                self.spade_level.cardframe = new_cardframe
                for app_to_load in app_keys:
                    pass
                try:
                    self.spade_level.viewer
                except AttributeError:
                    self.spade_level.viewer = MolecularViewer.MolecularViewer(self.spade_level.viewer_pane, None)
                self.spade_level.viewer.loadSystem(self.spade_level.system)
        if selection != 0:
            self.pdbSelectdialog.deactivate(selection)			
