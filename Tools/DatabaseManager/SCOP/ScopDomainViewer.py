import sys
import os
print os.path.abspath(os.path.join(os.getcwd(),os.path.normpath('../')))
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.normpath('../../'))))

from Tkinter import *
from Bio.SCOP import Node
from Bio import SCOP
import urllib
import MolecularSystem
from GUICards import *
from BlissMolecularViewer import *
from ScopSearch import *
from ScopResults import *
import parms
import Pmw
import re

#location of plusImage, minusImage
gifs_dir = parms.get('gifs_dir')
pdb_dir = parms.get('pdb_dir')

"""
This class represents each SCOP Node item in the treeWinText textbox.
Creates a frame with:
	 1.+/-button, 
	 2.Empty Label for padding, and  
	 3.Pointer to the actual Node in SCOP Data structure

"""
class ScopNode(Frame):
	def __init__(self, obj, parent=None):
		Frame.__init__(self, parent, bg='white')
		self.pack()
		self.toplevel= parent
		self.nameHierarchy = parms.get('nameHierarchy')
		self.thisNode = Node()
		self.thisNode = obj
		self.isDisplayed = 1
		self.isExpanded = 0
		self.isALeaf = 0
		self.createButton()
		self.createLabel()
		
	def createButton(self):

		self.plusImage = PhotoImage(file=os.path.join(gifs_dir,"plusnode.gif"))
		self.minusImage = PhotoImage(file=os.path.join(gifs_dir,"minusnode.gif"))
		
		#Create an empty Label add whitespace padding before the frame according to Hierarchy
		self.emptyLabel = Label(self, text='', bg='white')
		self.emptyLabel.pack(side=LEFT, padx= 2 * self.nameHierarchy[self.thisNode.type])
		if self.thisNode.children:
			self.button = Button(self, image=self.plusImage, command=(lambda node=self: self.toplevel.selectNode(node)))			
		else:
			self.button = Button(self, image=self.minusImage, command=(lambda node=self: self.toplevel.selectNode(node)))
			self.isALeaf = 1

		self.button.pack(side=LEFT)
		
		self.button.bind('<Button-1>', self.singleClick)

	def createLabel(self):
		self.label = Label(self, text='  '+self.thisNode.description+'  ', bg='white')
		self.label.pack(side=RIGHT)
		self.label.bind('<Double-1>', self.doubleClick)
		self.label.bind('<Button-1>', self.singleClick)

	def singleClick(self, event):
		self.toplevel.toggleSelection(self)
		if len(self.thisNode.children) == 0:
			self.toplevel.viewButton.config(state=NORMAL)
		else:
			self.toplevel.viewButton.config(state=DISABLED)
		self.displayLineage()

	#This function manages the lineageBox text box. Updates text box with the lineage of the selected node	
	def displayLineage(self):
		type = {'cf': 'Fold', 'cl': 'Class', 'dm': 'Protein', 'fa': 'Family', 'px': 'Domain', 'sf': 'Superfamily', 'sp': 'Species', '':''}  
		self.toplevel.lineageBox.config(state=NORMAL)
		self.toplevel.lineageBox.delete(1.0, END)

		currentNode = self.thisNode
		while len(currentNode.type) > 0:
			indx = '1.'+str(len(type[currentNode.type]))

			self.toplevel.lineageBox.insert('1.0', type[currentNode.type]+":  "+currentNode.description+'\n')
			currentNode = currentNode.parent

			self.toplevel.lineageBox.tag_add('label', '1.0', indx)
			self.toplevel.lineageBox.tag_config('label', font=('veranda', 8, 'bold'))
			
		self.toplevel.lineageBox.insert('1.0', '\t\tLineage\n')				
		self.toplevel.lineageBox.tag_add('label', '1.0', '2.0')
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
		
	
#This class handles all the functionality associated with the SCOPDomainViewer
class ScopViewer(Frame):
	#Arguments-
	#	parent-This is the paned window object
	#	viewer-This is the actual Viewer window (Toplevel()). Only needed to change the title of the window
	#	menSystem-This is the MenuSystem object. The only objective of this 
	#		  is to add newly created viewer/results windows 
	#		  to menSystem.scopViewer_windows list in the MenuSystem object.
	def __init__(self, parent, viewer, menSystem):
		Frame.__init__(self, parent.pane('top'), bg='white')
		self.focus()	
		self.nodesList = []  #List of all the nodes in the treeWinText textbox. (has ScopNode objects)
		self.currentNodeClaList = [] #Temp. variable used in UpdateCla funciton
		self.entireClaList = {} #Classification file is parsed and all domains are hashed by their sids
		
		#The bottom pane in the viewer represents the Molecular Viewer
		self.molnir_level = parent.pane('bottom')
		self.viewer_window = viewer	
		self.topMenuSystem = menSystem
		self.saveResultsin = ""      #directory path set by the user to save the profile
		self.currentProfileName = ""
		self.profilesPath = parms.get('profilesPath')
		self.type = parms.get('classifications_fullname')
		
		#Open the SCOP Parseable Text files
		clasi = file(parms.get('scopClassification'), 'r')
		descr = file(parms.get('scopDescription'), 'r')
		heira = file(parms.get('scopHierarchy'), 'r')
		
		#Load SCOP Structure and SCOP Domains
		self.scopDataStruct = SCOP.Scop(clasi, descr, heira)    	 
		self.domains = self.scopDataStruct.getDomains()
		
		self.nameHierarchy = parms.get('nameHierarchy')
		self.pdbPath = ""

		self.balloon = Pmw.Balloon(self)
		self.menuBar = Pmw.MenuBar(self, hull_relief=RAISED, hull_borderwidth=1, balloon=self.balloon)
		self.menuBar.pack(fill=X, expand=NO, anchor=N)
		
		self.menuBar.addmenu('Profiles', 'Search Profiles')
		self.menuBar.addmenuitem('Profiles', 'command', command=self.loadProfile, label='Load Profile')
		self.menuBar.addmenuitem('Profiles', 'command', command=self.editProfile, label='Edit Profile(s)')

		self.menuBar.addmenu('Save', 'Create a local DB Structure')
		self.menuBar.addmenuitem('Save', 'command', command=self.saveSelectedNode, label='Save selected Node as a Profile')

		self.menuBar.addmenu('Search', 'Search Entire SCOP Database')
		self.menuBar.addmenuitem('Search', 'command', command=self.loadSCOPSearch, label='Search SCOP')
		
		self.treeWinText = Pmw.ScrolledText(self, labelpos = 'n', usehullsize = 1, label_text='SCOP Domain Viewer',
				hull_width = 100, hull_height = 20, text_padx = 1, text_pady = 1, text_wrap='none', text_cursor='arrow')
		self.treeWinText.pack(side=LEFT, expand=YES, fill=BOTH)
		self.pack(expand=YES, fill=BOTH)	
		
		self.lineageBox = Text(self, relief=FLAT, cursor = "arrow")
		self.lineageBox.pack(side=TOP, anchor=N)
		self.lineageBox.config(height=12, width=38, state=DISABLED)
		
		self.viewButton = Button(self, text="View", state=DISABLED, height=1, command=self.displaySelected, width=5)
		self.viewButton.pack(side=LEFT, anchor=S, padx=10)
		
		self.lines = 1
		
		#Load all the toplevel nodes (all of type 'Class') into the textbox
		for item in self.scopDataStruct.root.children:
			self.nodesList.append(ScopNode(item, self))
			indx = "%0.1f"%(self.lines)
			self.treeWinText.window_create(indx, window=self.nodesList[len(self.nodesList)-1])
			indx = "%0.1f"%(self.lines+0.1)
			self.treeWinText.insert(indx, '\n')
			self.lines = self.lines + 1
		
		#self.currentSelection is updated everytime when a user clicks on a Node in the treeWinText text box
		self.currentSelection = self.nodesList[0]
		self.treeWinText.pack()
		self.treeWinText.configure(text_state='disabled')
		self.pack()		
	
	#This function initiates saving a node as a separate profile.
	def saveSelectedNode(self):
		self.nameProfileDlg= Pmw.PromptDialog(self,
			title = 'Save Current Profile As...',
			label_text = 'The selected node will be saved in the current profiles direcotry. \nEnter the Name of the Profile:',
			entryfield_labelpos = 'n',
			defaultbutton = 0,
			buttons = ('OK', 'Cancel'),
			command = self.validateName)
	
	#Check to see if a profile with the same name already exists in the Profiles Directory
	def validateName(self, result):
		if result is None or result == 'Cancel':
			self.nameProfileDlg.withdraw()
		elif result == 'OK':
			self.currentProfileName = self.nameProfileDlg.get()
			self.nameProfileDlg.withdraw()
			
			#If a profile with the entered name already exists, prompt again
			if path.isdir(parms.get('profilesPath')+self.currentProfileName):
				dialog = Pmw.MessageDialog(self, title = 'Profile Name Rejected',
					defaultbutton = 0,
					buttons = ('OK', ),
					message_text = 'A Profile with name \''+self.currentProfileName+'\' already exists!')
				dialog.activate()
				self.currentProfileName = ""
				self.saveSelectedNode()
			else:
			#otherwise, create a directory and initiate the SCOP Text file and database creattion
				self.saveResultsin = parms.get('profilesPath')+self.currentProfileName
				mkdir(self.saveResultsin)
				self.saveResultsin = self.saveResultsin + '\\'
				self.finishOffSaving()
	
	#After a valid profile name is selected, create SCOP Parseable Text files (_cla.txt, 
	#_hie.txt and _des.txt) for the node. With these files a separate SCOP Structure can be 
	#created with the node as a toplevel item.
	
	def finishOffSaving(self):
		if len(self.entireClaList) == 0:
			self.buildClaList() #If the original 'Cla' file is not loaded, load it
		
		node = self.currentSelection.thisNode 
		#if the node has a sunid:45678, then the parseable files are named 45678_hie.txt, etc
		cla = file(self.saveResultsin+str(node.sunid)+"_cla.txt", 'w')
		hie = file(self.saveResultsin+str(node.sunid)+"_hie.txt", 'w')
		des = file(self.saveResultsin+str(node.sunid)+"_des.txt", 'w')

		#From the newly created SCOP Parseable files, there is no way to retreive the
		#full lineage of the node. So, save the lineage of the selected node's parent
		#into the hierarchy file as comments in it. Later when a profile is loaded, this
		#information is used to display the full lineage of the node (even when the node's 
		#parent is not available)
		
		ling_item = node.parent
		lineage_list = []
		while len(ling_item.type) > 0:
			strin = '#'+self.type[ling_item.type]+': '+ling_item.description+'!&!'+ling_item.sunid+'\n'
			lineage_list.insert(0, strin)
			ling_item = ling_item.parent
		hie.writelines(lineage_list)
		
		#Create the SCOP Parseable text files for the node
		hieRec = str(node.toHieRecord())
		hie.write('0\t-\t'+node.sunid+'\n')
		hie.write(hieRec[:hieRec.index('\t')+1]+'0'+hieRec[hieRec.index('\t')+1:][hieRec[hieRec.index('\t')+1:].index('\t'):])
		des.write(str(node.toDesRecord()))
		self.updateClaList(node, cla)
		
		for ndes in node.children:
			hieRec = str(ndes.toHieRecord())
			#hie.write(hieRec[:hieRec.index('\t')+1]+'0'+hieRec[hieRec.index('\t')+1:][hieRec[(hieRec.index('\t')+1):].index('\t'):])
			hie.write(hieRec)
			des.write(str(ndes.toDesRecord()))
			self.updateClaList(ndes, cla)
			for nde in ndes.children:
				self.createHieCla(nde, hie, des, cla)
		cla.close()
		des.close()
		hie.close()
		#Create a folder structure based on the heirarchy of the node
		self.copyPDBs(self.saveResultsin)
	
	#Writes the hierarchy and description records to the parseable files
	def createHieCla(self, nde, hie, des, cla):
		hie.write(str(nde.toHieRecord()))
		des.write(str(nde.toDesRecord()))
		
		self.updateClaList(nde, cla)
		for ndes in nde.children:
			self.createHieCla(ndes, hie, des, cla)
	
	#If the node is a domain, it writes the classification for that domain
	def updateClaList(self, nde, cla):
		if len(nde.children) == 0:
			if self.currentNodeClaList.count(nde.sunid) == 0:
				self.currentNodeClaList.append(nde.sunid)
				cla.write(self.entireClaList[nde.sunid])
	
	#Builds a hashed list with the sunids as keys from the classification file
	def buildClaList(self):
		fullFile = file(parms.get('scopClassification'), 'r').readlines()
		rex = re.compile('\t\d\d\d\d\d\t')
		for line in fullFile:
			a = rex.search(line)
			if a:
				self.entireClaList[a.group().strip()] = line
	
	#This function creates an independent database for the profile.
	#Creates folder structure based on the hierarchy
	#Copies the associated .ent files from the SCOP Database to current profile location		
	def copyPDBs(self, savePath):
		cla = file(self.saveResultsin+str(self.currentSelection.thisNode.sunid)+"_cla.txt", 'r')
		hie = file(self.saveResultsin+str(self.currentSelection.thisNode.sunid)+"_hie.txt", 'r')
		des = file(self.saveResultsin+str(self.currentSelection.thisNode.sunid)+"_des.txt", 'r')
		scopStruct = SCOP.Scop(cla, des, hie)
		domsStruct = scopStruct.getDomains()
		
		for eachItem in domsStruct:
			hierarchy = []
			curItem = eachItem
			"""
			This loop creates a list of directories need to save the domain
			Eg: If the node hierarchy is of the form
				Globin->Globin_child->Globin_child_child->Globin_child_child_child
				and 'Globin_child_child_child' is selected
				then the list hierarchy will contain
					[Globin, Globin_child, Globin_child_child, Globin_child_child_child]
			"""		
			
			while curItem is not None:
				name = self.filterFoldername(curItem.description)
				#If the length of the node.description is > 20, the pad the end of it 
				#with its sunid 
				if len(name)>20:
					name = name[:20]+'_'+curItem.sunid
					
				hierarchy.append(name.strip())
				if self.currentSelection.thisNode.sunid.find(curItem.sunid) == 0:
					curItem = None	
				else:
					curItem = curItem.parent 
				
			curPath = savePath
			#From the 'hierarchy' list, from top down, each item is checked to see if the
			#folder exists by that name. If it doesnot exist, the folder structure is created.
			#If the folder by that name exists, then the next item in the 'heirarchy' list is checked
			while len(hierarchy) > 0:
				foldr = hierarchy.pop()
				curPath = curPath+'\\'+foldr
				if not path.isdir(curPath):
					mkdir(curPath)
					while len(hierarchy) > 0:
						curPath = curPath + '\\' + hierarchy.pop()
						mkdir(curPath)
			
			scopPdbsLoc = parms.get('scopPDBsLocation')
			#Copy the .ent file from the SCOP Database to the current profile location		
			if path.isdir(scopPdbsLoc+eachItem.sid[2:4]) and path.isfile(scopPdbsLoc+eachItem.sid[2:4]+'\\'+eachItem.sid+'.ent'):
				copy(scopPdbsLoc+eachItem.sid[2:4]+'\\'+eachItem.sid+'.ent', curPath+'\\'+eachItem.sid+'.ent')	
			else:
				print "Protein Domain: "+eachItem.sid+" doesnot exist!"
	
	#This function removes characters that cannot be accepted in the directory names
	def filterFoldername(self, name):
		dirNamesDonotIncludelist = parms.get('dirNamesDonotIncludelist')	
		for itm in dirNamesDonotIncludelist:
			name = name.replace(itm, '')
		return name.strip()	
	
	#This function is called from the SCOP Results Viewer Window. When a node in the 
	#results window is selected and 'View Selected Node in SCOP Domain Viewer is selected'
	#this funciton is called			
	def expandToSelectedNode(self, lineage):
		#Begin to update the treeWinText text box. 
		#Clear the text box and empty 'nodesList' list
		self.treeWinText.configure(text_state='normal')	
		self.treeWinText.delete('1.0', 'end')
		while len(self.nodesList) > 0:
			self.nodesList.pop()
		self.focus()

		self.lines = 1
		#Begin creating new Nodes in the treeWinText textbox
		for item in self.scopDataStruct.root.children:
			self.nodesList.append(ScopNode(item, self))
			indx = "%0.1f"%(self.lines)
			self.treeWinText.window_create(indx, window=self.nodesList[len(self.nodesList)-1])
			indx = "%0.1f"%(self.lines+0.1)
			self.treeWinText.insert(indx, '\n')
			self.lines = self.lines + 1
		
		#Lineage list has the list of the nodes with each node of format <descr>!&!<sunid>
		#Browse to the selected node from the top, by selecting one node at a time
		self.currentSelection = self.nodesList[0]
		while len(lineage) > 0:
			item = lineage.pop()
			sunid = item[(item.find('!&!')+3):]
			description = item[:item.find('!&!')]
			pos = 0
			while pos < len(self.nodesList):
				if self.nodesList[pos].thisNode.sunid.find(sunid) == 0:
					self.selectNode(self.nodesList[pos])
					self.nodesList[pos].displayLineage()
					pos = len(self.nodesList)
				pos = pos + 1		
	
	#This function is called when a node is selected for expansion or contraction	
	def selectNode(self, nde):
		#Select the 'nde' Node
		self.toggleSelection(nde)
		self.treeWinText.configure(text_state='normal')
		if nde.isExpanded == 1:
			pos = 0
			#Browse to the node in nodesList list
			while nde.thisNode.sunid != self.nodesList[pos].thisNode.sunid:
				pos = pos + 1

			self.nodesList[pos].contract()
			pos = pos + 1
			
			#Remove all the subnodes of the contrated node
			#This removal is done by: removing the nodes following the previously contracted node
			#until a node with higher hierarchy(type='cf', etc)  is reached
			while len(self.nodesList) > pos and  self.nameHierarchy[nde.thisNode.type] < self.nameHierarchy[self.nodesList[pos].thisNode.type]:
				self.nodesList.pop(pos)
				indx = "%0.1f"%(pos+1)
				indx1 = "%0.1f"%(pos+2)
				self.treeWinText.delete(indx, indx1)
							
		else:
			#Expand the selected node
			pos = 0
			len_nodesList = len(self.nodesList)
			while pos < len_nodesList:
				item = self.nodesList[pos]
				#Iterate until the selected node is found in the self.nodesList list
				if item.thisNode.sunid == nde.thisNode.sunid:
					chld_pos = 1
					self.nodesList[pos].expand()
					#Now display all its subnodes in the treeWinText text box
					if self.nodesList[pos].thisNode.children:
						for child in self.nodesList[pos].thisNode.children:
							self.nodesList.insert(pos+chld_pos, ScopNode(child, self))							
							indx = "%0.1f"%(pos+chld_pos+1)
							self.treeWinText.window_create(indx, window=self.nodesList[pos+chld_pos])
							indx = "%0.1f"%(pos+chld_pos+1+0.1)
							self.treeWinText.insert(indx, '\n')
							chld_pos = chld_pos + 1
													
				pos = pos + 1
				
		self.treeWinText.configure(text_state='disabled')

	#Change the highlighted color of the previously highlighted and current selected node.	
	def toggleSelection(self, nde):
		self.currentSelection.label.config(bg='white', fg='black')
		nde.label.config(bg='blue', fg='white')
		self.currentSelection = nde

	#Loads the SCOP search object
	def loadSCOPSearch(self):
		root = Toplevel()
		a = SearchFrame(root, self)
	
	#Follows the user request to load a profile. Lists all the existing profiles and prompts the user to select one
	def loadProfile(self):
		files = listdir(self.profilesPath)
		profiles = []
		for item in files:
			direc = self.profilesPath+item
			if path.isdir(direc):
				profiles.append(item)

		self.dropdownLevel = Toplevel()
	        self.dropdown = Pmw.ComboBox(self.dropdownLevel, label_text = 'Select a Profile:', labelpos = 'nw', selectioncommand=self.loadSelectedProfile,
	        	scrolledlist_items = profiles)
	        self.dropdown.pack(side = LEFT, anchor = N, fill = X, expand = 1, padx = 8, pady = 8)
		
	#Following the selection of a profile, initiate loading the profile into the ResultsFrame
	def loadSelectedProfile(self, selection):
	
		self.dropdownLevel.destroy() #Destroy the dropdown menu box
		new_window=Toplevel()
		
		pw = Pmw.PanedWidget(new_window)	
		pw.add('top', min=10)
		pw.add('bottom', min = 10)
		
		new_window.title("SCOP Profile")
		new_window.config(bg='white')
		
		geometry_string = "%dx%d%+d%+d" %(700,800,50,30) # width,height,x-offset,y-offset
		new_window.geometry(geometry_string)
		
		pw.pane('top').configure(background='white')
		a = ResultsFrame(self.profilesPath+selection, pw, new_window, self).pack(side=TOP)
		pw.pack(expand=1, fill=BOTH) 
		
		bottomPane = pw.pane('bottom')
		
		bottomPane.system = MolecularSystem.System(self, bottomPane)
		bottomPane.system.load_pdb("")
		bottomPane.system.color_ribbon_by_chain()
		
		bottomPane.viewer = MolnirMolecularViewer(bottomPane, bottomPane.system)
		bottomPane.viewer.loadSystem(bottomPane.system)
		bottomPane.cardframe = CardFrame(bottomPane, bottomPane.system)
		bottomPane.cardframe.pack(expand=NO, fill=X)
		self.topMenuSystem.parent.scopViewer_windows.append(new_window)
	
	#This function is not implemented. The idea is to build a GUI to help the user manage (rename, delete, etc) profiles 
	def editProfile(self):
		print "Not Implemented"	
	
	#This function is not used.
	#Useful when downloading a pdb for a selected domain from the online Protein Data Bank
	def downloadPDB(self):
		pos = 0
		while self.currentSelection.thisNode.sunid != self.domains[pos].sunid:
			pos = pos + 1
		
		res = self.domains[pos].residues.pdbid
		pdb_path = pdb_dir+res+".pdb"
		
		if not path.exists(pdb_path):
			url = 'http://www.rcsb.org/pdb/cgi/export.cgi/'+res+'.pdb?job=download;pdbId='+res+';page=0;opt=show;format=PDB;pre=1&compression=None'
			urlfile = urllib.urlopen(url)

			writeto = file(pdb_path, 'w')
			lines = urlfile.readlines()
			for line in lines:
				writeto.write(line)
			writeto.close()
		return pdb_path
	
	#This function returns the location of the '.ent' file for the selected node
	#If it doesnot exist then an empty string is returned
	def getPDBLocation(self):
		#check if the path includes the profiles directory path
		pos = 0
		while self.currentSelection.thisNode.sunid != self.domains[pos].sunid:
			pos = pos + 1
		
		pdb_path = parms.get('scopPDBsLocation') + self.domains[pos].sid[2:4] + '\\' + self.domains[pos].sid + '.ent'
		if path.isfile(pdb_path):
			return pdb_path
				
		return ""
	
	#This function displays the selected Node
	def displaySelected(self):			
		self.pdbPath = self.getPDBLocation()
		if len(self.pdbPath) > 0:
			self.molnir_level.cardframe.forget()
			self.molnir_level.viewer.screen.forget()			
			
			self.molnir_level.string = self.pdbPath 	#identifier
			self.molnir_level.system = MolecularSystem.System(self, self.molnir_level)
			self.molnir_level.system.load_pdb(self.pdbPath)
			self.molnir_level.system.color_ribbon_by_chain()
			self.viewer_window.title(self.molnir_level.system.header)
			self.molnir_level.viewer = MolnirMolecularViewer(self.molnir_level, self.molnir_level.system)
			self.molnir_level.viewer.loadSystem(self.molnir_level.system)
			self.molnir_level.cardframe = CardFrame(self.molnir_level, self.molnir_level.system)
			self.molnir_level.cardframe.pack(expand=NO, fill=X)
			
		else:
			print "Selected Node does not have a corresponding .ent file in the Local Database!"
		
		

