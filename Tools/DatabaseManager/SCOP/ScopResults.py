
from Tkinter import *
from Bio.SCOP import Node
from Bio import SCOP
from os import *
from shutil import copy
import re
import Pmw
import urllib
import parms
#from ScopFrame import ScopFrame

import MolecularSystem
from GUICards import *
from BlissMolecularViewer import *

#location of plusImage, minusImage
gifs_dir = "C:\\CourseWork\\CS499\\MyTest\\"
pdb_dir = "C:\\CourseWork\\CS499\\MyTest\\SCOP_Pdbs\\"

"""
This class represents each SCOP Node item in the treeWinText textbox.
Creates a frame with:
	 1.+/-button, 
	 2.Empty Label for padding, and  
	 3.Pointer to the actual Node in SCOP Data structure

KNOWN ISSUE: If the search returns a node of type 'CLASS' the '-' button 
	     for the toplevel node does not function as expected. 
	     It works as expected for all other types
"""
class ResultsScopNode(Frame):
	def __init__(self, obj, id, parent=None, toplevel=None):
		Frame.__init__(self, parent, bg='white')
		self.pack()
		self.toplevel = toplevel
		self.thisNode = Node()
		self.thisNode = obj
		self.isDisplayed = 1
		self.isExpanded = 0
		self.isALeaf = 0
		self.id = id   #### EACH SCOP Struct. in the profile is assigned an ID
		self.gifs_dir = parms.get('gifs_dir')
		self.pdb_dir = parms.get('pdb_dir')

		self.createButton()
		self.createLabel()


	def createButton(self):
		self.plusImage = PhotoImage(file=self.gifs_dir+"plusnode.gif")
		self.minusImage = PhotoImage(file=self.gifs_dir+"minusnode.gif")
		
		#Create an empty Label to pad according to Hierarchy
		self.emptyLabel = Label(self, text='', bg='white')
		self.emptyLabel.pack(side=LEFT, padx= 2 * self.toplevel.nameHierarchy[self.thisNode.type])
		if self.thisNode.children:
			self.button = Button(self, image=self.plusImage, command=(lambda node=self, id=self.id: self.toplevel.selectNode(node)))		
		else:
			self.button = Button(self, image=self.minusImage, command=(lambda node=self, id=self.id: self.toplevel.selectNode(node)))
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
	#If the ResultsNode's top level node is not of type 'CLASS' and when no information is available from the 
	#SCOP structure, the file _hie.txt is parsed for the toplevel lineage	
	def displayLineage(self):
		type = {'cf': 'Fold', 'cl': 'Class', 'dm': 'Protein', 'fa': 'Family', 'px': 'Domain', 'sf': 'Superfamily', 'sp': 'Species', '':''}  
		self.toplevel.lineageBox.config(state=NORMAL)
		self.toplevel.lineageBox.delete(1.0, END)
			
		currentNode = self.thisNode
		while len(currentNode.type) > 0:
			self.toplevel.lineageBox.insert('1.0', type[currentNode.type]+":  "+currentNode.description+'\n')
			indx = '1.'+str(len(type[currentNode.type]))
			currentNode = currentNode.parent
			
			self.toplevel.lineageBox.tag_add('label', '1.0', indx)
			self.toplevel.lineageBox.tag_config('label', font=('veranda', 8, 'bold'))
		#The list self.nodesLineage has the all the toplevel lineage. This list is created when loading the SCOP Files
		for lines in self.toplevel.nodesLineage[self.id-1]:
			indx = '1.'+str(lines.find(':'))
			lines = lines[:lines.find('!&!')]
			self.toplevel.lineageBox.insert('1.0', lines+'\n')
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
			


class ResultsFrame(Frame):
	#Arguments-
	#	pth-This the profile name. If the profile is not saved, then pth points to the  
	#	    location where searchResults are stored
	#	parent-This is the paned window object
	#	toplevel-This is the actual Viewer window (Toplevel()). Only needed to change the title of the window
	#	viewer-This is the SCOP Domain Viewer object. Note that a SCOP Domain viewer object should be loaded 
	#	       to load profiles

	def __init__(self, pth, parent=None, toplevel=None, viewer=None):
		Frame.__init__(self, bg='white')
		self.focus()	
		self.scopViewer = viewer
		#The top pane in the viewer represents the SCOP Browser frame
		self.topFrame = parent.pane('top')
		
		#The bottom pane in the viewer represents the Molecular Viewer
		self.toplevel = parent.pane('bottom')
		self.viewer_window = toplevel
		self.dirpath = pth
		
		self.isCurrentProfileSaved = 0
		self.currentProfileName = ""
		
		if self.dirpath.find(parms.get('profilesPath')) == 0 and self.dirpath.find(parms.get('saveResultsin')) != 0:
			pt = self.dirpath
			if pt[len(pt)-1] == '\\':
				self.currentProfileName = self.dirpath[self.dirpath[:-1].rindex('\\')+1:-1]
			else:
				self.currentProfileName = self.dirpath[self.dirpath.rindex('\\')+1:]
			self.isCurrentProfileSaved = 1
			print "HERE++++++++++++"+self.currentProfileName+"  --  "+self.dirpath
		self.nodesList = []
		self.nodesNames = []
		self.removedNames = []
		self.nodesLineage = []
		self.nodesScop = []
		self.domainsScop = []
		self.nameHierarchy = parms.get('nameHierarchy')
		
		if self.dirpath[len(self.dirpath)-1] != '\\':
			self.dirpath = self.dirpath + '\\'
		#Loads the profile/search results from the path given
		self.getAllNodes()
		print "Nodes SCOP: "+str(len(self.nodesScop))
		self.balloon = Pmw.Balloon(self.topFrame)
		self.menuBar = Pmw.MenuBar(self.topFrame, hull_relief=RAISED, hull_borderwidth=1, balloon=self.balloon)
		self.menuBar.pack(fill=X, expand=NO, anchor=N)
		
		self.menuBar.addmenu('Profiles', 'Search Profiles')
		self.menuBar.addmenuitem('Profiles', 'command', command=self.saveCurrentProfile, label='Save Profile')
		self.menuBar.addmenuitem('Profiles', 'command', command=self.loadProfile, label='Load Profile')
		self.menuBar.addmenuitem('Profiles', 'command', command=self.editProfile, label='Edit Profile(s)')

		self.menuBar.addmenu('Edit', 'Edit Results')
		self.menuBar.addmenuitem('Edit', 'command', command=self.removeSelectedNode, label='Remove Node')

		self.menuBar.addmenu('SCOP', 'Search Entire SCOP Database')
		self.menuBar.addmenuitem('SCOP', 'command', label='View Selected Node in SCOP Domain Viewer', command=self.dispInScopDomViewer)

		self.treeWinText = Pmw.ScrolledText(self.topFrame, labelpos = 'n', usehullsize = 1, label_text='SCOP Search Results',
				hull_width = 100, hull_height = 20, text_padx = 1, text_pady = 1, text_wrap='none', text_cursor='arrow')
		self.treeWinText.pack(side=LEFT, expand=YES, fill=BOTH)
		self.pack(expand=YES, fill=BOTH)

		self.lineageBox = Text(self.topFrame, relief=FLAT, cursor = "arrow")
		self.lineageBox.pack(side=TOP, anchor=N)
		self.lineageBox.config(height=12, width=38, state=DISABLED)
		
		self.viewButton = Button(self.topFrame, text="View", state=DISABLED, height=1, command=self.displaySelected, width=5)
		self.viewButton.pack(side=LEFT, anchor=S, padx=10)

		self.scopNodeId = 0
		self.lines = 1
		#Load all the toplevel nodes (all of type 'Class') into the textbox
		for item in self.nodesScop:
			self.scopNodeId = self.scopNodeId + 1
			self.nodesList.append(ResultsScopNode(item.root, self.scopNodeId, self.topFrame, self))
			indx = "%0.1f"%(self.lines)
			self.treeWinText.window_create(indx, window=self.nodesList[len(self.nodesList)-1])
			indx = "%0.1f"%(self.lines+0.1)
			self.treeWinText.insert(indx, '\n')
			self.lines = self.lines + 1	
		
		#self.currentSelection is updated everytime when a user clicks on a Node in the treeWinText text box
		self.currentSelection = self.nodesList[0]	
		self.treeWinText.pack()
		self.treeWinText.configure(text_state='disabled')
	
	#When the user selects to 'View the Selected Domain in SCOP Viewer', the following function
	#builds the lineage list and calls expandToSelectedNode in the Viewer
	def dispInScopDomViewer(self):
		curNode = self.currentSelection.thisNode
		ling = []
		if curNode.sunid == '0':
			print "Cannot Select Identification Node"
		else:
			ling.append(curNode.description+'!&!'+curNode.sunid)
			while curNode.parent:
				curNode = curNode.parent
				ling.append(curNode.description+'!&!'+curNode.sunid)	
			
			ling.pop()
			for item in self.nodesLineage[self.currentSelection.id-1]:
				item = item[(item.find(':')+1):].strip()
				ling.append(item)
			self.scopViewer.expandToSelectedNode(ling)
	
	#This function initiates saving a node as a separate profile.
	def saveCurrentProfile(self):
	    self.nameProfileDlg= Pmw.PromptDialog(self,
	        title = 'Save Current Profile As...',
	        label_text = 'Enter the Name of the Profile:',
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
			print self.currentProfileName+'   IN get'
			if path.isdir(parms.get('profilesPath')+self.currentProfileName):
				dialog = Pmw.MessageDialog(self, title = 'Profile Name Rejected',
					defaultbutton = 0,
					buttons = ('OK', ),
					message_text = 'A Profile with name \''+self.currentProfileName+'\' already exists!')
				dialog.activate()
				self.currentProfileName = ""
				self.saveCurrentProfile()
			else:
				self.finishOffSaving()
	#After a valid profile name is selected, create SCOP Parseable Text files (_cla.txt, 
	#_hie.txt and _des.txt) for the node. With these files a separate SCOP Structure can be 
	#created with the node as a toplevel item.
		
	def finishOffSaving(self):
		print 'In save: '+self.currentProfileName
		if len(self.currentProfileName) > 0:
			#Create a profile directory,
			try:
				mkdir(parms.get('profilesPath')+self.currentProfileName)
			except OSError:
				dialog = Pmw.MessageDialog(self, title = 'Unable to Create Directory',
					defaultbutton = 0,
					buttons = ('OK', ),
					message_text = 'Unable to create the directory: '+parms.get('profilesPath')+self.currentProfileName)
				dialog.activate()				

				dialog = Pmw.MessageDialog(self, title = 'Profile Not Saved',
					defaultbutton = 0,
					buttons = ('OK', ),
					message_text = 'Unable to Save Profile')
				dialog.activate()
				return
			
			#Copy all the results from the temperory directory to the actual profile directory
			files = listdir(parms.get('saveResultsin'))
			for name in self.nodesNames:
				if self.removedNames.count(self.nodesNames.index(name)+1) == 0 and files.count(name+'_hie.txt') > 0 and files.count(name+'_cla.txt') > 0 and files.count(name+'_des.txt') > 0:
					copy(parms.get('saveResultsin')+name+'_hie.txt', parms.get('profilesPath')+self.currentProfileName+'\\'+name+'_hie.txt')
					copy(parms.get('saveResultsin')+name+'_des.txt', parms.get('profilesPath')+self.currentProfileName+'\\'+name+'_des.txt')
					copy(parms.get('saveResultsin')+name+'_cla.txt', parms.get('profilesPath')+self.currentProfileName+'\\'+name+'_cla.txt')
					remove(parms.get('saveResultsin')+name+'_hie.txt')
					remove(parms.get('saveResultsin')+name+'_des.txt')
					remove(parms.get('saveResultsin')+name+'_cla.txt')
				else:
					print "Unable to save: "+name+" under the profile name: "+self.currentProfileName
			self.copyPDBs(parms.get('profilesPath')+self.currentProfileName)	
			self.isCurrentProfileSaved = 1
			self.treeWinText.configure(label_text = 'Current Profile: '+self.currentProfileName)
	
	#This function creates an independent database for the profile.
	#Creates folder structure based on the hierarchy
	#Copies the associated .ent files from the SCOP Database to current profile location					
	
	def copyPDBs(self, savePath):
		#dirNamesDonotIncludelist = ['(', ')', ':', '-', '.', ',', ';']	
		for ScopItem in self.nodesScop:
			indx = self.nodesNames.index(ScopItem.root.description[len(ScopItem.root.description)-5:])+1
			if self.removedNames.count(indx) == 0:
				
				pos = 0
				#parse the scop structure
				for eachItem in self.domainsScop[indx-1]:
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
						curItem = curItem.parent	
					
					curPath = savePath
					#From the 'hierarchy' list, from top down, each item is checked to see if the
					#folder exists by that name. If it doesnot exist, the folder structure is created.
					#If the folder by that name exists, then the next item in the 'heirarchy' list is checked
					hierarchy.pop()
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
							
					
	def loadProfile(self):
		files = listdir(parms.get('profilesPath'))
		profiles = []
		for item in files:
			direc = parms.get('profilesPath')+item
			if path.isdir(direc):
				profiles.append(item)

		self.neww = Toplevel()
	        dropdown = Pmw.ComboBox(self.neww, label_text = 'Select a Profile:', labelpos = 'nw', selectioncommand=self.loadSelectedProfile,
	        	scrolledlist_items = profiles)
	        dropdown.pack(side = LEFT, anchor = N, fill = X, expand = 1, padx = 8, pady = 8)
		
	
	def loadSelectedProfile(self, selection):
		self.neww.destroy()
		new_window=Toplevel()

		pw = Pmw.PanedWidget(new_window)	
		pw.add('top', min=10)
		pw.add('bottom', min = 10)
			
		new_window.title("SCOP Profile")
		new_window.config(bg='white')
		
		geometry_string = "%dx%d%+d%+d" %(700,800,50,30) # width,height,x-offset,y-offset
		new_window.geometry(geometry_string)		
		
		pw.pane('top').configure(background='white')
		a = ResultsFrame(parms.get('profilesPath')+selection+'\\', pw, new_window, self.scopViewer).pack(side=TOP)
		pw.pack(expand=1, fill=BOTH) 
				
		bottomPane = pw.pane('bottom')
		bottomPane.system = MolecularSystem.System(self, bottomPane)
		bottomPane.system.load_pdb("")
		bottomPane.system.color_ribbon_by_chain()
		bottomPane.viewer = MolnirMolecularViewer(bottomPane, bottomPane.system)
		bottomPane.viewer.loadSystem(bottomPane.system)
		bottomPane.cardframe = CardFrame(bottomPane, bottomPane.system)
		bottomPane.cardframe.pack(expand=NO, fill=X)
		
		self.scopViewer.topMenuSystem.parent.scopViewer_windows.append(new_window)
		self.topFrame.config(bg='white')
		
	#This function is not implemented. The idea is to build a GUI to help the user manage (rename, delete, etc) profiles 
	def editProfile(self):
		print "Not Implemented"	
	
	#Removes the selected node from the treeWinText text box and in self.nodesList 
	#the removedNode acts just as a placeholder
	def removeSelectedNode(self):
		if self.currentSelection.thisNode.parent is None:
			self.removedNames.append(self.currentSelection.id)
			self.treeWinText.configure(text_state='normal')
			pos = 0
			
			while self.nodesList[pos].id != self.currentSelection.id:
				pos = pos + 1
			
			while len(self.nodesList) > pos and self.nodesList[pos].id == self.currentSelection.id:
				self.nodesList.pop(pos)
				indx = "%0.1f"%(pos+1)
				indx1 = "%0.1f"%(pos+2)
				self.treeWinText.delete(indx, indx1)

			self.treeWinText.configure(text_state='disabled')
			self.currentSelection = self.nodesList[0]

	#Given the path to the profile, this function loads all the SCOP Structures and 
	#also creates a list with lineage lists for each SCOP Structure.			
	def getAllNodes(self):
		files = listdir(self.dirpath)
		for item in files:
			if re.compile('\d\d\d\d\d_hie.txt').match(item):
				if path.isfile(self.dirpath+item[:5]+'_des.txt') and path.isfile(self.dirpath+item[:5]+'_cla.txt'):
					self.nodesNames.append(item[:5])
					clasi = file(self.dirpath+item[:5]+'_cla.txt', 'r')
					descr = file(self.dirpath+item[:5]+'_des.txt', 'r')
					hiera = file(self.dirpath+item[:5]+'_hie.txt', 'r')
					self.nodesScop.append(SCOP.Scop(clasi, descr, hiera))
					self.domainsScop.append(self.nodesScop[len(self.nodesScop)-1].getDomains())
					descr.seek(0)
					self.nodesScop[len(self.nodesScop)-1].root.description = descr.readline().split('\t').pop().strip()+'_'+item[:5]			
					hiera.seek(0)
					line = hiera.readline()
					lineage = []
					while line[0] == '#':
						lineage.append(line[1:].strip())
						line = hiera.readline()
					lineage.reverse()
					self.nodesLineage.append(lineage)
					#print self.nodesScop[len(self.nodesScop)-1].root.description
					#print lineage

	def toggleSelection(self, nde):
		self.currentSelection.label.config(bg='white', fg='black')
		nde.label.config(bg='blue', fg='white')
		self.currentSelection = nde
		
	#This function is called when a node is selected for expansion or contraction			
	def selectNode(self, nde):
		self.toggleSelection(nde)
		self.treeWinText.configure(text_state='normal')
		if nde.isExpanded == 1:
			pos = 0
			#Browse to the corresponding SCOP Structure in nodesList list
			while self.nodesList[pos].id != nde.id:
				pos = pos + 1
			
			#Browse to the corresponding Node in the nodesList
			while nde.thisNode.sunid != self.nodesList[pos].thisNode.sunid:
				pos = pos + 1

			self.nodesList[pos].contract()
			pos = pos + 1
			#Remove all the subnodes of the contrated node
			while len(self.nodesList) > pos and  nde.id == self.nodesList[pos].id and self.nameHierarchy[nde.thisNode.type] < self.nameHierarchy[self.nodesList[pos].thisNode.type]:
				self.nodesList.pop(pos)
				indx = "%0.1f"%(pos+1)
				indx1 = "%0.1f"%(pos+2)
				self.treeWinText.delete(indx, indx1)

		else:
			#Expand the selected node		
			pos = 0
			len_nodesList = len(self.nodesList)
			
			while self.nodesList[pos].id != nde.id:
				pos = pos + 1
			
			while pos < len_nodesList:
				#Iterate until the selected node is found in the self.nodesList list
				item = self.nodesList[pos]
				if item.thisNode.sunid == nde.thisNode.sunid and item.id == nde.id:
					chld_pos = 1
					self.nodesList[pos].expand()
					if self.nodesList[pos].thisNode.children:
						for child in self.nodesList[pos].thisNode.children:
							self.nodesList.insert(pos+chld_pos, ResultsScopNode(child, self.nodesList[pos].id, self.topFrame, self))							
							indx = "%0.1f"%(pos+chld_pos+1)
							self.treeWinText.window_create(indx, window=self.nodesList[pos+chld_pos])
							indx = "%0.1f"%(pos+chld_pos+1+0.1)
							self.treeWinText.insert(indx, '\n')
							#self.nodesList[pos+1].pack(anchor=W)
							chld_pos = chld_pos + 1

				pos = pos + 1

		print "Current Size of Nodes List: %d\n"%(len(self.nodesList))
		self.treeWinText.configure(text_state='disabled')

	#This function is not used.
	#Useful when downloading a pdb for a selected domain from the online Protein Data Bank
	def downloadPDB(self):
		
		pdb_path=""
		pos = 0
		while pos < len(self.domainsScop[self.currentSelection.id-1]) and self.currentSelection.thisNode.sunid != self.domainsScop[self.currentSelection.id-1][pos].sunid:
			pos = pos + 1
		
		if pos < len(self.domainsScop[self.currentSelection.id-1]):				
			res = self.domainsScop[self.currentSelection.id-1][pos].residues.pdbid
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
		pdb_path=""
		pos = 0
		while pos < len(self.domainsScop[self.currentSelection.id-1]) and self.currentSelection.thisNode.sunid != self.domainsScop[self.currentSelection.id-1][pos].sunid:
			pos = pos + 1
			
		if self.dirpath.find(parms.get('saveResultsin')) >= 0: #if the domain is to be viewed from a profile 
			pdb_path = parms.get('scopPDBsLocation') + self.domainsScop[self.currentSelection.id-1][pos].sid[2:4] + '\\' + self.domainsScop[self.currentSelection.id-1][pos].sid + '.ent'
			if path.isfile(pdb_path):
				return pdb_path
			print "Selected Domain doesnot exist in the SCOP Database"

		elif self.dirpath.find(parms.get('profilesPath')) >= 0: #if the domain is to be viewed from the ScopDatabase 
			curItem = self.domainsScop[self.currentSelection.id-1][pos]
			tmpPath = ""
			while curItem is not None:
				name = self.filterFoldername(curItem.description)
				if len(name) > 20:
					name = name[:20]+'_'+curItem.sunid
				tmpPath = name.strip() + '\\' + tmpPath
				curItem = curItem.parent
			tmpPath = tmpPath[tmpPath.index('\\')+1:]
			print "TmpPath: "+tmpPath
			pdb_path = parms.get('profilesPath')+ self.currentProfileName + '\\' + tmpPath + self.domainsScop[self.currentSelection.id-1][pos].sid + '.ent'	
			print "PDB PATH: "+pdb_path
			if path.isfile(pdb_path):
				return pdb_path
			print "Selected Domain doesnot exist in the profile"				
		return ""
		
	#This function removes characters that cannot be accepted in the directory names				
	def filterFoldername(self, name):
		dirNamesDonotIncludelist = parms.get('dirNamesDonotIncludelist')
		for itm in dirNamesDonotIncludelist:
			name = name.replace(itm, '')
		return name
	
	#This function displays the selected Node
	def displaySelected(self):			
		self.pdbPath = self.getPDBLocation()
		if len(self.pdbPath) > 0:
			self.toplevel.cardframe.forget()
			self.toplevel.viewer.screen.forget()
			#self.topMenuSystem.loadSystem(self.pdbPath, self.toplevel)			
			
			self.toplevel.string = self.pdbPath 	#identifier
			self.toplevel.system = MolecularSystem.System(self, self.toplevel)
			self.toplevel.system.load_pdb(self.pdbPath)
			self.toplevel.system.color_ribbon_by_chain()
			self.viewer_window.title(self.toplevel.system.header)
			self.toplevel.viewer = MolnirMolecularViewer(self.toplevel, self.toplevel.system)
			self.toplevel.viewer.loadSystem(self.toplevel.system)
			self.toplevel.cardframe = CardFrame(self.toplevel, self.toplevel.system)
			self.toplevel.cardframe.pack(expand=NO, fill=X)
			#self.parent.system_windows.append(self.toplevel) 


		
if __name__ == '__main__':

	root = Toplevel()
	root.title("SCOP Search Results")
	root.config(bg='white')
	
	new_window = root
	
	geometry_string = "%dx%d%+d%+d" %(700,800,50,30) # width,height,x-offset,y-offset
	new_window.geometry(geometry_string)
	
	tp = Frame(root)
	tp.pack(side=TOP)

	pth = 'C:\\CourseWork\\CS499\\MyTest\\SCOPSearchProfiles\\1223'
	a = ResultsFrame(pth, tp, root)

	tp.config(bg='white')
	new_window.system = MolecularSystem.System(new_window, new_window)
	new_window.system.load_pdb("")
	new_window.system.color_ribbon_by_chain()
	new_window.viewer = MolnirMolecularViewer(new_window, new_window.system)
	new_window.viewer.loadSystem(new_window.system)
	new_window.cardframe = CardFrame(new_window, new_window.system)
	new_window.cardframe.pack(expand=NO, fill=X)
	
	
	root.config(bg='white')

	a.mainloop()        
	
		
			