from Tkinter import *
import urllib
import re
import Pmw

from os import *
import MolecularSystem
from GUICards import *
from MolnirMolecularViewer import *
#from Molnir import CardFrame

#global plusImage, minusImage
gifs_dir = "C:\\CourseWork\\CS499\\MyTest\\"
root_path = "C:\\CourseWork\\CS499\\MyTest\\"


###########################################################
class DirNode(Frame):
	def __init__(self, path, obj, parent):
		Frame.__init__(self, parent, bg='white')
		self.pack()
		self.toplevel = parent
		self.thisPath = path
		self.thisNode = obj
		self.isDisplayed = 1
		self.isExpanded = 0
		self.isALeaf = 0
		self.createButton()
		self.createLabel()

		
	def createButton(self):
		self.plusImage = PhotoImage(file=gifs_dir+"plusnode.gif")
		self.minusImage = PhotoImage(file=gifs_dir+"minusnode.gif")
		
		if path.isdir(self.thisPath+self.thisNode):
			self.button = Button(self, image=self.plusImage, command=(lambda node=self: self.toplevel.selectNode(node)))			
		else:
			self.button = Button(self, image=self.minusImage, command=(lambda node=self: self.toplevel.selectNode(node)))
			self.isALeaf = 1
		regex = re.compile('\\\\')
		padding = len(regex.split(self.thisPath))- len(regex.split(self.toplevel.root))
		print "root"+self.toplevel.root
		print "chil"+self.thisPath
		self.button.pack(side=LEFT, padx = padding * 5)
		
		self.button.bind('<Button-1>', self.singleClick)

	def createLabel(self):
		self.label = Label(self, text='  '+self.thisNode+'  ', bg='white')
		self.label.pack(side=RIGHT)
		self.label.bind('<Double-1>', self.doubleClick)
		self.label.bind('<Button-1>', self.singleClick)

	def singleClick(self, event):
		self.toplevel.toggleSelection(self)
		if path.isfile(self.thisPath+self.thisNode):
			self.toplevel.viewButton.config(state=NORMAL)
		else:
			self.toplevel.viewButton.config(state=DISABLED)
		
	def doubleClick(self, event):
		self.toplevel.selectNode(self)
		
	def destroyDisplay(self):
		isDisplayed = 0
		self.button.destroy()
		self.label.destroy()
		
	def createDisplay(self):
		isDisplayed = 1
		self.createButton()
		self.createLabel()
	
	def expand(self):
		self.button.config(image=self.minusImage)
		self.isExpanded = 1
	
	def contract(self):
		self.button.config(image=self.plusImage)
		self.isExpanded = 0
		
	def handl(self):
		self.button.destroy()
				
	
	
###########################################################

class UserDbViewer(Frame):
	def __init__(self, parent, root_path):
		Frame.__init__(self, parent, bg='white')
		self.nodesList = []
		self.molnir_level = parent
		self.root = root_path
		
		self.treeWinText = Pmw.ScrolledText(self, labelpos = 'n', usehullsize = 1, label_text='SCOP Nodes',
				hull_width = 100, hull_height = 20, text_padx = 1, text_pady = 1, text_wrap='none', text_cursor='arrow')
		self.treeWinText.pack(side=LEFT, expand=YES, fill=BOTH)
		self.pack(expand=YES, fill=BOTH)	
		
		self.viewButton = Button(self, text="View", state=DISABLED, height=1, command=self.displaySelected, width=5)
		self.viewButton.pack(side=RIGHT, anchor=W, padx=10)
		self.lines = 1

		#Parse the TopLevel Tree Structure Here		
		root_list = self.getDirectoryListing(self.root)
		if len(root_list) > 0:
			for item in root_list:
				self.nodesList.append(DirNode(self.root, item, self))
				indx = "%0.1f"%(self.lines)
				self.treeWinText.window_create(indx, window=self.nodesList[len(self.nodesList)-1])
				indx = "%0.1f"%(self.lines+0.1)
				self.treeWinText.insert(indx, '\n')
				self.lines = self.lines + 1
		
		self.currentSelection = self.nodesList[0]
		#	self.nodesList[len(self.nodesList)-1].pack(anchor=W)

		self.headerLabel = Label(self, text=" ")
		self.headerLabel.pack(side=BOTTOM, anchor=S)
		self.treeWinText.pack()
		self.treeWinText.configure(text_state='disabled')
		
		#self.config(width=400, height=500)

		self.pack(expand=YES, fill=BOTH)
	
	def getDirectoryListing(self, dest_path):
		objs = listdir(dest_path)
		directories = []
		files = []
		
		for item in objs:
			if path.isdir(dest_path+item):
				directories.append(item)
			else:
				if item.find('.pdb', len(item)-4) > -1:
					files.append(item)
	
		for item in files:
			directories.append(item)
			
		return directories
			
	def selectNode(self, nde):
		self.toggleSelection(nde)
		self.treeWinText.configure(text_state='normal')
		if nde.isExpanded == 1:
			pos = 0
			while nde.thisNode != self.nodesList[pos].thisNode:
				pos = pos + 1

			self.nodesList[pos].contract()
			
			pos = 0
			while len(self.nodesList) > pos:
#				verify = self.nodesList[pos].thisPath.split(nde.thisPath)				
#				if len(regex.split(nde.thisPath)) < len(regex.split(self.nodesList[pos].thisPath)):
#				if len(verify) == 2 and verify[1] == nde.thisNode+"\\":
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
			
		elif nde.isALeaf == 1:
			print "It is a leaf"
			
		else:
			pos = 0
			len_nodesList = len(self.nodesList)
			while pos < len_nodesList:
				item = self.nodesList[pos]
				if item.thisNode == nde.thisNode and item.thisPath == nde.thisPath:
					chld_pos = 1
					self.nodesList[pos].expand()
					selectedNode = self.nodesList[pos].thisPath+self.nodesList[pos].thisNode+"\\"
					if path.isdir(selectedNode):
						subdirs = self.getDirectoryListing(selectedNode)
						for child in subdirs:
							self.nodesList.insert(pos+chld_pos, DirNode(selectedNode, child, self))							
							indx = "%0.1f"%(pos+chld_pos+1)
							self.treeWinText.window_create(indx, window=self.nodesList[pos+chld_pos])
							indx = "%0.1f"%(pos+chld_pos+1+0.1)
							self.treeWinText.insert(indx, '\n')
							#self.nodesList[pos+1].pack(anchor=W)
							chld_pos = chld_pos + 1
						pos = len_nodesList
				pos = pos + 1
				
				
		print "Current Size of Nodes List: %d\n"%(len(self.nodesList))
		self.treeWinText.configure(text_state='disabled')
		
	
	def toggleSelection(self, nde):
#		self.currentSelection.label.config(bg='white', fg='black')
#		nde.label.config(bg='blue', fg='white')
		self.currentSelection.label.config(bg='white', fg='black')
		nde.label.config(bg='gray', fg='white')
		self.currentSelection = nde

	def displaySelected(self):
		pdbPath = self.currentSelection.thisPath+self.currentSelection.thisNode
				
		new_system = MolecularSystem.System(pdbPath, self, self.molnir_level)
		
		self.headerLabel.config(text=new_system.header)
		self.headerLabel.pack()
		new_system.color_ribbon_by_chain()
		self.molnir_level.system = new_system
		self.molnir_level.title(self.molnir_level.system.header)
		
		new_viewer = MolnirMolecularViewer(self.molnir_level, self.molnir_level.system)
		self.molnir_level.viewer = new_viewer
		new_cardframe = CardFrame(self.molnir_level, self.molnir_level.system)
		new_cardframe.pack(expand=NO, fill=X)
		new_window.cardframe = new_cardframe
		self.parent.system_windows.append(self.molnir_level)		
		


if __name__ == '__main__':

	root = Tk()
	frme = Frame(root, height=300, width=400)
	frme.pack(expand=NO)
	root_path = "C:\\CourseWork\\CS499\\MyTest\\SCOP_Pdbs"
	a = UserDbViewer(root, root_path)
	a.mainloop()        

