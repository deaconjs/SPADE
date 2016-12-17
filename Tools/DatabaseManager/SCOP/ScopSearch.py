

from Tkinter import *
from Bio import SCOP
from os import *
import Pmw
import re
import parms
from ScopResults import *

import MolecularSystem
from GUICards import *
from BlissMolecularViewer import BlissMolecularViewer

class SearchFrame(Frame):
	#Arguments-
	#parent-it is the Toplevel() window
	#viewer-it is the ScopViewer object
	def __init__(self, parent, viewer):
		Frame.__init__(self, parent, bg='white')
		self.pack(expand=YES, fill=BOTH)
		self.scopViewer = viewer
		self.scopStruct = viewer.scopDataStruct
		self.tier1_frame = Frame(self, bg='white')
		self.currentExpr = re.compile('')
		
		self.toplevel = parent
				
		geometry_string = "%dx%d%+d%+d" %(700,500,150,150) # width,height,x-offset,y-offset
		self.toplevel.geometry(geometry_string)
		self.toplevel.title('SCOP Search')		
		
		self.keywords=[]
		self.removedKwds=[]
		self.includeRegex = []
		self.excludeRegex = []
		self.entireClaList = {}
		self.commsList = {}
		self.buildClaList()   #Build a hash table with sunids as kes from the classification.txt file
		self.currentNodeClaList = []		
		self.loadKeywords()
		self.types_state = []
		self.listOfresultsSunid = []
		self.saveResultsin = parms.get('saveResultsin')

		#########
		self.add_txt_frame = Frame(self.tier1_frame, bg='white')
		self.enter_txt_label = Label(self.add_txt_frame, bg='white', text='Enter Keyword:', justify=LEFT)
		self.enter_txt_label.pack(anchor=W, pady=5)
		self.kwd_entry = Entry(self.add_txt_frame, width=15)
		self.kwd_entry.pack()
		self.kwd_entry.focus()
		self.kwd_entry.bind('<Return>', (lambda event: self.addKeywordToList()))
		self.add_kwd_button = Button(self.add_txt_frame, text='Add Keyword', command=self.addKeywordToList)
		self.add_kwd_button.pack(pady=5)

		self.exclude_txt_label = Label(self.add_txt_frame, bg='white', text='Exclude Keyword:', justify=LEFT)
		self.exclude_txt_label.pack(anchor=W, pady=5)
		self.exclude_kwd_entry = Entry(self.add_txt_frame, width=15)
		self.exclude_kwd_entry.pack()
		self.exclude_kwd_entry.bind('<Return>', (lambda event: self.addExcludeKeywordToList()))
		self.exclude_kwd_button = Button(self.add_txt_frame, text='Add Keyword', command=self.addExcludeKeywordToList)
		self.exclude_kwd_button.pack(pady=5)
		
		self.add_txt_frame.pack(side=LEFT, anchor=N, padx = 25, pady=75)
		
		##########
		self.keyword_lstbox_frame = Frame(self.tier1_frame, bg='white')
		self.kwd_listbox_label = Label(self.keyword_lstbox_frame, text='Keywords List', bg='white')
		self.kwd_listbox_label.pack(pady=5)

		self.kwd_listbox = Pmw.ScrolledListBox(self.keyword_lstbox_frame, listbox_selectmode=SINGLE, items=self.keywords,
			listbox_height=15, vscrollmode='static', usehullsize=1, hull_width=200, hull_height=250)
		self.kwd_listbox.config(background='white')

		self.kwd_listbox.pack(side=LEFT)
		self.keyword_lstbox_frame.pack(side=LEFT, anchor=N)
				
		##########
		self.add_remove_frame = Frame(self.tier1_frame, bg='white')
		self.add_from_remove = Button(self.add_remove_frame, text='<-- Add', width=10, command=self.addFromRemoveList)
		self.add_from_remove.pack(padx=15, pady=15)
		self.remove_from_keyword = Button(self.add_remove_frame, text='Remove -->', width=10, command=self.removeFromKeywordList)
		self.remove_from_keyword.pack(padx=15)		
		self.add_remove_frame.pack(side=LEFT, anchor=N, pady=75)
		
		##########
		self.remove_lstbox_frame = Frame(self.tier1_frame, bg='white')
		self.rm_listbox_label = Label(self.remove_lstbox_frame, text='Removed Keywords List', bg='white')
		self.rm_listbox_label.pack(pady=5)
		self.rm_listbox = Pmw.ScrolledListBox(self.remove_lstbox_frame, listbox_selectmode=SINGLE, items = self.removedKwds,
			listbox_height=15, vscrollmode='static', usehullsize=1, hull_width=200, hull_height=250)
		self.rm_listbox.config(background='white')
		self.rm_listbox.pack(side=LEFT)
		self.remove_lstbox_frame.pack(side=LEFT, anchor=N)
		
		self.tier1_frame.pack()
		
		self.type = parms.get('classifications_fullname')
		self.indx = parms.get('classifications')
		self.classi_index = parms.get('nameHierarchy')
				
		self.tier2_frame = Frame(self, bg='white')
		self.toInlcludeSublevels = IntVar(self)
		self.include_sublevels_checkbox = Checkbutton(self.tier2_frame, text='Search sublevels of Hits', bg='white', variable=self.toInlcludeSublevels)
		self.include_sublevels_checkbox.pack(side=LEFT)
		self.searchComments = IntVar(self)
		self.searchComments_checkbox = Checkbutton(self.tier2_frame, text='Search Comments', bg='white', variable=self.searchComments)
		self.searchComments_checkbox.pack(side=LEFT)
		self.tier2_frame.pack()	
		
		self.tier3_frame = Frame(self, bg='white')
		self.chkBox_frame_label = Label(self.tier3_frame, text='Refine Search by selecting the TYPES', bg='white', justify=LEFT)
		self.chkBox_frame_label.pack(side=TOP, anchor=W, fill=X)
		
		self.chk_vars=[]
		for sel in self.indx:
			var = IntVar(self)
			chk = Checkbutton(self.tier3_frame, text=self.type[sel], bg='white', variable=var)
			chk.pack(side=LEFT, anchor=W)
			if sel == 'cl': #Donot search the Class Level of the hierarchy, because of a BUG while 
					#displaying the corresponding hit in SearchResults
				chk.forget()
			else:
				chk.select()
			self.chk_vars.append(var)
		
		self.tier3_frame.pack(pady=30)
		self.tier4_frame = Frame(self, bg='white')
		self.ok_button = Button(self.tier4_frame, text='OK', command=self.okbutton)
		self.ok_button.pack(side=LEFT, padx=100)
		self.cancel_button = Button(self.tier4_frame, text='Cancel', command=self.cancelbutton)
		self.cancel_button.pack(side=LEFT, padx=100)
		self.tier4_frame.pack()
	
	#A 'keywords.ini' file exists in which the keywords are saved for later use.
	#This ini file is loaded in the following function.
	def loadKeywords(self):
		if path.isfile(parms.get('keywordsFile')):
			wordsFile = file(parms.get('keywordsFile'), 'r')
			strg = wordsFile.readline()
			while strg.find('#keywords') == -1:
				strg = wordsFile.readline()
			strg = wordsFile.readline()
			while strg.find('#removed keywords') == -1:
				for words in strg.split(','):
					if len(words.strip()) > 0:
						self.keywords.append(words.strip())
				strg = wordsFile.readline()		
					
			strg = wordsFile.readline()
			while strg:
				for words in strg.split(','):
					if len(words.strip()) > 0:
						self.removedKwds.append(words.strip())
				strg = wordsFile.readline()
	
	#After the search is finished, update the 'keywords.ini' file	
	def updateKeywords(self):
		wordsFile = file(parms.get('keywordsFile'), 'w')
		wordsFile.write('#keywords\n')
		for word in self.keywords:
			if self.keywords.index(word) < (len(self.keywords) - 1):				
				wordsFile.write(word+', ')
			else:
				wordsFile.write(word)
		
		wordsFile.write('\n#removed keywords\n')
		for word in self.removedKwds:
			if self.removedKwds.index(word) < (len(self.removedKwds) - 1):				
				wordsFile.write(word+', ')
			else:
				wordsFile.write(word)
	
	def addKeywordToList(self):
		if len(self.kwd_entry.get()) > 0:
			self.keywords.append('(+)'+self.kwd_entry.get())		
			if self.keywords.index('(+)'+self.kwd_entry.get()) == (len(self.keywords) - 1): 
				self.kwd_listbox.setlist(self.keywords)
			else:
				self.keywords.pop()
			
			self.kwd_entry.delete(0, len('(+)'+self.kwd_entry.get()))

	def addExcludeKeywordToList(self):
		if len(self.exclude_kwd_entry.get()) > 0:
			self.keywords.append('(-)'+self.exclude_kwd_entry.get())		
			if self.keywords.index('(-)'+self.exclude_kwd_entry.get()) == (len(self.keywords) - 1): 
				self.kwd_listbox.setlist(self.keywords)
			else:
				self.keywords.pop()
			
			self.exclude_kwd_entry.delete(0, len('(-)'+self.exclude_kwd_entry.get()))
		
	def removeFromKeywordList(self):
		if len(self.kwd_listbox.getvalue()) != 0:
			sel = self.kwd_listbox.getvalue()[0]
			self.keywords.pop(self.keywords.index(sel))
			self.kwd_listbox.setlist(self.keywords)
			
			self.removedKwds.append(sel)
			if self.removedKwds.index(sel) == (len(self.removedKwds) - 1):
				
				self.rm_listbox.setlist(self.removedKwds)
			else:
				self.removedKwds.pop()
		
	def addFromRemoveList(self):
		if len(self.rm_listbox.getvalue()) != 0:
			sel = self.rm_listbox.getvalue()[0]
			self.removedKwds.pop(self.removedKwds.index(sel))
			self.rm_listbox.setlist(self.removedKwds)
			
			self.keywords.append(sel)
			if self.keywords.index(sel) == (len(self.keywords) - 1):
				
				self.kwd_listbox.setlist(self.keywords)
			else:
				self.keywords.pop()
		
	def okbutton(self):	
		if len(self.keywords) > 0:
			self.types_state = []
			#Get the state of Check Boxes
			for var in self.chk_vars:
				self.types_state.append(var.get())
			
			#Delete the previous search results
			tempFiles = listdir(self.saveResultsin)
			if len(tempFiles) > 0:
				dialog = Pmw.MessageDialog(self, title = 'Temp Results Directory',
					defaultbutton = 0,
					buttons = ('OK', ),
					message_text = 'Deleting Results of previous Search!')
				dialog.activate()
	
			for rfile in tempFiles:
				os.remove(self.saveResultsin+rfile)
			
			#Parse the comments file and build a hash
			allComms = file(parms.get('scopComments'), 'r').readlines()
			for line in allComms:
				self.commsList[line[:5]] = line[8:]				
			
			#Build lists of regular expressions 
			for item in self.keywords:
				if item.find('(+)') == 0:
					self.includeRegex.append(re.compile(item[3:], re.IGNORECASE))
				elif item.find('(-)') == 0:
					self.excludeRegex.append(re.compile(item[3:], re.IGNORECASE))
					
			for kwrd in self.keywords:
				self.currentNodeClaList = []
				self.findKeyword(self.scopStruct.root)
				
			self.updateKeywords()
			print "KEY WORDS ARE:"
			print self.keywords
			if len(self.listOfresultsSunid) > 0:
				#If there are any search results, display them in ResultsFrame
				self.destroy()
				self.toplevel.title("SCOP Search Results")
				self.toplevel.config(bg='white')

				pw = Pmw.PanedWidget(self.toplevel)	
				pw.add('top', min=10)
				pw.add('bottom', min = 10)
				
				geometry_string = "%dx%d%+d%+d" %(700,800,50,30) # width,height,x-offset,y-offset
				self.toplevel.geometry(geometry_string)
				
				pw.pane('top').configure(background='white')
				a = ResultsFrame(self.saveResultsin, pw, self.toplevel, self.scopViewer).pack(side=TOP)
				pw.pack(expand=1, fill=BOTH) 

				bottomPane = pw.pane('bottom')
				bottomPane.system = MolecularSystem.System(self, bottomPane)
				bottomPane.system.load_pdb("")
				bottomPane.system.color_ribbon_by_chain()
				bottomPane.viewer = MolnirMolecularViewer(bottomPane, bottomPane.system)
				bottomPane.viewer.loadSystem(bottomPane.system)
				bottomPane.cardframe = CardFrame(bottomPane, bottomPane.system)
				bottomPane.cardframe.pack(expand=NO, fill=X)				
				self.scopViewer.topMenuSystem.parent.scopViewer_windows.append(self.toplevel)
			else:
				dialog = Pmw.MessageDialog(self, title = 'Search Results',
					defaultbutton = 0,
					buttons = ('OK', ), command='self.destroy()',
					message_text = 'Search returned 0 results!')
				dialog.activate()

			
		else:
			print "No Keywords Entered"			
	
	def cancelbutton(self):
		self.toplevel.destroy()
	
	#For all matching keywords, a separate SCOP Structure is created (by creating the parseable files)	
	def findKeyword(self, node):
		if  self.keywordNodeMatch(node) == 1 and len(node.type) > 0 and self.types_state[self.classi_index[node.type]-1] == 1:
			self.listOfresultsSunid.append(node.sunid)
			print "Search Found! - "+str(node.description)
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

			if self.toInlcludeSublevels.get() == 1:
				for ndes in node.children:
					self.currentNodeClaList = []
					self.findKeyword(ndes)
		else:
			for ndes in node.children:
				self.currentNodeClaList = []
				self.findKeyword(ndes)

	#Implements the actual search functionality. An 'and' operation is performed to find SCOP Nodes 
	#which match all the keywords entered
	
	def keywordNodeMatch(self, node):
		
		commStr = ""
		if self.searchComments.get() == 1 and self.commsList.has_key(node.sunid):
			commStr = self.commsList[node.sunid]
			
		for expr in self.excludeRegex:
			if expr.search(node.description) is not None or expr.search(commStr) is not None:
				return (0)
		
		for expr in self.includeRegex:
			if expr.search(node.description) is None and expr.search(commStr) is None:
				return (0)
		
		return (1)

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
			#else:
			#	print "Unable to index the line: "+line


if __name__ == '__main__':
	main_window = Tk()	
	clasi = file(parms.get('scopClassification'), 'r')
	descr = file(parms.get('scopDescription'), 'r')
	heira = file(parms.get('scopHierarchy'), 'r')
	
	scopStruct = SCOP.Scop(clasi, descr, heira)

	a = SearchFrame(main_window, scopStruct.root)
	a.mainloop()


