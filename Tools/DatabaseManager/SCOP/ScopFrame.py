

from Tkinter import *
import Pmw

class ScopFrame(Frame):
	def __init__(self, parent=None):
		Frame.__init__(self, bg='white')
		self.toplevel = parent
		
		self.menuBar = Pmw.MenuBar(self, hull_relief=RAISED,	hull_borderwidth=1)
		self.menuBar.pack(fill=X, expand=YES, anchor=N)
		
		self.menuBar.addmenu('Profiles', 'Search Profiles')
		self.menuBar.addmenuitem('Profiles', 'command', label='Load Profiles')
		self.menuBar.addmenuitem('Profiles', 'command', label='Edit Profiles')

		self.menuBar.addmenu('Search', 'Search Entire SCOP Database')
		self.menuBar.addmenuitem('Search', 'command', label='Search SCOP')

		self.treeWinText = Text(self, cursor="arrow")		
		self.sbar = Scrollbar(self)
		self.sbar.config(command=self.treeWinText.yview)
		self.treeWinText.config(yscrollcommand=self.sbar.set)
		self.sbar.pack(side=LEFT, fill=Y) 
		self.treeWinText.config(height=15, width=60)
		self.treeWinText.pack(side=LEFT, expand=YES, fill=BOTH)

		self.lineageBox = Text(self, relief=FLAT, cursor = "arrow")
		self.lineageBox.pack(side=TOP, anchor=N, expand=YES, fill=BOTH)
		self.lineageBox.config(height=8, width=38, state=DISABLED)
		
		self.viewButton = Button(self, text="View", state=DISABLED, height=1, command=lambda: self.toplevel.displaySelected, width=5)
		self.viewButton.pack(side=LEFT, anchor=S, padx=10)
		
		self.pack(expand=YES, fill=BOTH)
		
		
if __name__ == '__main__':
	root = Tk()
	ScopFrame().mainloop()