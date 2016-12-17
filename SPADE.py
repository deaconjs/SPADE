# python imports
import os
import os.path
import sys
import string
import re
from Tkinter import *
from tkFileDialog import *
# dependency imports
sys.path.append('./Dependencies')
import vtk
#import vtkTkImageViewerWidget
import Pmw
# internal imports
import MolecularSystem
import MolecularViewer
import GUICards
sys.path.append(os.getcwd())
import SystemMemory
sys.path.append('./Tools/DatabaseManager')
sys.path.append('./Tools/ConservationTools')
import ConservationTools
import DatabaseDomainViewer
import DatabaseHandlers
sys.path.append('./Applications/Modification')
import ModificationViewer
sys.path.append('./Applications/AlignmentEditor')
import TreeSystem
sys.path.append('./Applications/Dynamics')
import DynamicsViewer

def load_application(parent, app, selected_title, active_window, active_window_key, chain=None):
    if app == 'Modification':
        reload(ModificationViewer)
        #modification_system = MolecularSystem.ModificationSystem(self.item)
        CM_viewer = ModificationViewer.CMViewer(systems[active_window_key], 
                                                active_window.application_pages[selected_title],
                                                active_window.viewer)
    elif app == 'AlignmentEditor':
        print 'loading alignment editor for chain %s'%(chain)
        reload(TreeSystem)
        window = active_window.application_pages[selected_title]
        viewer = active_window.viewer
        alignment_viewer = TreeSystem.TreeSystem(window, systems[active_window_key], None, 'editor', chain, viewer)
        alignment_viewer.pack(expand=YES, fill=BOTH)
            
    elif app == 'Dynamics':
        active_window.pw.configurepane('info', size=0.4)
        active_window.pw.configurepane('viewer', size=0.6)
        reload(DynamicsViewer)
        window = active_window.application_pages[selected_title]
        viewer = active_window.viewer
        active_system = systems[active_window_key]
        dynamics_viewer = DynamicsViewer.MDViewer(active_system, window, viewer)

class SPADEGUI(Frame):
    def __init__(self, parent=None):
        Frame.__init__(self, parent)
        self.master.title('SPADE')
        #self.run_mode = 'demo'      # just a thought
        self.run_mode = 'work'
        self.menu = MenuSystem(self)
        self.pack(expand=NO, fill=X)
        self.applicationbox = ApplicationBox(spade)
        self.objectbox = ObjectBox(self.run_mode)
        self.openCommandWindow()
        self.menu.fill_menu()
        self.system_windows = {}
        self.database_windows = {}

    def openCommandWindow(self):
        self.toplevel = Toplevel(spade)
        self.toplevel.title('Command Window')
        border=8
        self.toplevel.top_width=int(0.8*spade.winfo_screenwidth()-2*border)
        self.toplevel.top_height=int(0.182*spade.winfo_screenheight())
        geometry_string = "%dx%d%+d%+d"%(self.toplevel.top_width,
                                         self.toplevel.top_height,
                                         int(0.2*spade.winfo_screenwidth()+border),
                                         int(0.8*0.91*spade.winfo_screenheight())) # width,height,x-offset,y-offset
        self.toplevel.geometry(geometry_string)

        self.codebox  = CodeBox(self.toplevel)
        
class MenuSystem(Frame):
    def __init__(self, parent):
        """ create the main menu """
        Frame.__init__(self, parent)
        self.pack(side=TOP, fill=X, expand=NO)
        self.parent = parent
        # File operations
        file_button=Menubutton(self, text='File', underline=0)
        file_button.pack(side=LEFT, fill=X, expand=NO)
        file=Menu(file_button)
        file.add_command(label='Open a PDB', command=self.openNewPDB, underline=0)
        file.add_command(label='Open a System', command=self.openNewSystem, underline=0)
        file_button.config(menu=file)
        c_lambda = lambda: self.openDatabase("user_defined")
        file.add_command(label='Open a Database', command=c_lambda, underline=0)
        c_lambda = lambda: self.openDatabase("SCOP")
        file.add_command(label='Open the SCOP Database', command=c_lambda, underline=0)
        c_lambda = lambda: self.openDatabase("new")
        file.add_command(label='Create a Database', command=c_lambda, underline=0)
        c_lambda = lambda: self.openDatabase("transform")
        file.add_command(label='Transform a Database', command=c_lambda, underline=0)
        
    def fill_menu(self):
        # window controls
        file_button=Menubutton(self, text='Window', underline=0)
        file_button.pack(side=LEFT, fill=X, expand=NO)
        file=Menu(file_button)
        file.add_command(label='Command Window', command=self.parent.openCommandWindow, underline=0)
        file_button.config(menu=file)
        #file.add_command(label='Object Window', command=self.parent.openObjectWindow, underline=0)
        #file_button.config(menu=file)
        # script controls
        #file_button=Menubutton(self, text='Scripts', underline=0)
        #file_button.pack(side=LEFT, fill=X, expand=NO)
        #file=Menu(file_button)
        #file.add_command(label='Execute', command=self.parent.codebox.executeCurrentScript, underline=0)
        #file_button.config(menu=file)
        #file.add_command(label='New Script', command=self.parent.codebox.addNewScript, underline=0)
        #file_button.config(menu=file)
        #file.add_command(label='Close Current Script', command=self.parent.codebox.closeCurrentScript, underline=0)
        #file_button.config(menu=file)
        #file.add_command(label='Save Current Script', command=self.parent.codebox.saveCurrentScript, underline=0)
        #file_button.config(menu=file)
        
    def openDatabase(self, type):
        if type == "new" or type == "transform":
            # a database transformation draws from one source and supplies a user-defined database.
            # Create a window for choosing a target database name and a database source.
            tnsf_win = Toplevel(spade)
            geometry_string = "%dx%d%+d%+d" %(250,100,400,200) # width,height,x-offset,y-offset
            tnsf_win.geometry(geometry_string)
            # add an entry form for the target database name
            tnsf_win.target_entry = Pmw.EntryField(tnsf_win, labelpos = 'w', label_text = 'Target Database:', validate = None, value='default_database')
            tnsf_win.target_entry.pack(side=TOP,fill='x', expand=1, padx=10, pady=5)
            # add a button box for the different sources
            tnsf_win.source_buttonbox = Pmw.ButtonBox(tnsf_win, labelpos = 'w',label_text='Data Sources:',frame_borderwidth=2, orient='vertical', frame_relief = 'groove')
            c_lambda = lambda: DatabaseHandlers.download_from_rcsb(tnsf_win.target_entry.getvalue())
            tnsf_win.source_buttonbox.add('PDBs from the RCSB', command = c_lambda)
            c_lambda = lambda: DatabaseHandlers.transform_from_pdbs(tnsf_win.target_entry.getvalue())
            tnsf_win.source_buttonbox.add('From local PDB files', command = c_lambda)
            tnsf_win.source_buttonbox.pack(side=TOP,fill='x', expand=1, padx=10, pady=5)
        elif type == "SCOP":
            self.parent.objectbox.database_listbox.insert('end',"SCOP")
            # reset the selection
            for i in range(self.parent.objectbox.database_listbox.size()):
                self.parent.objectbox.database_listbox.select_clear(i)
            self.parent.objectbox.database_listbox.select_set(len(databases)-1)
        elif type == "user_defined":
            """ load a new database """
            root_path = os.path.normpath("./Databases")
            db_path = os.path.normpath(askdirectory(initialdir=root_path, title="Select the Database Directory", mustexist=1))
            if len(db_path) > 0:
                soln = re.search('.*SPADE.*',db_path)
                if not soln:
                    print "The directory must be within SPADE' Databases subdirectory\n"
                    self.openNewDatabase()  # if it fails, recurse through this function
                    return
                # now remove everything up to 'Databases'. This assumption allows short and meaningful names.
                db_token = os.path.normpath(os.path.abspath("./"))
                rep_db_path = string.split(db_path, db_token, 1)
                databases.append(rep_db_path[1])        # slice out the slash that gets left behind
                self.parent.objectbox.database_listbox.insert('end',"%s"%(rep_db_path[1]))
                # reset the selection
                for i in range(self.parent.objectbox.database_listbox.size()):
                    self.parent.objectbox.database_listbox.select_clear(i)
                self.parent.objectbox.database_listbox.select_set(len(databases)-1)
                
    def openNewSystem(self):
        # these two functions should be condensed and create a new system directory from an arbitrary pdb or sps file
        pass
            
    def openNewPDB(self):
        pass
            
    def closeSystem(self, which_window):
        """ destroy the which_window'th system window and delete it from the list
            -- superseded by objectbox.close_window()""" 
        spade.ui.system_windows[which_window].destroy()
        if len(spade.ui.system_windows) > which_window+1:
              spade.ui.system_windows = spade.ui.system_windows[:which_window] + spade.ui.system_windows[which_windows+1:]
        else:
              spade.ui.system_windows = spade.ui.system_windows[:which_window]

class ApplicationBox:
    def __init__(self, parent):
        """ SPADE' ApplicationBox is a simple construct for launching
            applications into active windows
        """
        self.parent = parent
        self.toplevel = spade
        screenwidth=spade.top_width
        screenheight=spade.top_height
        self.applications_listbox = Pmw.ScrolledListBox(self.toplevel,
                                                        listbox_exportselection=0,
                                                        labelpos='nw',
                                                        label_text='applications',
                                                        usehullsize=1,
                                                        selectioncommand=self.select_target_system)
        self.applications_listbox.pack(side=TOP,expand=1,fill='both', anchor=N)
        self.applications_available = os.listdir('./Applications')
        for application in self.applications_available:
            self.applications_listbox.insert('end', '%s'%(application))
        self.active_window_key = None
            
    def select_target_system(self, event=None):
        # if more than one system_window is open, ask which one to apply to
        # can't use a buttonbox here because it disallows use of '_' keys in titles
        active_window = None
        if len(spade.ui.system_windows.keys()) + len(spade.ui.database_windows.keys()) > 1:
            buttonlist = ['Ok']
            self.dialog = Pmw.Dialog(self.parent,
                                     buttons = buttonlist,
                                     buttonbox_orient = 'vertical',
                                     defaultbutton = buttonlist[0],
                                     title = 'Select a System Window',
                                     command = self.return_active_window)
            self.dialog.withdraw()
            # Add some contents to the dialog.
            w = Label(self.dialog.interior(),text = 'Select a Window')
            w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
            syss = spade.ui.system_windows.keys()
            dbs  = spade.ui.database_windows.keys()
            self.dialog.listbox = Listbox(self.dialog.interior(), exportselection=0)
            for sys in syss:
                self.dialog.listbox.insert(END, sys)
            for db in dbs:
                self.dialog.listbox.insert(END, db)
            self.dialog.listbox.pack(expand=1, fill='both')
            self.dialog.show()
        elif len(spade.ui.system_windows.keys()) == 1 or len(spade.ui.database_windows.keys()) == 1:
            buttonlist = ['Ok', 'Cancel']
            self.launch_decision = 0
            self.dialog = Pmw.Dialog(self.parent,
                                     buttons = buttonlist,
                                     buttonbox_orient = 'horizontal',
                                     defaultbutton = buttonlist[0],
                                     title = 'Confirm Application Launch',
                                     command = self.decide_on_launch)
            self.dialog.withdraw()
            # Add some contents to the dialog.
            if len(spade.ui.system_windows.keys()) == 1:
                w = Label(self.dialog.interior(),text = 'Launching an application on %s'%(spade.ui.system_windows.keys()[0]))
            else:
                w = Label(self.dialog.interior(),text = 'Launching an application on %s'%(spade.ui.database_windows.keys()[0]))
            w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
            self.dialog.show()
        else:
            print 'open a system or database first'
                
            
    def decide_on_launch(self, result):
        if result == 'Ok':
            if len(spade.ui.system_windows.keys()) == 1:
                self.active_window_key = spade.ui.system_windows.keys()[0]
            elif len(spade.ui.database_windows.keys()) == 1:
                self.active_window_key = spade.ui.database_windows.keys()[0]
            self.launchApplication()
            self.dialog.destroy()
        else:
            self.dialog.destroy()
            
    def return_active_window(self, result):
        i = self.dialog.listbox.curselection()
        if len(i) > 0:
            idx = int(i[0])
            if idx <= len(spade.system_windows.keys()):
                self.active_window_key = spade.ui.system_windows.keys()[idx]
            else:   # the listbox curselection content order is system keys first database keys second, so
                self.active_window_key = spade.ui.database_windows.keys()[idx-len(spade.ui.system_windows.keys())]
        self.dialog.destroy()
        self.launchApplication()
        
    def launchApplication(self):
        """ for now, write code here or write a function and access it from here. Keep It Simple """
        if self.active_window_key == None or len(self.applications_listbox.curselection())==0:
            return
        x = self.applications_listbox.curselection()[0]
        y = self.applications_listbox.get(int(x))
        selected_title = string.split(y)[0]
        cntr = 2
        keys = string.split(selected_title, ' ')
        if keys[0] == 'AlignmentEditor':
            active_window = spade.ui.system_windows[self.active_window_key]
            for pchain in systems[self.active_window_key].ProteinList:
                selected_title = 'AlignmentEditor %s'%(pchain.chain_name)
                active_window.application_pages[selected_title] = active_window.application_notebook.add(selected_title)
                active_window.pw.configurepane('info', size=0.5)    # default sizes
                active_window.pw.configurepane('viewer', size=0.5)
                active_window.application_notebook.selectpage(selected_title)
                # open a new notebook page
                load_application(self, keys[0], selected_title, active_window, self.active_window_key, pchain.chain_name)
                active_window.pw.setnaturalsize()
            active_window.viewer.display('volume', 1)
        elif keys[0] == 'Procast':
            active_window = spade.ui.database_windows[self.active_window_key]
            # preprocess here by performing the alignment from a selected database.
            #   This will eventually be replaced with a new keys option in this function, that will opt
            #   that Procast is called to interact with a previously opened, interactive (or at least
            #   re-runnable) multiply positioned alignment Application
            # collect all systems from the database
            selected_index = self.parent.ui.objectbox.database_listbox.curselection()
            if len(selected_index) == 0:
                print "none selected"
                return
            selected_title = string.strip(self.parent.ui.objectbox.database_listbox.get(selected_index)[1:])
            print 'title %s'%(selected_title)
            db_dir = './Databases/' + selected_title
            
            # open a new System window
            self.parent.ui.objectbox.launchWindow('system', 'empty')
            # call the Procast UI construct to open the superimposed systems
            
        elif self.active_window_key in spade.ui.system_windows.keys(): # all other System Window Applications
            active_window = spade.ui.system_windows[self.active_window_key]
            while selected_title in active_window.application_pages.keys():
                if cntr == 2:
                    selected_title = '%s %s'%(selected_title, '2')
                else:
                    selected_title = '%s %s'%(selected_title[:-2], '%s'%(cntr))
                cntr += 1
            active_window.application_pages[selected_title] = active_window.application_notebook.add(selected_title)
            active_window.pw.configurepane('info', size=0.5)    # default sizes
            active_window.pw.configurepane('viewer', size=0.5)
            active_window.application_notebook.selectpage(selected_title)
            # open a new notebook page
            load_application(self, keys[0], selected_title, active_window, self.active_window_key)
            active_window.pw.setnaturalsize()

class ObjectBox:
    def __init__(self, run_mode):
        """SPADE' ObjectBox has three panels. A list of open system objects and databases are on the top two,
        and the bottom contains a MolecularViewer. Click the 'O' to open a viewer. Click 'X' to close
        the selected item and its viewer.
        """
        self.run_mode = run_mode
        if run_mode == 'demo':
            self.systems_directory = 'Examples'
        elif run_mode == 'work':
            self.systems_directory = 'Systems'
        self.databases_directory = 'Databases'
        self.toplevel = spade
        # First create the pane widget
        self.pw = Pmw.PanedWidget(self.toplevel,orient='vertical')
        self.object_pane  = self.pw.add('systems', min=.1,max=.9,size=.35)
        self.database_pane = self.pw.add('databases', min=.1,max=.9,size=.3)
        self.viewer_pane   = self.pw.add('Viewer',  min=.1,max=.9)
        self.pw.pack(side=TOP,expand=1,fill='both')
        # make two scrolled listboxes, for containing objects and databases
        self.system_buttons = Pmw.ButtonBox(self.object_pane,
                                     orient='vertical')
        self.system_buttons.add('O',command=self.launchSystemWindow)
        self.system_buttons.pack(side=RIGHT, expand=0,fill=Y)
        self.system_listbox = Pmw.ScrolledListBox(self.object_pane,
                                                  labelpos='nw',
                                                  label_text='systems',
                                                  usehullsize=1,
                                                  selectioncommand=self.loadMiniView)
        self.system_listbox.pack(side=RIGHT,expand=1,fill='both')
        self.database_buttons = Pmw.ButtonBox(self.database_pane,
                                              orient='vertical')
        self.database_buttons.add('O',command=self.launchDatabaseWindow)
        #self.database_buttons.add('X',command=self.closeDatabase)
        self.database_buttons.pack(side=RIGHT, expand=0,fill=Y)
        self.database_listbox = Pmw.ScrolledListBox(self.database_pane,
                                                    labelpos='nw',
                                                    label_text='databases',
                                                    usehullsize=1)
        self.database_listbox.pack(side=RIGHT,expand=1,fill='both')
        w = int(round(spade.top_width,0))
        h = int(round(spade.top_height/4.0,0))
        #self.miniviewer = vtkTkImageViewerWidget.vtkTkImageViewerWidget(self.viewer_pane, width=w, height=h)
        #self.miniviewer.pack()
        # systems_not_loaded is a list of directories in the Systems directory
        self.systems_not_loaded = []
        path = os.path.join(os.getcwd(), self.systems_directory)
        self.systems_not_loaded = os.listdir(path)
        # system_filenames gets all of these filenames
        for sys_file in self.systems_not_loaded:
            for system_key in systems.keys():
                if system_key == sys_file:
                    self.system_listbox.insert('end',"O %s"%(sys_file))
                    break
            else:
                self.system_listbox.insert('end',".  %s"%(sys_file))

        self.databases_not_loaded = []
        path = os.path.join(os.getcwd(), self.databases_directory)
        self.databases_not_loaded = os.listdir(path)
        # database_filenames gets all of these filenames
        for database_file in self.databases_not_loaded:
            for database_key in databases.keys():
                if database_key == database_file:
                    self.database_listbox.insert('end',"O %s"%(database_file))
                    break
            else:
                self.database_listbox.insert('end',".  %s"%(database_file))

    def loadMiniView(self, event='None'):
        """ Load system and view 
        """
        selected_title = self.system_listbox.get(self.system_listbox.curselection())
        # if the system is not previously loaded, load it
        keys = string.split(selected_title, ' ', 1)
        keys[1] = string.strip(keys[1])
        internal_location = keys[1] #os.path.normpath(os.path.join(self.systems_directory, keys[1]))
        print 'internal loc %s'%(internal_location)
        location = None
        open_type = None
        if keys[0] == '.':                                                          # if the system at the given path is not loaded
            descriptor = string.split(internal_location, os.path.normpath('/'))[-1] # get the filename from the directory
            if string.find(descriptor, '.') == -1:                                  # if no descriptor, its a directory, search the subdirectory for system
                full_path = os.path.normpath(os.path.join(os.getcwd(), self.systems_directory, keys[1])) # create a full path string for the directory
                d_list = os.listdir(full_path)                                      # get the directory's listing to look for systems
                for descriptor2 in d_list:                                          # first look for pickled system (.sps) files    
                    filenametokens = string.split(descriptor2, '.')                   # figure out if its a file or directory
                    if len(filenametokens) > 0:                                       # skip directories
                        if filenametokens[-1] == 'sps':                               # if its a sps file
                            location = os.path.normpath(os.path.join(full_path, descriptor2))
                            open_type = 'sps'                                       # and the type to load
                            break
                else:                                                               # if no sps files, search for .pdbs instead
                    for descriptor2 in d_list:                                      # first look for pickled system (.sps) files    
                        filenametokens = string.split(descriptor2, '.')               # figure out if its a file or directory
                        if len(filenametokens) > 0:                                   # skip directories
                            if filenametokens[-1] in ['pdb','PDB']:                           # if its a sps file
                                location = os.path.normpath(os.path.join(full_path, descriptor2))
                                open_type = 'pdb'                                   # and the type to load
                                break
                    else:                                                       # else fail
                        print 'No .sps or .pdb files found in the directory'
                        if len(d_list) == 0:
                            print 'In fact, directory %s is empty'%(full_path)
                        return
            else:                                                               # if it is a file
                split_filename = string.split(descriptor, '.')
                if split_filename[-1] == 'sps':
                    location = os.path.normpath(os.path.join(os.getcwd(), internal_location))
                    open_type = 'sps'
                elif split_filename[-1] == 'pdb':
                    location = os.path.normpath(os.path.join(os.getcwd(), internal_location))
                    open_type = 'pdb'

            if internal_location not in systems.keys():
                new_sys = MolecularSystem.System(location)
                systems[internal_location] = new_sys               # systems may be retrieved by location or by index
                print 'system available at systems[%s]'%(internal_location)

            # set the selection of the listbox
            for i in range(self.system_listbox.size()):
                self.system_listbox.select_clear(i)
            for i in range(self.system_listbox.size()):
                if string.strip(string.split(self.system_listbox.get(i), ' ', 1)[1]) == keys[1]:
                    self.system_listbox.select_set(i)
                    break

        system = systems[internal_location]
        #self.loadIcon(system)
    
    def loadIcon(self, system):
        if os.path.exists(system.get_filename_by_extension('.jpg')):
            reader = vtk.vtkJPEGReader()
            reader.SetFileName(system.get_filename_by_extension('jpg'))
            reader.Update()
            extent = reader.GetDataExtent()
            resizer = vtk.vtkImageResample()
            resizer.SetInterpolationModeToCubic()
            resizer.SetInput(reader.GetOutput())
            scale = min((spade.top_width+0.0)/extent[1], (spade.top_height/2.0)/extent[3])
            resizer.SetAxisMagnificationFactor(0,scale)
            resizer.SetAxisMagnificationFactor(1,scale)
            self.miniviewer.GetImageViewer().SetColorWindow(200)
            self.miniviewer.GetImageViewer().SetColorLevel(100)
            self.miniviewer.GetImageViewer().SetInput(resizer.GetOutput())
            
        else:
            print 'no %s image available'%(system.get_filename_by_extension('.jpg'))
        self.miniviewer.Render()
    
    def load3DMiniView(self, system):
        graphics_memory = SystemMemory.GraphicsMemory(system)
        selection_memory = SystemMemory.SelectionMemory(system)
        SystemMemory.QuickRenderVisitor(system)
        SystemMemory.ObjectSelectionVisitor(system)
        self.miniviewer.loadSystem(system)
        graphics_memory.restore_system()
        selection_memory.restore_system()
        
    def closeDatabase(self):
        """ close the selected database """
        # which database is it?
        selected_index = self.database_listbox.curselection()
        if len(selected_index)==0:
            return
        selected_title = self.database_listbox.get(selected_index)
        # destroy the system window (if present) and remove it from the list
        for db in spade.ui.database_windows:
            if db.title() == selected_title:
                spade.ui.database_windows.remove(db)
                db.destroy()
        # maybe the next bit should be implemented, but i'm not so sure
        # remove the system from the databases list
        #for dbkey in databases.keys:
        #    db = databases[dbkey]
        #    if db.header == selected_title:
        #        save = db
        #databases.remove(save)
        #self.database_listbox.delete(selected_index)
        
    def closeSystem(self):
        """ close the selected system """
        # which system is it?
        selected_index = self.system_listbox.curselection()
        if len(selected_index)==0:
            return
        selected_title = self.system_listbox.get(selected_index)
        # destroy the system window and remove it from the list
        for sys in spade.ui.system_windows.keys():
            if sys.title() == selected_title:
                del spade.ui.system_windows[sys]
                sys.destroy()
        # maybe the next bit should be implemented, but i'm not so sure
        # remove the system from the systems list
        #for syskey in systems.keys:
        #    sys = systems[syskey]
        #    if sys.header == selected_title:
        #        save = sys
        #    del systems.remove(save)
        #self.system_listbox.delete(selected_index)
        #self.miniviewer.closeSystem()

    def launchDatabaseWindow(self):
        self.launchWindow('database')
        return

    def launchSystemWindow(self):
        """ create a window to view the selected system object """
        self.launchWindow('system')
        return
        
    def launchWindow(self, type, input_system=None):
        # the 'empty' keyword can be used in input_system to produce an empty system window for manipulation
        new_window = Toplevel()
        # adjust screen geometry
        border = 8
        screenwidth=5*spade.top_width-border
        screenheight=spade.top_height
        geometry_string = "%dx%d%+d%+d" %(int(0.8*screenwidth),
                                          int(0.9*screenheight),
                                          int(0.2*screenwidth+border),
                                          int(1)) # width,height,x-offset,y-offset
        new_window.geometry(geometry_string)

        if type == 'database':
            # which database is it?
            selected_index = self.database_listbox.curselection()
            if len(selected_index)==0:
                print "none selected"
                return
            selected_title = string.strip(self.database_listbox.get(selected_index)[1:])
            # if its a user_defined database
            if selected_title != "SCOP":
                # create the true path
                pathname = os.path.join('./Databases', selected_title)
        elif type == 'system':
            # which system is it?
            if input_system == 'empty':
                selected_index = None
                selected_title = ""
                system_key     = None
                system_index   = None
            else:
                selected_index = self.system_listbox.curselection()
                if len(selected_index)==0:
                    return
                selected_title = string.strip(self.system_listbox.get(selected_index)[1:])
                i = 0
                system_key = None
                for sys_key in systems.keys():
                    if string.find(systems[sys_key].filename, selected_title) != -1:
                        system_key   = sys_key
                        system_index = i
                        break
                    i=i+1

        if type == 'system':
            if input_system == 'empty':
                sys_obj = None
            elif input_system == None:
                new_window.title('%s -- %s'%(sys_key, systems[sys_key].header))
                sys_obj = systems[sys_key]
            else:
                sys_obj = input_system
        else:
            new_window.title("%s Database"%(selected_title))
            sys_obj = None
            
        c_lambda = lambda win=new_window: self.close_window(win)
        new_window.protocol("WM_DELETE_WINDOW",c_lambda)
        # store it in the window
        new_window.system = sys_obj
        new_window.parent = spade
        
        if type == 'database':
            # add a few panes
            new_window.panes = Pmw.PanedWidget(new_window)
            topPane    = new_window.panes.add('top', min=10)
            new_window.bottompane = new_window.panes.add('bottom', min=20)
            new_window.panes.pack(expand=1, fill='both')
            # add the Db viewer
            if selected_title == "SCOP":
                new_window.dbViewer = ScopDomainViewer.ScopViewer(new_window.panes, new_window.panes, 'None')
            else:
                print 'opening db'
                print new_window.panes
                print pathname
                print new_window
                new_window.dbViewer = DatabaseDomainViewer.UserDbViewer(new_window.panes, pathname, new_window)
            new_window.dbViewer.pack(expand=1, fill=BOTH)
            # add an empty system
            new_window.bottompane.system = MolecularSystem.System(None)
            # add a molecular viewer
            target_pane = new_window.bottompane
        elif type == 'system':
            target_pane = new_window
            
        new_window.pw = Pmw.PanedWidget(target_pane,orient='vertical')
        new_window.viewer_pane  = new_window.pw.add('viewer', min=.1,max=.9,size=.75)
        new_window.viewer_pane.parent = new_window
        new_window.info_pane = new_window.pw.add('info', min=.1,max=.9,size=.25)
        new_window.info_pane.parent = new_window
        new_window.pw.pack(side=TOP,expand=1,fill='both')
        # open a molecular viewer
        print 'here 1'
        new_window.viewer = MolecularViewer.MolecularViewer(new_window.viewer_pane, new_window.system)
        print 'here 2'
        new_window.application_notebook = Pmw.NoteBook(new_window.info_pane, borderwidth=1, tabpos='n')
        new_window.application_notebook.pack(expand=YES, fill=BOTH)
        new_window.application_pages = {}
        new_window.application_pages['Info'] = new_window.application_notebook.add('Info')
        # new_window.application_notebook.setnaturalsize()
        new_cardframe = GUICards.CardFrame(new_window.application_pages[new_window.application_pages.keys()[0]], new_window.system, new_window)
        new_cardframe.pack(expand=NO, fill=X)
        new_window.cardframe = new_cardframe
        if type == 'database':
            spade.ui.database_windows[selected_title] = new_window
        elif type == 'system':
            spade.ui.system_windows[selected_title] = new_window
        new_window.pw.setnaturalsize()

        
    def close_window(self, win):
        """ 1) closes system windows
        2) removes them from the system_window list
        3) reinitializes molecule and graphics selection states
        """
        for key in spade.ui.system_windows.keys():
            if spade.ui.system_windows[key] == win:
                del spade.ui.system_windows[key]
                break
        else:
            for key in spade.ui.database_windows.keys():
                if spade.ui.database_windows[key] == win:
                    del spade.ui.database_windows[key]
        win.destroy()
        
class CodeBox:
    def __init__(self, parent_window, alias_system=None, alias_viewer=None):
        """SPADE' CodeBox has two panels, where output is presented on the left or bottom, and the right 
        or top contains a textbox where code can be simply written and evaluated, and where script files can 
        be managed. A PMW Notepad allows multiple files to be edited simultaneously.
        """
        print 'Codebox Next and Prev commands dont yet capture w/out mouse interaction'
        self.alias_system = alias_system         # used when codebox is attached to a MolecularViewer
        self.alias_viewer = alias_viewer
        self.toplevel = parent_window
        # controls label
        labelfont = ('times', 8)
        self.controls_label = Label(self.toplevel, text='NEW:ctrl+N   EXECUTE:ctrl+X   SAVE:ctrl+S   CLOSE:ctrl+shift+Q   NEXT:ctrl+right arrow   PREV:ctrl+left arrow', font=labelfont)
        self.controls_label.pack(side=TOP,fill='x', expand=0, padx=0,pady=0)
        # First create the paned widget
        self.pw = Pmw.PanedWidget(self.toplevel,orient='horizontal',
                                  hull_borderwidth = 1,
                                  hull_relief = 'sunken',
                                  hull_height=100)
        self.toplevel.bind("<Control-x>", self.executeCurrentScript)
        self.toplevel.bind("<Control-n>", self.addNewScript)
        self.toplevel.bind("<Control-Shift-q>", self.closeCurrentScript)
        self.toplevel.bind("<Control-s>", self.saveCurrentScript)
        self.toplevel.bind("<Control-Right>", self.nextScript)
        self.toplevel.bind("<Control-Left>", self.lastScript)
        self.output_pane  = self.pw.add('Output', min=.1,max=.9,size=.5)
        self.input_pane   = self.pw.add('Input',  min=.1,max=.9,size=.5)
        self.pw.pack(side=BOTTOM, anchor=S,expand=0,fill='x')
        # The left gets a scrolled text box
        self.output_textbox = Pmw.ScrolledText(self.output_pane, borderframe=5,text_wrap='none')
        self.output_textbox.configure(text_state='disabled')
        self.output_textbox.pack(expand=1, fill='both', padx=2, pady=2)
        # The right gets a Pmw Notebook
        self.file_notebook = Pmw.NoteBook(self.input_pane, borderwidth=1,pagemargin=2)
        self.file_notebook.pack(fill='both', expand = 1, padx = 2, pady = 2)
        self.file_pages = {}
        self.addNewScript()
        # create a default starting page
        #self.pw.setnaturalsize()
        self.pw.updatelayout()
        
    def destroy(self):
        self.controls_label.destroy()
        self.pw.destroy()

    def addNewScript(self, event=None):
        """ make a new filename of the format tmp#.mpy, where # is some value not present in the current list
         of opened mpy scripts."""
        names = self.file_notebook.pagenames()
        print names
        for i in range(0,10):
            file_name = 'tmp%d.mpy'%(i)
            if file_name not in names:
                break
        else:
            print "Too many files open (>10)\n"
        print file_name
        self.file_pages[file_name] = self.file_notebook.add(file_name)
        self.file_pages[file_name].text = Pmw.ScrolledText(self.file_pages[file_name], borderframe=5,text_wrap='none')
        self.file_pages[file_name].text.pack()

    def nextScript(self, event):
        self.file_notebook.nextpage()
        active_page = self.file_pages[self.file_notebook.getcurselection()]
        active_page.text.focus_force()
        active_page.text.tag_add(SEL, 1.0, END)

    def lastScript(self, event):
        self.file_notebook.previouspage()
        active_page = self.file_pages[self.file_notebook.getcurselection()]
        active_page.text.focus_force()
        active_page.text.tag_add(SEL, 1.0, END)

    def execute_script_callback(self, event):
        self.executeCurrentScript()

    def executeCurrentScript(self, event=None):
        """ execute the currently selected script """
        selected_title = self.file_notebook.getcurselection()
        self.saveCurrentScript()
        if self.alias_system:
            globals = {'system':self.alias_system, 'viewer':self.alias_viewer}
        else:
            globals = {'systems':systems, 'spade':spade, 'system_windows':spade.ui.system_windows}
        locals  = {}
        exec open(os.path.join('Scripts', selected_title)).read() in globals, locals
        
    def saveCurrentScript(self, event=None):
        """ save the currently selected script """
        selected_title = self.file_notebook.getcurselection()
        text = self.file_pages[selected_title].text.get(1.0,'end')
        self.executing_file = open('Scripts/%s'%(selected_title), 'w')
        self.executing_file.writelines(text)
        self.executing_file.close()

    def closeCurrentScript(self, event=None):
        """ close the currently selected script """
        selected_title = self.file_notebook.getcurselection()
        self.file_notebook.delete(selected_title)
        
    def write(self, text):
        """ write text to the outbox of the codebox. This functions is present to override stdout. Use
        python's normal print() command instead."""
        self.output_textbox.configure(text_state='normal')
        self.output_textbox.insert(END, text)
        self.output_textbox.configure(text_state='disabled')
        self.output_textbox.see(END)
        self.output_textbox.update_idletasks()
        
if __name__ == '__main__':
    
    # begin by initializing the Tk and manipulating its geometry
    spade = Tk()
    spade.title('objectbox')
    screenwidth=spade.winfo_screenwidth()
    screenheight=spade.winfo_screenheight()
    #toplevel_geometry_string = "%dx%d%+d%+d" %(192,int((screenheight/2)-((screenheight/10)+15)),1,int(screenheight/2)) # width,height,x-offset,y-offset
    spade.top_width=int(0.2*screenwidth)
    spade.top_height=int(0.91*screenheight)
    geometry_string = "%dx%d%+d%+d" %(spade.top_width,spade.top_height,1,1) # width,height,x-offset,y-offset
    spade.geometry(geometry_string)

    systems   = {}              # stores the MolecularSystem objects
    databases = {}              # just stores the directory (a string) of the database
    # redirect print statements to the codebox out box

    spade.ui = SPADEGUI()
    #sys.stdout = spade.ui.codebox
    #sys.stderr = spade.ui.codebox
    spade.ui.mainloop()
    
