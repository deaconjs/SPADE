# this program is a component of RAVE that handles visualization of
# sets of Reaction_Profiles.
import sys
import os
sys.path.append(os.path.abspath('../../Tools/Plotting'))
import PlottingTools
import Tkinter
import MolecularSystem
import MolecularViewer
import Pmw
import pickle

class Viewer(Tkinter.Frame):
    def __init__(self, parent, reaction_profiles=None, start_where=None):
        print reaction_profiles
        self.reaction_profiles = reaction_profiles
        self.current_profile_index = 0
        self.current_fragment_index = 0
        if start_where:    # start_where is a select peak profile. search for it.
            i = 0
            for rp in self.reaction_profiles:
                if rp.get_peak == start_where:
                    self.current_profile_index = i
                i += 1
        
        self.parent = parent
        self.experiment = parent.parent.parent.experiment
        # create a toplevel
        new_window = Tkinter.Toplevel()
        # adjust screen geometry
        #border = 8
        #screenwidth=5*spade.top_width-border
        #screenheight=spade.top_height
        #geometry_string = "%dx%d%+d%+d" %(int(0.8*screenwidth),
        #                                  int(0.9*screenheight),
        #                                  int(0.2*screenwidth+border),
        #                                  int(1)) # width,height,x-offset,y-offset
        #new_window.geometry(geometry_string)
        new_window.pw = Pmw.PanedWidget(new_window,orient='vertical')
        new_window.viewer_pane  = new_window.pw.add('viewer', min=.1,max=.9,size=.67)
        new_window.viewer_pane.parent = new_window
        new_window.info_pane = new_window.pw.add('info', min=.1,max=.9,size=.33)
        new_window.info_pane.parent = new_window
        new_window.pw.pack(side="top",expand=1,fill='both')
        # open the system
        str_filename = self.experiment.get_structure_filename()
        if len(str_filename) == 0:
            str_filename = askopenfilename(title = 'Select a structure file', defaultextension='.pdb', filetypes=[("pdb format", "*.pdb"),("system format", "*.sys"),("all files", "*")])
        new_window.system = MolecularSystem.System(str_filename)
        # adjust the view and open a molecular viewer
        for res in new_window.system.ProteinList[0].residues:
            if res.res_number in range(fragment_results[0]['site1'], fragment_results[0]['site2']):
                if res.res_type1 in self.experiment.get_modifying_reagent_specificity():
                    res.vtk_arg_list['trace']['color'] = [0.85, 0.05, 0.05]
                res.vtk_arg_list['trace']['color'] = [0.85,0.85,0.05]
        new_window.viewer = MolecularViewer.MolecularViewer(new_window.viewer_pane, new_window.system)
        # add some plots
        
         
        

