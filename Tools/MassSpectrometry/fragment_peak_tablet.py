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

class FPTablet(Frame):
    def __init__(self, parent, reaction_profiles=None, structure_filename=None, start_where=0):
        # open a pickled reaction_profiles object
        filename = './Tools/MassSpectrometry/protease_dict.pkl'
        if os.path.exists(filename):
            rp_file = open(filename, 'rb')
            self.reaction_profiles = pickle.load(protease_file)
            rp_file.close()
        else:
            print 'no protease file %s found'%(filename)

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
        # insert the molecular viewer
        new_window.pw = Pmw.PanedWidget(new_window,orient='vertical')
        new_window.viewer_pane  = new_window.pw.add('viewer', min=.1,max=.9,size=.67)
        new_window.viewer_pane.parent = new_window
        new_window.info_pane = new_window.pw.add('info', min=.1,max=.9,size=.33)
        new_window.info_pane.parent = new_window
        new_window.pw.pack(side=TOP,expand=1,fill='both')
        # open the system
        filename = askopenfilename(title = 'Select a structure file', defaultextension='.pdb', filetypes=[("pdb format", "*.pdb"),("system format", "*.sys"),("all files", "*")])
        new_window.system = MolecularSystem.System(filename)
        # open a molecular viewer
        for res in new_window.system.ProteinList[0].residues:
            if res.res_number in range(fragment_results[0]['site1'], fragment_results[0]['site2']):
                print self.parent.parent.experiment.get_modifying_reagent_specificity()
                if res.res_type1 in self.parent.parent.experiment.get_modifying_reagent_specificity():
                    print res.res_type1
                    res.vtk_arg_list['trace']['color'] = [0.85, 0.05, 0.05]
                res.vtk_arg_list['trace']['color'] = [0.85,0.85,0.05]
        new_window.viewer = MolecularViewer.MolecularViewer(new_window.viewer_pane, new_window.system)
        # adjust the view
