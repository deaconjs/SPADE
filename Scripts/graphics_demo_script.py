import MolecularViewer
import MolecularSystem
import Tkinter
import sys
sys.path.append('./Dependencies')
import Pmw

system = MolecularSystem.System(None)      # open a molecular structure
system.load_pdb('otcA3.pdb')                        # www.rcsb.org format
window = Tkinter.Tk()
viewer = MolecularViewer.MolecularViewer(window, system)

chain = system.ProteinList[0]
chain.vtk_arg_list['trace']['representation'] = 'tube'
chain.vtk_arg_list['trace']['width'] = 0.8
chain.vtk_arg_list['trace']['splines'] = 10
chain.vtk_arg_list['trace']['specular'] = 0.8
chain.vtk_arg_list['trace']['specular_power'] = 50
for resnum in range(len(chain.residues)):
    v = resnum/float(len(chain.residues))
    chain.residues[resnum].vtk_arg_list['trace']['color'] = [v,v,1.0]

viewer.update_view()

window.mainloop()

