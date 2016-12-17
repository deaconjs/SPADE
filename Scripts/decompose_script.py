import sys
sys.path.append('Tools/DetectDomains')
import DetectDomains
import MolecularSystem
import Tkinter
import MolecularViewer
x = MolecularSystem.System(None)
x.load_pdb('4MDH.pdb')
window = Tkinter.Tk()
viewer = MolecularViewer.MolecularViewer(window, None)
viewer.load_system(x)
DetectDomains.auto_decompose(x.ProteinList[0], viewer)
window.mainloop()
