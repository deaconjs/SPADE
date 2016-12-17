import sys
sys.path.append("/mnt/work/Research/Bliss") # Hard coded path
import MolecularSystem

ms = MolecularSystem.System ("test_pdbs/1c7k.pdb")

"""
Test the protonation and hbond functionality
"""
ms.create_hbonds(None,None,True)
ms.print_hbonds()
ms.save_pdb("output.pdb")
