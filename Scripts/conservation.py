# python imports
import sys
import os
import os.path
sys.path.append(os.getcwd())
import copy
import math
import time
import string
import random
from Tkinter import *
# tool imports
import SPADE
sys.path.append('./Tools/DetectDomains')
import DetectDomains
sys.path.append('./Tools/Selection')
import SystemSelectionDialog
sys.path.append('./Tools/Aligner')
import SequenceAligner
sys.path.append('./Tools/ConservationTools')
import ConservationTools
# handle command line the lame way
homedir = '.'
from sys import argv         
if len(argv) == 2:
    homedir = 'C:\Documents and Settings\Deacon\Desktop\SPADE'

sys.path.append(homedir)

# dependency imports
from vtk import *
import vtkpython
import vtkTkRenderWidget
import tkFileDialog
#from vtk.tk.vtkTkRenderWindowInteractor import *
from MolecularComponents.classFutamuraHash import FutamuraHash 
import vtk

# internal imports
sys.path.append(os.getcwd())
sys.path.append(os.path.join(homedir, './Dependencies'))
import Pmw
import MolecularSystem
import parms
import SystemMemory
import MolecularViewer
#sys.path.append(os.path.join(homedir, './Tools/SequenceFetcher'))
#import SequenceFetcher

verbose = 0


if __name__ == '__main__':
    pdb_list = ['1a00', '1a0o', '1a22', '1a50', '1a6d', '1a97', '1atn', '1cf1', '1cqi', '1ct9', '1eq2', '1fm2', '1ira', '1kut', '1ewy',
                '1a9n', '1aa1', '1ab8', '1acm', '1ad1', '1adu', 
                '1agr', '1all', '1am4', '1ay4', '1b3q', '1b8g', 
                '1c1b', '1cb7', '1cc1', '1cd9', '1cqz', '1dce', 
                '1djs', '1e2h', '1e7w', '1e51', '1e9h', '1efp', 
                '1eg9', '1ep1', '1evu', '1fpp', '1g5h', '1gmn', 
                '1h0m', '1h1l', '1ha3', '1hfw', '1hxp', '1hyh', 
                '1i7q', '1ib1', '1ipj', '1j94', '1kkm', '1ko8', 
                '1kq4', '1ksh', '1lx7', '1qkr', '1qs0', '1scf', 
                '1sct', '1sky', '1tub', '1wap']
    
    #  

    for pdb in pdb_list:
        outfile = open('./conservation_out.txt', 'a')
        sys.stdout = outfile
        sys.stderr = outfile
        pdb_loc = './Systems/'+pdb+'/'
        pdb_file= pdb_loc+pdb+'.pdb'
        new_system = MolecularSystem.System(pdb_file)
        ConservationTools.calculate_conservation(new_system, 'sidechain_asa', 0)
        new_system.calculate_interface_significance('sidechain_asa')
        outfile.close()

