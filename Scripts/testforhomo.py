import sys
import os
sys.path.append('./Tools/Aligner')
import SequenceAligner

from MolecularSystem import System

list = os.listdir('.')
for pdb in list:
    sys = System(pdb)
    print pdb
    for i in range(len(sys.ProteinList)-1):
        for j in range(i+1,len(sys.ProteinList)):
            a = SequenceAligner.SequenceAligner(0, 'global')
            a.add_target(sys.ProteinList[i])
            a.add_template(sys.ProteinList[j])
            pid = a.align_sequences()
            line = '%5.3f percent identity %s %s\n'%(pid, sys.ProteinList[i].chain_name, sys.ProteinList[j].chain_name)
            print line
