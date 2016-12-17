import sys
import time
sys.path.append('../..')
import MolecularSystem
from Superimpose import multi_superimpose

if __name__ == '__main__':
    start_time = time.clock()
    system_list = ['./test/1KAW_test_frag.pdb', './test/1L1OA_test_frag.pdb', './test/1L1OB_test_frag.pdb', './test/1L1OC_test_frag.pdb']
    systems = []
    for system in system_list:
        systems.append(MolecularSystem.System(None))
        systems[-1].load_pdb(system)
    for system in systems:
        system.ProteinList[0].fill_pseudo_sidechains(1)
        system.ProteinList[0].fill_neighbors_lists(0.35,15.0)
    multi_superimpose(systems)
    for ind in range(len(systems)):
        systems[ind].save_pdb('atoms%s.pdb'%(ind))
    end_time = time.clock()
    
    print 'program done in %5.3f seconds'%(end_time - start_time)
