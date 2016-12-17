from MolecularSystem import System
import sys
import os
sys.path.append('./Tools/Aligner')
import SequenceAligner

het_insert_dir = './Databases/heterodimers'
#hom_insert_dir = './Databases/homodimers'
het_extras_dir  = './Databases/old sweeney hetero/het'
#hom_extras_dir  = './Databases/dimer_pdb'

het_insert_dirs = os.listdir(het_insert_dir)
#hom_insert_dirs = os.listdir(hom_insert_dir)
het_extras_pdbs = os.listdir(het_extras_dir)
#hom_extras_pdbs = os.listdir(hom_extras_dir)

print 'hetero files'
for extra_name in het_extras_pdbs:
    line = 'scanning %s'%(extra_name) 
    print '\n' + line
    outf = open('sweeney_insertion_out.txt', 'a')
    outf.write(line + '\n')
    extra_pdb_name = os.path.join(het_extras_dir, extra_name, extra_name) + '.pdb'
    print 'opening %s'%(extra_pdb_name)
    extra_system = System(extra_pdb_name)
    for insert_dir_name in het_insert_dirs:
        insert_pdb_name = os.path.join(het_insert_dir, insert_dir_name, insert_dir_name) + '.pdb'
        print 'opening %s'%(insert_pdb_name)
        insert_system = System(insert_pdb_name)
        for extra_pchain in extra_system.ProteinList:
            for insert_pchain in insert_system.ProteinList:
                # do a sequence comparison
                a = SequenceAligner.SequenceAligner(0, 'global')
                a.add_target(extra_pchain)
                a.add_template(insert_pchain)
                pid = a.align_sequences()
                if pid > 0.25:
                    line = '%5.3f percent identity between %s %s %s and %s %s %s'%(pid, insert_dir_name, insert_pchain.chain_name, len(insert_pchain.residues), extra_name, extra_pchain.chain_name, len(extra_pchain.residues))
                    print line
                    outf.write(line+'\n')
    outf.close()
"""
print '\n\n homo files'
for extra_name in hom_extras_pdbs:
    print 'scanning %s\n'%(extra_name)
    extra_pdb_name = os.path.join(hom_extras_dir, extra_name)
    extra_system = System(extra_pdb_name)
    for insert_dir_name in hom_insert_dirs:
        insert_pdb_name = os.path.join(hom_insert_dir, insert_dir_name, insert_dir_name) + '.pdb'
        insert_system = System(insert_pdb_name)
        for extra_pchain in extra_system.ProteinList:
            for insert_pchain in insert_system.ProteinList:
                # do a sequence comparison
                a = SequenceAligner.SequenceAligner(0, 'global')
                a.add_target(extra_pchain)
                a.add_template(insert_pchain)
                pid = a.align_sequences()
                if pid > 0.25:
                    line = '%5.3f percent identity between %s %s %s and %s %s %s\n'%(pid, insert_dir_name, insert_pchain.chain_name, len(insert_pchain.residues), extra_name, extra_pchain.chain_name, len(extra_pchain.residues))
                    print line
                    outf.write(line)
"""
print 'done'
                 