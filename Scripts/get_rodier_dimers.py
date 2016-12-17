import os
from ftplib import FTP
import sys
sys.path.append('./Applications')
from AlignmentEditor import Sequencer
import MolecularSystem

lists = {}
types = ['heterodimers', 'homodimers']
# homodimers (116)
# replaced 1tcl with 1tc1 (from Bahadur 2003)
# 1alo >> 1vlb
# monomers: 1a3c, 1amk, 1b68, 1b69, 1bif, 1buo, 1csh, 1ctt,
#           1czj, 1e98, 1gvp, 1jhg, 1jsg, 1kpf, 1mor, 1nox,
#           1opy, 1rpo, 1tox #? # , 1uby, 1utg, 2ilk, 2tct,
#           2tgi, 3grs, 3ssi, 5csm, ?? 9wga ??, 
# multimers: 1cg2, 1fro, 1hjr, 1hss, 1vlb, 4cha, 

lists['homodimers'] = ['12as', '1a3c', '1a4i', '1a4u', '1aa7', '1ad3', '1ade', '1afw', '1ajs', '1vlb',
                       '1aor', '1aq6', '1auo', '1b3a', '1b5e', '1b67', '1b68', '1b69', '1b70',
                       '1bbh', '1bd0', '1bif', '1biq', '1bis', '1bjw', '1bkp', '1bmd', '1brw', '1bsl',
                       '1bsr', '1buo', '1bxg', '1bxk', '1cdc', '1cg2', '1chm', '1cmb', '1cnz', '1coz',
                       '1csh', '1ctt', '1cvu', '1czj', '1daa', '1dor', '1dpg', '1dqs', '1dxg', '1e98',
                       '1ebh', '1f13', '1fip', '1fro', '1gvp', '1hjr', '1hss', '1hxp', '1icw', '1imb',
                       '1isa', '1ivy', '1jhg', '1jsg', '1kba', '1kpf', '1m6p', '1mkb', '1mor', '1nox',
                       '1nse', '1nsy', '1oac', '1opy', '1pgt', '1qfh', '1qhi', '1qr2', '1r2f', '1reg',
                       '1rpo', '1ses', '1slt', '1smn', '1sox', '1tc1', '1tox', '1trk', '1uby', '1utg',
                       '1vfr', '1vok', '1wtl', '1xso', '2arc', '2ccy', '2hdh', '2ilk', '2lig', '2mcg',
                       '2nac', '2ohx', '2spc', '2sqc', '2tct', '2tgi', '3dap', '3grs', '3sdh', '3ssi',
                       '4cha', '5csm', '5rub', '8prk', '9wga']

# heterodimers (46)
lists['heterodimers'] = ['1a2k', '1acb', '1ak4', '1avw', '1brs', '1bth', '1cbw', '1cho', '1cse', '1dan',
                         '1dfj', '1dhk', '1dvf', '1efn', '1efu', '1fin', '1fle', '1gg2', '1got', '1gua',
                         '1hia', '1hwg', '1igc', '1kb5', '1mct', '1mel', '1mlc', '1nca', '1nmb', '1osp',
                         '1ppf', '1stf', '1tgs', '1tx4', '1vfb', '1ycs', '1ydr', '2kai', '2pcc', '2ptc',
                         '2sic', '2sni', '2trc', '3sgb', '3tpi', '4htc']


# download the pdb files (remove statement included)
"""
ftp = FTP('ftp.rcsb.org', 'anonymous', '')   # connect to host, default port
ftp.cwd("pub/pdb/data/structures/all/pdb")
for type in types:
    for pdb in lists[type]:
        dirloc = './Databases' + '/' + type + '/' + pdb
        fileloc = dirloc + '/' + pdb + '.pdb.Z'
        if os.path.exists(fileloc):
            print 'skipping file %s'%(fileloc)
            #os.remove(fileloc)
        else:
            request = 'pdb' + pdb + '.ent.Z'
            print 'putting %s at %s'%(request, fileloc)
            ftp.retrbinary('RETR ' + request, open(fileloc, 'wb').write)
ftp.quit()
"""

# construct initial sequence alignments
for type in types:
    for pdb in lists[type]:
        dirloc = './Databases' + '/' + type + '/' + pdb
        pdbloc = dirloc + '/' + pdb + '.pdb'
        fstloc = dirloc + '/' + pdb + '.fst'
        if os.path.exists(fstloc) and os.path.getsize(fstloc) != 0:
            print 'skipping file %s'%(pdbloc)
        else:
            system = MolecularSystem.System(pdbloc)
            for pchain in system.ProteinList:
                sequence = pchain.get_sequence()
                x = Sequencer.Sequence_Through_Clustalw( sequence, system )
                print 'output %s '%( x)
                break
        break
    break

