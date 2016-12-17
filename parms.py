import os
import pickle
import string


peakrec_fraction = 2
home_dir = 'C:\Users\Dude\Desktop\SPADE'

parms_dat = os.path.join(home_dir, 'parms.dat')

def get(query_parm=None):
    execfile(parms_dat, globals())
    if query_parm == None:
        return parms_list.keys()
    else:
        return parms_list[query_parm]

def set(query_parm, value):
    parms_file = open(parms_dat, 'r')
    lines = parms_file.readlines()
    for lind in range(len(lines)):
        line = string.strip(lines[lind])
        if query_parm in line:         # if already present, modify, else, append at the end
            lines[lind:lind+1] = '    \'%s\':%s,\n'%(query_parm, repr(value))
            break
    else:
        lineslen = len(lines)
        lines[lineslen-1:lineslen-1] = '    \'%s\':%s,\n'%(query_parm, repr(value))
    parms_file.close()
    parms_file = open(parms_dat, 'w')
    parms_file.writelines(lines)
    parms_file.close()

def what():
    execfile(parms_dat, globals())
    print parms_list.keys()
    