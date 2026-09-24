import os
import pickle

peakrec_fraction = 2
parms_dat = os.path.join(".", 'parms.dat')

def get(query_parm=None):
    with open(parms_dat, 'r') as f:
        exec(f.read(), globals())
    if query_parm == None:
        return parms_list.keys()
    else:
        return parms_list[query_parm]

def set(query_parm, value):
    parms_file = open(parms_dat, 'r')
    lines = parms_file.readlines()
    for lind in range(len(lines)):
        line = lines[lind].strip()
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
    with open(parms_dat, 'r') as f:
        exec(f.read(), globals())
    print(parms_list.keys())
    
