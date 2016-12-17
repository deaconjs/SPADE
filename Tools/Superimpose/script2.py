import sys
import math
sys.path.append('../..')
from MolecularSystem import System
x = System(None)
y = System(None)
z = System(None)
x.load_pdb('1KAW.pdb')
y.load_pdb('1L1OA.pdb')
z.load_pdb('1L1OB.pdb')
x.res_list = [8,14,32,34,57,59,71,77,79,101,109]
y.res_list = [25,31,41,43,53,55,63,67,69,79,84]
z.res_list = [74,80,92,94,104,106,120,126,128,141,146]
dsf = 0.2
for res_num1 in range(len(x.res_list)):
    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    for res_num2 in range(0,len(x.res_list)):
        #print x.ProteinList[0].residue_dict[x.res_list[res_num2]].res_number, x.ProteinList[0].residue_dict[x.res_list[res_num2]].central_atom.atom_number
        #print y.ProteinList[0].residue_dict[y.res_list[res_num2]].res_number, y.ProteinList[0].residue_dict[y.res_list[res_num2]].central_atom.atom_number
        #print z.ProteinList[0].residue_dict[z.res_list[res_num2]].res_number, z.ProteinList[0].residue_dict[z.res_list[res_num2]].central_atom.atom_number
        #continue
        dist1 = x.ProteinList[0].residue_dict[x.res_list[res_num1]].central_atom.dist(x.ProteinList[0].residue_dict[x.res_list[res_num2]].central_atom)
        dist2 = y.ProteinList[0].residue_dict[y.res_list[res_num1]].central_atom.dist(y.ProteinList[0].residue_dict[y.res_list[res_num2]].central_atom)
        dist3 = z.ProteinList[0].residue_dict[z.res_list[res_num1]].central_atom.dist(z.ProteinList[0].residue_dict[z.res_list[res_num2]].central_atom)
        sum1 += dist1
        sum2 += dist2
        sum3 += dist3
        #if res_num2 > res_num1:
        str = '%5s %5s %5.2f %5s %5s %5.2f %5s %5s %5.2f %2s %2s %2s'%(x.res_list[res_num1], x.res_list[res_num2], dist1*dsf, 
                                                                    y.res_list[res_num1], y.res_list[res_num2], dist2*dsf, 
                                                                    z.res_list[res_num1], z.res_list[res_num2], dist3*dsf,
                                                                    int(math.floor(dist1*dsf)),
                                                                    int(math.floor(dist2*dsf)),
                                                                    int(math.floor(dist3*dsf))),
        addon = ""
        if int(math.floor(dist1*dsf)) == int(math.floor(dist2*dsf)):
            if int(math.floor(dist2*dsf)) == int(math.floor(dist3*dsf)):
              addon = " *"
        print '%s%s'%(str[0], addon)
    #print 'ave dist = %5.3f, %5.3f, %5.3f'%(sum1/(len(x.res_list)+0.0),sum2/(len(x.res_list)+0.0),sum3/(len(x.res_list)+0.0))



