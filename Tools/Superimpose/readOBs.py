import math
import sys
import string
sys.path.append("../..")
from MolecularSystem import System
x = System(None)
y = System(None)
z = System(None)
x.load_pdb('1KAW.pdb')
y.load_pdb('1L1OA.pdb')
z.load_pdb('1L1OB.pdb')
for prot in [x,y,z]:
    prot.ProteinList[0].fill_pseudo_sidechains(1)
    prot.ProteinList[0].fill_neighbors_lists(0.35,15.0)
x.res_list = [8, 14,31,32,33,34,35,57, 58, 59, 71, 77, 78, 79, 109]
y.res_list = [25,31,40,41,42,43,44,53, 54, 55, 63, 67, 68, 69, 84 ]
z.res_list = [74,80,91,92,93,94,95,104,105,106,120,126,127,128,146]

dsf = 0.15

do_replicate = 1
replicate_thresh = 0.05

shell_start = 7.0
shell_end   = 13.0
p_lo_hash = {}
p_hi_hash = {}
combinations = 0
beta_dist_sum = 0.0
beta_dist_cnt = 0.0
beta_dist_lzst = []

p_cnt = -1
for p in [x,y,z]:
    p_cnt += 1
    for rn1 in range(len(p.res_list)-3):
        b1 = p.ProteinList[0].residue_dict[p.res_list[rn1]].central_atom
        c1 = p.ProteinList[0].residue_dict[p.res_list[rn1]].pseudo_sidechain
        x1,y1,z1 = c1.x,c1.y,c1.z
        xb1,yb1,zb1 = b1.x,b1.y,b1.z
        for rn2 in range(rn1+1,len(p.res_list)-2):
            b2 = p.ProteinList[0].residue_dict[p.res_list[rn2]].pseudo_sidechain
            c2 = p.ProteinList[0].residue_dict[p.res_list[rn2]].central_atom
            d2 = c1.dist(c2)
            for rn3 in range(rn2+1,len(p.res_list)-1):
                b3 = p.ProteinList[0].residue_dict[p.res_list[rn3]].pseudo_sidechain
                c3 = p.ProteinList[0].residue_dict[p.res_list[rn3]].central_atom
                d3 = c1.dist(c3)
                for rn4 in range(rn3+1,len(p.res_list)):
                    b4 = p.ProteinList[0].residue_dict[p.res_list[rn4]].pseudo_sidechain
                    c4 = p.ProteinList[0].residue_dict[p.res_list[rn4]].central_atom
                    d4 = c1.dist(c4)
                    dist_list = [d2, d3, d4]
                    for d in dist_list:
                        if d<=shell_start or d>=shell_end:
                            break
                    else:
                        atom_list = [c2,c3,c4]
                        beta_list = [b2,b3,b4]
                        atom_num_list = [c2.atom_number, c3.atom_number, c4.atom_number]
                        sorted_list = [c2.atom_number, c3.atom_number, c4.atom_number]
                        sorted_list.sort()
                        f = [0,0,0]
                        for i in range(len(sorted_list)):
                           for j in range(len(dist_list)):
                                if atom_num_list[j] == sorted_list[i]:
                                    f[i] = j
                        xs = [atom_list[f[0]].x, atom_list[f[1]].x, atom_list[f[2]].x]
                        ys = [atom_list[f[0]].y, atom_list[f[1]].y, atom_list[f[2]].y]
                        zs = [atom_list[f[0]].z, atom_list[f[1]].z, atom_list[f[2]].z]

                        xbs = [beta_list[f[0]].x, beta_list[f[1]].x, beta_list[f[2]].x]
                        ybs = [beta_list[f[0]].y, beta_list[f[1]].y, beta_list[f[2]].y]
                        zbs = [beta_list[f[0]].z, beta_list[f[1]].z, beta_list[f[2]].z]

                        new_distance_list = [math.sqrt(((x1-   xs[0])**2) + ((y1-   ys[0])**2) + ((z1-   zs[0])**2)),
                                             math.sqrt(((x1-   xs[1])**2) + ((y1-   ys[1])**2) + ((z1-   zs[1])**2)),
                                             math.sqrt(((x1-   xs[2])**2) + ((y1-   ys[2])**2) + ((z1-   zs[2])**2)),
                                             math.sqrt(((xs[0]-xs[1])**2) + ((ys[0]-ys[1])**2) + ((zs[0]-zs[1])**2)),
                                             math.sqrt(((xs[0]-xs[2])**2) + ((ys[0]-ys[2])**2) + ((zs[0]-zs[2])**2)),
                                             math.sqrt(((xs[1]-xs[2])**2) + ((ys[1]-ys[2])**2) + ((zs[1]-zs[2])**2))]

                        bet_distance_list = [math.sqrt(((xb1-   xbs[0])**2) + ((yb1-   ybs[0])**2) + ((zb1-   zbs[0])**2)),
                                             math.sqrt(((xb1-   xbs[1])**2) + ((yb1-   ybs[1])**2) + ((zb1-   zbs[1])**2)),
                                             math.sqrt(((xb1-   xbs[2])**2) + ((yb1-   ybs[2])**2) + ((zb1-   zbs[2])**2)),
                                             math.sqrt(((xbs[0]-xbs[1])**2) + ((ybs[0]-ybs[1])**2) + ((zbs[0]-zbs[1])**2)),
                                             math.sqrt(((xbs[0]-xbs[2])**2) + ((ybs[0]-ybs[2])**2) + ((zbs[0]-zbs[2])**2)),
                                             math.sqrt(((xbs[1]-xbs[2])**2) + ((ybs[1]-ybs[2])**2) + ((zbs[1]-zbs[2])**2))]
                                
                        hires_distances  = [new_distance_list[0], new_distance_list[1], new_distance_list[2], new_distance_list[3], new_distance_list[4], new_distance_list[5]]
                        lowres_dl_bins   = [[],[],[],[],[],[]]
                        lowres_dlstrings = []
                        for i in range(len(new_distance_list)):
                            lowres_dl_bins[i].append(math.floor(dsf*new_distance_list[i]))
                            if do_replicate:
                                if (new_distance_list[i]*dsf)%1.0 <= replicate_thresh:              # if the distance is just over an integer change
                                    lowres_dl_bins[i].append((math.floor(dsf*new_distance_list[i]))-1)
                                elif (new_distance_list[i]*dsf)%1.0 >= (1.0-replicate_thresh):
                                    lowres_dl_bins[i].append((math.floor(dsf*new_distance_list[i]))+1)

                        if do_replicate:
                            for i0 in lowres_dl_bins[0]:
                                for i1 in lowres_dl_bins[1]:
                                    for i2 in lowres_dl_bins[2]:
                                        for i3 in lowres_dl_bins[3]:
                                            for i4 in lowres_dl_bins[4]:
                                                for i5 in lowres_dl_bins[5]:
                                                    lowres_dlstrings.append('%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_'%(i0,i1,i2,i3,i4,i5))
                        else:
                            lowres_dlstrings.append('%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_%2.0f_'%(lowres_dl_bins[0][0],
                                                                                            lowres_dl_bins[1][0],
                                                                                            lowres_dl_bins[2][0],
                                                                                            lowres_dl_bins[3][0],
                                                                                            lowres_dl_bins[4][0],
                                                                                            lowres_dl_bins[5][0]))
                        index_key = '%s %s %s %s'%(rn1,rn2,rn3,rn4)
                        try:
                            p_lo_hash[index_key]
                        except KeyError:
                            p_lo_hash[index_key] = [{'lowstr':lowres_dlstrings, 'betas':bet_distance_list, 'hilist':hires_distances}]
                        else:
                            p_lo_hash[index_key].append({'lowstr':lowres_dlstrings, 'betas':bet_distance_list, 'hilist':hires_distances})




keys = p_lo_hash.keys()
keys.sort()
good_count = 0

print '%s combinations'%(combinations)
print '%s keys'%(len(keys))

sum_alpha1 = 0.0
sum_beta1  = 0.0
sum_alpha2 = 0.0
sum_beta2  = 0.0
cnt_alpha = 0.0
values = []
distance_count = 0

for key in keys:
    if len(p_lo_hash[key]) == 3:
        bail = 0
        for s1 in range(len(p_lo_hash[key][0]['lowstr'])):
            for s2 in range(len(p_lo_hash[key][1]['lowstr'])):
                for s3 in range(len(p_lo_hash[key][2]['lowstr'])):
                    if p_lo_hash[key][0]['lowstr'][s1] == p_lo_hash[key][1]['lowstr'][s2]:
                         if p_lo_hash[key][1]['lowstr'][s2] == p_lo_hash[key][2]['lowstr'][s3]:
                            dist1, dist2, dist3, dist4 = 0.0, 0.0, 0.0, 0.0
                            # accumulate the squared distance
                            for d_ind in range(len(p_lo_hash[key][0]['hilist'])):
                                d1 = (p_lo_hash[key][0]['hilist'][d_ind] - p_lo_hash[key][1]['hilist'][d_ind])**2
                                d2 = (p_lo_hash[key][0]['hilist'][d_ind] - p_lo_hash[key][2]['hilist'][d_ind])**2
                                d3 = (p_lo_hash[key][0]['betas'][d_ind] -  p_lo_hash[key][1]['betas'][d_ind])**2
                                d4 = (p_lo_hash[key][0]['betas'][d_ind] -  p_lo_hash[key][2]['betas'][d_ind])**2
                                dist1      += d1
                                dist2      += d2
                                dist3      += d3
                                dist4      += d4
                                distance_count += 1
                            ln = len(p_lo_hash[key][0]['hilist'])
                            dist1,dist2,dist3,dist4 = math.sqrt(dist1/ln), math.sqrt(dist2/ln), math.sqrt(dist3/ln), math.sqrt(dist4/ln)
                            values.append([dist1, dist2, dist3, dist4])
                            sum_alpha1 += dist1
                            sum_alpha2 += dist2
                            sum_beta1  += dist3
                            sum_beta2  += dist4
                            bail = 1
                            print '\n*\n',
                            print '1 alpha %4.2f, beta %4.2f'%(dist1, dist3)
                            print '2 alpha %4.2f, beta %4.2f'%(dist2, dist4)
                            good_count += 1.0
                            break
                if bail:
                    break
            if bail:
                break
        key_tokens = string.split(key)
        key_tokens[0] = int(key_tokens[0])
        key_tokens[1] = int(key_tokens[1])
        key_tokens[2] = int(key_tokens[2])
        key_tokens[3] = int(key_tokens[3])
        print '%s\n[[%3s,%3s,%3s,%3s], [%3s,%3s,%3s,%3s], [%3s,%3s,%3s,%3s]]'%(p_lo_hash[key][0]['lowstr'], x.res_list[key_tokens[0]],x.res_list[key_tokens[1]],x.res_list[key_tokens[2]],x.res_list[key_tokens[3]],y.res_list[key_tokens[0]],y.res_list[key_tokens[1]],y.res_list[key_tokens[2]],y.res_list[key_tokens[3]],z.res_list[key_tokens[0]],z.res_list[key_tokens[1]],z.res_list[key_tokens[2]],z.res_list[key_tokens[3]])

# calculate the standard deviation of the different core analogies
sum = [0.0, 0.0, 0.0, 0.0]
for value in values:
    sum[0] += (value[0] - (sum_alpha1/good_count))**2
    sum[1] += (value[1] - (sum_alpha2/good_count))**2
    sum[2] += (value[2] - (sum_beta1/good_count))**2
    sum[3] += (value[3] - (sum_beta2/good_count))**2

for i in range(len(sum)):
    sum[i] /= (len(values)-1.0)

for i in range(len(sum)):
    sum[i] = math.sqrt(sum[i])

print '%s of %s good (%s)'%(good_count, len(keys), good_count/(len(keys)+0.0))
print 'averages - a1 %4.2f a2 %4.2f b1 %4.2f b2 %4.2f'%(sum_alpha1/good_count, sum_alpha2/good_count, sum_beta1/good_count, sum_beta2/good_count)
print 'deviatio -    %4.2f    %4.2f    %4.2f    %4.2f'%(sum[0], sum[1], sum[2], sum[3])







