import string
infile = open('./conservation_analysis.out', 'r')
lines = infile.readlines()
types = ['normalized_0D_conservation', 'normalized_noactsit_0D_conservation',
         'normalized_1D_conservation', 'normalized_3D_conservation',
         'normalized_ms3D_conservation', 'normalized_noactsit_ms3D_conservation',
         'normalized_nobadloop_ms3D_conservation']

data = {}
for type in types:
    data[type] = []
    for line in lines:
        if type in line:
            ind = int(string.find(line, ' percent of the time'))
            pwin = float(string.strip(line[ind-4: ind]))
            data[type].append(pwin)

easy_indices, easy_thresh = [], 100.0
hard_indices, hard_thresh = [], 30.0
for i in range(len(data['normalized_0D_conservation'])):
    sum = 0.0
    for type in types:
        sum += data[type][i]
    avg = sum/float(len(types))
    if avg >= easy_thresh:
        easy_indices.append(i)
    elif avg <= hard_thresh:
        hard_indices.append(i)

print '%s easy indices, %s hard'%(len(easy_indices), len(hard_indices))        

for type in types:
    sum = 0.0
    count = 0.0
    i = 0
    for datum in data[type]:
        if i not in easy_indices and i not in hard_indices:
            sum += datum
            count += 1
        i += 1
    print 'average %s, count %s - %s '%(sum/float(count), count, type)
infile.close()
print '>>>'
