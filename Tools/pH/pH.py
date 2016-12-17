import math
pKs = {'A': [2.4,  9.9],
       'C': [1.9, 10.8,   8.3],
       'D': [2.0,  9.9,   3.9],
       'E': [2.1,  9.5,   4.1],
       'F': [2.2,  9.2],
       'G': [2.4, 9.8],
       'H': [1.8,  9.2,   6.0],
       'I': [2.3,  9.8],
       'K': [2.2,  9.2,  10.8],
       'L': [2.3,  9.7],
       'M': [2.1,  9.3],
       'N': [2.1,  8.8],
       'P': [2.0, 10.6],
       'Q': [2.2,  9.1],
       'R': [1.8,  9.0,  12.5],
       'S': [2.2,  9.2,  13.0],
       'T': [2.1,  9.1,  13.0],
       'V': [2.2,  9.7],
       'W': [2.4,  9.4],
       'Y': [2.2,  9.1,  10.1]}

def calculate_charge_from_sequence(sequence, pH=7.0):
    charge = 0.0
    for res in sequence:
        if res == 'D':
            if pH > pKs['D'][2]:
                charge -= 1.0
        elif res == 'E':
            if pH > pKs['E'][2]:
                charge -= 1.0
        elif res == 'K':
            if pH < pKs['K'][2]:
                charge += 1.0
        elif res == 'R':
            if pH < pKs['R'][2]:
                charge += 1.0
        elif res == 'H':
            if pH < pKs['H'][2]:
                charge += 1.0
    return charge

def calculate_pI_from_sequence(sequence):
    charge = calculate_charge_from_sequence(sequence)
    counts = {'D':0,'E':0,
              'K':0,'R':0,'H':0}
    for res in sequence:
        if res in counts.keys():
            counts[res] += 1
    my_pks = []
    my_pks.append(pKs[sequence[0]][1])
    my_pks.append(pKs[sequence[-1]][0])
    for res in sequence:
        if res in counts.keys():
            my_pks.append(pKs[res][2])
    my_pks.sort()
    pk1_index = int(math.floor(len(my_pks)/2.0)+charge)-1
    if pk1_index >= len(my_pks)-1:
        PI = my_pks[-1]
    elif pk1_index <= 0:
        PI = my_pks[0]
    else:
        PI = (my_pks[pk1_index] + my_pks[pk1_index+1])/2.0
    return PI
        
                