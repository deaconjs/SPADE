infile  = open('./conservation_out.txt', 'r')
outfile = open('./conservation_analysis.out', 'w')
lines = infile.readlines()
outlines = []
for line in lines:
    if 'Opening pdb file' in line:
        outlines.append('\n' + line)
    elif line[:5] == 'chain':
        outlines.append(line)
    elif 'angstroms under' in line:
        outlines.append(line)
        
outfile.writelines(outlines)
outfile.close()
