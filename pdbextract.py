import sys

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','UNK':'X'}

seq = ''
prv = ''
structure = []

chain = ''
model = ''
if sys.argv[2] != 'p' and sys.argv[2] != 'f':
    print ("Second command line argument have to be f (for fasta format output)")
    print ("or p (for pdb format output).")
    sys.exit()
if len(sys.argv) > 3: chain = str(sys.argv[3])
else: print ("WARNING: only first chain of first model will be picked in this run")
if len(sys.argv) > 4: model = str(sys.argv[4])
else: print ("WARNING: only selected chain of the first model will be picked in this run")
if model == '0': model = ''
rightmodel = 'n'

for line in open(sys.argv[1],'r'):

    if model != '' and line.startswith('MODEL'):
        if line.split()[1] == model: rightmodel = 'y'
    if model != '' and line.startswith('ENDMDL'): rightmodel = 'n'
    if model != '' and rightmodel == 'n': continue

    
    if line.startswith('ATOM'):
        if chain != '': 
            if line[21] != chain: continue 
        if line[22:27].strip() != prv:
            seq += three2one[line[17:20]]
            prv = line[22:27].strip()
        structure.append(line)

    if chain == '' and model == '': 
        if line.startswith('TER'): break
    if model == '':
        if line.startswith('ENDMDL'): break

if len(structure) == 0 or len(seq) == 0:
    print ("Chain "+chain+", Model "+model+" not found!")
    print ("Command line arguments have to be in the order: pdb_file, out_type(p/f) [chain,] [model]")
    sys.exit()

outname = sys.argv[1].split('/')[-1].split('.')[0]
if chain != '': outname += '_'+chain
if model != '': outname += '_'+model

if sys.argv[2] == 'f':
    with open(outname+'.fasta', 'w') as outfile:
        outfile.write('>'+outname+'\n')
        for char in seq: outfile.write(char)
        
elif sys.argv[2] == 'p': 
    with open(outname+'.pdb', 'w') as outfile:
        for line in structure: outfile.write(line)
