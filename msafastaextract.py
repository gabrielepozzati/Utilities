import sys

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','UNK':'X'}

seq = ''
prv = ''
for line in open(sys.argv[1],'r'):
    if line.startswith('ATOM') and line[22:27].strip() != prv:
        if line[17:20] in three2one: seq += three2one[line[17:20]]
        else: seq += 'X'
        prv = line[22:27].strip()

with open(sys.argv[2].rstrip('/')+'/'+sys.argv[1].rstrip('.pdb').split('/')[-1]+'.fasta', 'w') as outfile:
    outfile.write('>'+sys.argv[1].rstrip('.pdb')+'\n')
    for char in seq: outfile.write(char)
