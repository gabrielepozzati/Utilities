import sys

outpath = sys.argv[2].rstrip('/')+'/'+sys.argv[1].split('/')[-1].rstrip('.pdb')

outfile1 = open(outpath+'_1.pdb', 'w')
outfile2 = open(outpath+'_2.pdb', 'w')

with open(sys.argv[1], 'r') as f: 
    swap='n'
    for line in f: 

        if line.startswith('TER'):
            swap ='y'
            continue

        if swap == 'n': outfile1.write(line)
        else: outfile2.write(line)

outfile1.close()
outfile2.close()
