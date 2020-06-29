import sys


new = 1
prevorig = '' 
outfile = open(sys.argv[2],'w')
with open(sys.argv[1],'r') as f:
    for line in f:                 

        if line.startswith('ATOM'): 
            if prevorig != '' and line[22:27] != prevorig: new+=1
            renumb = ' '*(4-len(str(new)))+str(new)+' '
            outline = line[:22]+renumb+line[27:].rstrip()
            outfile.write(outline+'\n')
        else: outfile.write(line.rstrip()+'\n')

        prevorig = line[22:27]
        
outfile.close()
        
