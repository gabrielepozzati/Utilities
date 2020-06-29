import sys

outfile = open(sys.argv[2],'w')
with open(sys.argv[1],'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            outline = line[:21]+sys.argv[3]+line[22:].rstrip()
            outfile.write(outline+'\n')
            outline = ''
outfile.close() 
        
