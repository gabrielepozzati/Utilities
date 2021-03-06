import os
import os.path
import sys
import Bio
import time
import numpy as np
import subprocess
from Bio.PDB import StructureBuilder as SB
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select


def get_pdbchains(uniprotcode):
    pdblist = []
    chaindic = {}
    if uniprotcode in uniprot_dic:
        for line in uniprot_dic[uniprotcode]:
 
            if line.startswith('DR   PDB; '):
                fields = line.split('; ')
                fields[-1] = fields[-1].rstrip().rstrip('.')

                pdbcode = fields[1]

                maxlength = 0
                bestchaingroup = ''
                for chaingroup in fields[-1].split(','):
                    start = int(chaingroup.split('=')[1].split('-')[0])
                    end = int(chaingroup.split('=')[1].split('-')[1])
                    if (end-start)+1 > maxlength: 
                        bestchaingroup = chaingroup.split('=')[0]
                        maxlength = (end-start)+1

                pdblist.append(pdbcode)
                chaindic[pdbcode] = [bestchaingroup, maxlength]

        return pdblist, chaindic

    else: 
        print ("Unavailable data on UniprotKB!")
        return pdblist, chaindic

def derive_assembly(pdbpath, allowed_chains):

    rotation = []
    translation = []
    matrixgroup = []
    coordinates = {}
    rototranslations = {}
    chainlist = ''
    for line in open(pdbpath, 'r'):
        chainsrow = False
        matrixrow = False
        if line.startswith('REMARK 350'): r350 = True
        else: r350 = False

        if r350:
            if line.startswith('REMARK 350 APPLY THE FOLLOWING'): 
                chainsrow = True
                chainlist = []

            if chainsrow:
                line = line.split(': ')[-1].strip().rstrip(',')
                chainlist += line.split(',')
                for pos in range(len(chainlist)): chainlist[pos] = chainlist[pos].strip()

            if line.startswith('REMARK 350   BIOMT'):
                chainsrow = False
                matrixrow = True

            if matrixrow:
                model = line.split()[3]
                rototranslations[model] = rototranslations.get(model, [])
                rotation.append(line.split()[-4:-1])
                for pos in range(3): rotation[-1][pos] = float(rotation[-1][pos])
                translation.append(float(line.split()[-1]))

                if len(rotation) == 3:
                    rt_matrix = [list(rotation), list(translation)]
                    added = False
                    for pos in range(len(rototranslations[model])):
                        if rt_matrix in rototranslations[model][pos]:
                            rototranslations[model][pos][0] += chainlist
                            added = True
                    if not added: rototranslations[model] += [[chainlist, rt_matrix]]
                    translation = []
                    rotation = []

        
        if line.startswith('ATOM') or line.startswith('TER'): 
            coordinates[line[21]] = coordinates.get(line[21], [])
            coordinates[line[21]].append(line)

    rt_structures = []

    #print (pdbpath)
    #print (rototranslations)
    tmppath = pdbpath.split('/')[-1]+'tmpfile.pdb'
    original_structure = p.get_structure('', pdbpath)
    identity_rtmatrix = [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0]]

    for assembly in rototranslations:
        for rtmatrix in rototranslations[assembly]:
            if rtmatrix[1] == identity_rtmatrix:
                for chain in rtmatrix[0]:
                    if chain not in allowed_chains: continue
                    rt_structures.append(original_structure[0][chain])
                    continue
            else:
                for chain in rtmatrix[0]:
                    structure = get_rtstructure(tmppath, rtmatrix[1], coordinates[chain])
                    rt_structures.append(structure[0][chain])
                    subprocess.run('rm '+tmppath, shell=True)

    return rt_structures

def get_rtstructure(path, matrix, chain):

    outfile = open(path,'w')
    rotation_matrix = np.array(matrix[0])
    translation_matrix = np.array(matrix[1])

    for line in chain:
        if line.startswith('TER'): 
            outfile.write(line)
            continue
        atom_coordinates = np.array([float(line[30:38].strip()), 
                                     float(line[38:46].strip()), 
                                     float(line[46:54].strip())])
        rtcoords = np.add(rotation_matrix.dot(atom_coordinates), translation_matrix)
        rtX = str(np.around(rtcoords[0], decimals=3))
        rtY = str(np.around(rtcoords[1], decimals=3))
        rtZ = str(np.around(rtcoords[2], decimals=3))
        while len(rtX) < 8: rtX = ' '+rtX
        while len(rtY) < 8: rtY = ' '+rtY
        while len(rtZ) < 8: rtZ = ' '+rtZ
        outfile.write(line[:30]+rtX+rtY+rtZ+line[54:])
    outfile.close()

    return p.get_structure('', path)        

def extinterface(ichain1, ichain2):

    iatoms = []
    int_res = []
    for residue in ichain2:
        for atom in residue: iatoms.append(atom)

    ns = Bio.PDB.NeighborSearch(iatoms)
    for residue in ichain1:
        for atom in residue:
            neighbours = ns.search(atom.coord, 6, level='R')
            if len(list(neighbours)) != 0:
                int_res.append(residue)
                break

    return int_res

def model_chain_select(structure, targetchain):
    chainhit = None

    for model in structure:
        m = model.get_id()
        break

    for chain in structure[m]:
        if chain.get_id() == targetchain:
            chainhit = chain.get_id()
            break

    return m, chainhit

class ChainSelect(Select):

     def accept_chain(self, chain):
         if chain.get_id() == accept:
             return True
         else:
             return False

     def accept_model(self, model):
         if model.get_id() == m:
             return True
         else:
             return False

io=PDBIO()
p = PDBParser(QUIET=True)
if len(sys.argv) < 5:
    print ('Usage: python unipdbmap.py [codelist] [uniprot_txt_file] [outpath] [mode=single/pair/linkedpair]')
    print ('mode single:')
    print ('---downloads the uniprot id related pdb with the longest available chain.')
    print ('---codelist must contain one uniprot id per line')
    print ('mode pair:')
    print ('---downloads the uniprot-related pdbs with the longest chain if any pdb is available for both uniprot ids.')
    print ('---codelist must contain two uniprot ids for line separated by _')
    print ('mode linkedpair:')
    print ('---downloads only pdbs containing one interface between the uniprot ids in a pair.')
    print ('---codelist must contain two uniprot ids for line separated by _')
    sys.exit()

upcode = ''
uniprot_dic = {'':[]}
for line in open(sys.argv[2],'r'):
    if line.startswith("AC   "): 
        upcode = line[5:-1].split(';')
        for code in upcode: uniprot_dic[code.strip()] = []
    else: 
        for code in upcode: uniprot_dic[code.strip()].append(line.rstrip())

outpath = sys.argv[3]
mode = sys.argv[4]
try: pdbdb_path = os.environ['PDBDB']
except: 
    print ('WARNING! Local PDB database path not declared!') 
    print ('To use a local version of the PDB database, declare the path containing all the pdb subfolders (00/ 27/ 51/ .. ) by running:') 
    print ('export PDBDB=\'local_pdbfolder_path\'')
    print ('Starting anyway in 5 seconds ...')
    time.sleep(5)

if mode == 'single':
    outfile = open('largestructures','w')
    for code in open(sys.argv[1],'r'):
        pdblist, chaindict = get_pdbchains(code.rstrip())
        
        if len(pdblist) == 0: 
            print ('No PDB found!')
            continue

        maxlength = 0 
        longestpdb = ''
        longestchain = ''
        for pdb in pdblist:
            if chaindict[pdb][1] > maxlength:
                maxlength = chaindict[pdb][1]
                longestchain = chaindict[pdb][0]
                longestpdb = pdb
        
        if os.path.exists(pdbdb_path+'/'+longestpdb[1:3].lower()+'/'+'pdb'+longestpdb.lower()+'.ent'): 
            structure = p.get_structure('', pdbdb_path+'/'+longestpdb[1:3].lower()+'/'+'pdb'+longestpdb.lower()+'.ent')
        else:
            subprocess.run('wget https://files.rcsb.org/download/'+longestpdb+'.pdb', shell=True)
            if os.path.exists(longestpdb+'.pdb'): structure = p.get_structure('', longestpdb+'.pdb')
            else:
                outfile.write(code+'_'+longestpdb+'_'+longestchain.split('/')[0]+'\n')
                continue
        continue

        m, accept = model_chain_select(structure, longestchain[0]) 
        if accept == None: continue

        print ('Writing! PDB: '+longestpdb+' Length: '+str(maxlength)+'('+longestchain+')')

        io.set_structure(structure)
        io.save(outpath+code.rstrip()+'_'+longestpdb+'_'+accept+'.pdb', 
                ChainSelect())
    outfile.close()


if mode == 'pair':
    
    for codecouple in open(sys.argv[1],'r'):
        code1 = codecouple.split('_')[0]
        code2 = codecouple.split('_')[1].rstrip()
        pdblist1, chaindict1 = get_pdbchains(code1)
        pdblist2, chaindict2 = get_pdbchains(code2)

        maxlength1 = 0
        longestpdb1 = ''
        longestchain1 = ''
        for pdb in pdblist1:
            if chaindict1[pdb][1] > maxlength1:
                maxlength1 = chaindict1[pdb][1]
                longestchain1 = chaindict1[pdb][0]
                longestpdb1 = pdb

        maxlength2 = 0
        longestpdb2 = ''
        longestchain2 = ''
        for pdb in pdblist2:
            if chaindict2[pdb][1] > maxlength2:
                maxlength2 = chaindict2[pdb][1]
                longestchain2 = chaindict2[pdb][0]
                longestpdb2 = pdb

        if longestpdb1 == '' or longestpdb2 == '': continue

        try: structure1 = p.get_structure('', pdbdb_path+'/'+longestpdb1[1:3].lower()+'/'+'pdb'+longestpdb.lower()+'.ent')
        except: 
            subprocess.run('wget https://files.rcsb.org/download/'+longestpdb1+'.pdb', shell=True)
            if os.path.exists(longestpdb1+'.pdb'): structure1 = p.get_structure('', longestpdb1+'.pdb')
            else: continue

        try: structure2 = p.get_structure('', pdbdb_path+'/'+longestpdb2[1:3].lower()+'/'+'pdb'+longestpdb.lower()+'.ent')
        except: 
            subprocess.run('wget https://files.rcsb.org/download/'+longestpdb2+'.pdb', shell=True)
            if os.path.exists(longestpdb2+'.pdb'): structure2 = p.get_structure('', longestpdb2+'.pdb')
            else: continue

        
        m1, accept1 = model_chain_select(structure1, longestchain1[0])
        m2, accept2 = model_chain_select(structure2, longestchain2[0])

        lengthchain1 = len(structure1[m1][accept1])
        lengthchain2 = len(structure2[m2][accept2])
        if lengthchain1 > lengthchain2:

            io.set_structure(structure1)
            m = m1
            accept = accept1
            io.save(outpath+codecouple.rstrip()+'_'+longestpdb1+accept1+'_'+longestpdb2+accept2+'_u1.pdb',
                    ChainSelect())

            io.set_structure(structure2)
            m = m2
            accept = accept2
            io.save(outpath+codecouple.rstrip()+'_'+longestpdb1+accept1+'_'+longestpdb2+accept2+'_u2.pdb',
                    ChainSelect())
        else:

            io.set_structure(structure2)
            m = m2
            accept = accept2
            io.save(outpath+codecouple.rstrip()+'_'+longestpdb2+accept2+'_'+longestpdb1+accept1+'_u1.pdb',
                    ChainSelect())

            io.set_structure(structure1)
            m = m1
            accept = accept1
            io.save(outpath+codecouple.rstrip()+'_'+longestpdb2+accept2+'_'+longestpdb1+accept1+'_u2.pdb',
                    ChainSelect())

        subprocess.run('rm '+longestpdb1+'.pdb', shell=True)
        subprocess.run('rm '+longestpdb2+'.pdb', shell=True)

if mode == 'linkedpair':
    
    for codecouple in open(sys.argv[1],'r'):
        print ('\n\n\n'+codecouple.rstrip())
        code1 = codecouple.split('_')[0]
        code2 = codecouple.split('_')[1].rstrip()
        pdblist1, chaindict1 = get_pdbchains(code1)
        pdblist2, chaindict2 = get_pdbchains(code2)

        pdblist1=set(pdblist1)
        pdblist2=set(pdblist2)
        commonpdb = list(pdblist1.intersection(pdblist2))

        maxlength = 0
        longestpdb = ''
        longestchain = ''
        for pdb in commonpdb:
            if chaindict1[pdb][1]+chaindict2[pdb][1] > maxlength: 
                maxlength = chaindict1[pdb][1]+chaindict2[pdb][1]
                longestchain = chaindict1[pdb][0]+'-'+chaindict2[pdb][0]
                longestpdb = pdb

        if len(longestchain.split('-')) < 2: continue
        chains1 = longestchain.split('-')[0]
        chains2 = longestchain.split('-')[1]

        print ('PDB: '+longestpdb+' Length: '+str(maxlength)+'('+longestchain+')')
        if os.path.exists(pdbdb_path+'/'+longestpdb[1:3].lower()+'/'+'pdb'+longestpdb.lower()+'.ent'):
            reference_path = pdbdb_path+'/'+longestpdb[1:3].lower()+'/'+'pdb'+longestpdb.lower()+'.ent'
        else:
            subprocess.run('wget https://files.rcsb.org/download/'+longestpdb+'.pdb', shell=True)
            if os.path.exists(longestpdb+'.pdb'): reference_path = longestpdb+'.pdb'
            else: continue

        #print ('chain1:', chains1.split('/'), 'chain2:',chains2.split('/') )
        chain_list1 = derive_assembly(reference_path, chains1.split('/'))
        chain_list2 = derive_assembly(reference_path, chains2.split('/'))
        #for chain in chain_list: 
        #    print (chain.get_id())
        #continue

        maxlength = 0
        maxinterface1 = ''
        maxinterface2 = ''
        for c1 in range(len(chain_list1)):
            for c2 in range(len(chain_list2)):

                firstatom1 = list(list(chain_list1[c1])[0])[0].get_coord()
                firstatom2 = list(list(chain_list2[c2])[0])[0].get_coord()
                if np.all(firstatom1==firstatom2): continue

                int_residues1 = extinterface(chain_list1[c1], chain_list2[c2])
                if len(int_residues1) < 5: continue
                int_residues2 = extinterface(chain_list2[c2], chain_list1[c1])
                if len(int_residues2) < 5: continue
                int_residues = int_residues1 + int_residues2
                if len(int_residues) > maxlength: 
                    maxlength = len(int_residues)
                    maxinterface1 = chain_list1[c1]
                    maxinterface2 = chain_list2[c2]

        if maxinterface1 == '' or maxinterface2 == '': continue
        chaincouple = maxinterface1.get_id()+maxinterface2.get_id()

        print ('Wrinting '+outpath+codecouple.rstrip()+'_'+longestpdb+'_'+chaincouple+'....')
        lengthchain1 = len(maxinterface1)
        lengthchain2 = len(maxinterface2)
        if lengthchain1 < lengthchain2:
            swapinterface = maxinterface2
            maxinterface2 = maxinterface1
            maxinterface1 = swapinterface

        io.set_structure(maxinterface1)
        io.save(outpath+codecouple.rstrip()+'_'+longestpdb+'_'+chaincouple+'_u1.pdb')

        io.set_structure(maxinterface2)
        io.save(outpath+codecouple.rstrip()+'_'+longestpdb+'_'+chaincouple+'_u2.pdb')
    
        if os.path.exists(longestpdb+'.pdb'): subprocess.run('rm '+longestpdb+'.pdb', shell=True)

