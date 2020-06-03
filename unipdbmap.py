import sys
import Bio
import subprocess
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBIO import Select


def get_pdbchains(uniprotcode):
    pdblist = []
    chaindic = {}
    subprocess.run('wget https://www.uniprot.org/uniprot/'+uniprotcode+'.txt', shell=True)
    try:
        for line in open(uniprotcode+'.txt', 'r'):
 
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

        subprocess.run('rm '+uniprotcode+'.txt', shell=True)
        return pdblist, chaindic

    except Exception as e:
        print(e)
        subprocess.run('rm '+uniprotcode+'.txt', shell=True)
        return pdblist, chaindic

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
p = PDBParser()
if len(sys.argv) < 4:
    print ('Usage: python unipdbmap.py [codelist] [outpath] [mode=single/pair/linkedpair]')
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

outpath = sys.argv[2]
mode = sys.argv[3]

if mode == 'single':
    for code in open(sys.argv[1],'r'):
        pdblist, chaindict = get_pdbchains(code)

        maxlength = 0 
        longestpdb = ''
        longestchain = ''
        for pdb in pdblist:
                maxlength = chaindict[pdb][1]
                longestchain = chaindict[pdb][0]
                longestpdb = pdb

        subprocess.run('wget https://files.rcsb.org/download/'+longestpdb+'.pdb', shell=True)

        try:
            structure = p.get_structure('', longestpdb+'.pdb')
            m, accept = model_chain_select(structure, longestchain[0])

        except Exception as e: 
            print(e)
            subprocess.run('rm '+longestpdb+'.pdb', shell=True)
            continue

        io.set_structure(structure)
        io.save(outpath+code.rstrip()+'_'+longestpdb+'_'+accept+'.pdb', 
                ChainSelect())

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
        subprocess.run('wget https://files.rcsb.org/download/'+longestpdb1+'.pdb',
                       shell=True)
        subprocess.run('wget https://files.rcsb.org/download/'+longestpdb2+'.pdb',
                       shell=True)

        try:
            structure1 = p.get_structure('', longestpdb1+'.pdb')
            m1, accept1 = model_chain_select(structure1, longestchain1[0])
            
            structure2 = p.get_structure('', longestpdb2+'.pdb')
            m2, accept2 = model_chain_select(structure2, longestchain2[0])

        except Exception as e:
            print(e)
            subprocess.run('rm '+longestpdb1+'.pdb', shell=True)
            subprocess.run('rm '+longestpdb2+'.pdb', shell=True)
            continue

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
    
        subprocess.run('wget https://files.rcsb.org/download/'+longestpdb+'.pdb', 
                       shell=True)

        try:
            structure = p.get_structure('', longestpdb+'.pdb')
            for model in structure:
                m = model.get_id()
                break
            
            maxlength = 0
            maxinterface = ''
            for c1 in chains1.split('/'):
                for c2 in chains2.split('/'):
                    int_residues1 = extinterface(structure[m][c1], structure[m][c2])
                    if len(int_residues1) < 5: continue
                    int_residues2 = extinterface(structure[m][c2], structure[m][c1])
                    if len(int_residues2) < 5: continue
                    int_residues = int_residues1 + int_residues2
                    if len(int_residues) > maxlength: 
                        maxlength = len(int_residues)
                        maxinterface = c1+c2

        except Exception as e: 
            print(e)
            subprocess.run('rm '+longestpdb+'.pdb', shell=True)
            continue
    
        if len(maxinterface)<2: 
            subprocess.run('rm '+longestpdb+'.pdb', shell=True)
            continue
    
        lengthchain1 = len(structure[m][maxinterface[0]])
        lengthchain2 = len(structure[m][maxinterface[1]])
        if lengthchain1 < lengthchain2: maxinterface = maxinterface[1]+maxinterface[0]

        io.set_structure(structure)

        accept = maxinterface[0]
        io.save(outpath+codecouple.rstrip()+'_'+longestpdb+'_'+maxinterface+'_u1.pdb', 
                ChainSelect())

        accept = maxinterface[1]
        io.save(outpath+codecouple.rstrip()+'_'+longestpdb+'_'+maxinterface+'_u2.pdb', 
                ChainSelect())
    
        subprocess.run('rm '+longestpdb+'.pdb', shell=True)

