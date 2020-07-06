import sys
import subprocess
from Bio.PDB import *


def interface_extract(target, coupled):

    coupled_atoms = []
    target_interface = []
    for residue in coupled:
        for atom in residue: 
            coupled_atoms.append(atom)

    ns = NeighborSearch(coupled_atoms)

    for residue in target:
        for atom in residue:
            neighbours = ns.search(atom.coord, 12, level='R')
            if len(list(neighbours)) != 0:
                target_interface.append(residue)
                break

    return target_interface

def get_residues(structure):

    for model in structure:
        m = model.get_id()
        for chain in model:
            c = chain.get_id()
            break
        break

    return Selection.unfold_entities(structure[m][c], 'R')

class Select(Select):
     def accept_residue(self, residue):
         if residue in accept:
             return True
         else:
             return False


countB = 0
badpose = []
goodpose = []

for line in open(sys.argv[1], 'r'):

    if line.startswith('#'): continue
    pose = line.split()[0]
    quality = line.split()[1]

    if quality == 'I': 
        if countB>=10: continue
        badpose.append(pose)
        countB +=1
    else: goodpose.append(pose)

io = PDBIO()
p = PDBParser()
complex_name = sys.argv[1].rstrip('.rank')
tmp_folder = '/scratch/gramm_ranks/'+complex_name+'/'

ligand_name = complex_name.split('-')[1]
receptor_name = complex_name.split('-')[0]
receptor = p.get_structure('', receptor_name+'.pdb')
receptor_residues = get_residues(receptor)

for pose in goodpose:
    pose_name = tmp_folder+complex_name+'_'+pose+'.pdb'
    ligand = p.get_structure('', pose_name)
    ligand_residues = get_residues(ligand)

    receptor_interface = interface_extract(receptor_residues, ligand_residues)
    ligand_interface = interface_extract(ligand_residues, receptor_residues)

    accept = receptor_interface
    io.set_structure(receptor)
    io.save(tmp_folder+'good/'+receptor_name+'_i'+pose, Select())
    
    accept = ligand_interface
    io.set_structure(ligand)
    io.save(tmp_folder+'good/'+ligand_name+'_i'+pose, Select())

for pose in badpose:
    pose_name = tmp_folder+complex_name+'_'+pose+'.pdb'
    ligand = p.get_structure('', pose_name)
    ligand_residues = get_residues(ligand)

    receptor_interface = interface_extract(receptor_residues, ligand_residues)
    ligand_interface = interface_extract(ligand_residues, receptor_residues)

    accept = receptor_interface
    io.set_structure(receptor)
    io.save(tmp_folder+'bad/'+receptor_name+'_i'+pose, Select())
    
    accept = ligand_interface
    io.set_structure(ligand)
    io.save(tmp_folder+'bad/'+ligand_name+'_i'+pose, Select())

    


