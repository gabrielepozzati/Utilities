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

class SelectInterface(Select):
     def accept_residue(self, residue):
         if residue in accept:
             return True
         else:
             return False


if len(sys.argv) < 4:
    print ('Missing arguments! Usage:')
    print ('python interface_extract.py structure1_path structure2_path output_path')
    sys.exit()

io = PDBIO()
p = PDBParser()
complex_name = sys.argv[1].rstrip('.rank')
tmp_folder = '/scratch/gramm_ranks/'+complex_name+'/'

receptor = p.get_structure('', sys.argv[1])
receptor_residues = get_residues(receptor)

ligand = p.get_structure('', sys.argv[2])
ligand_residues = get_residues(ligand)

receptor_interface = interface_extract(receptor_residues, ligand_residues)
ligand_interface = interface_extract(ligand_residues, receptor_residues)

accept = receptor_interface
io.set_structure(receptor)
io.save(sys.argv[3]+sys.argv[1].split('/')[-1].rstrip('.pdb')+'_i1', SelectInterface())

accept = ligand_interface
io.set_structure(ligand)
io.save(sys.argv[3]+sys.argv[2].split('/')[-1].rstrip('.pdb')+'_i2', SelectInterface())

