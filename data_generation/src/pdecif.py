import os
import shutil
import math
from rdkit import Chem
from rdkit import RDLogger

# Disable RDKit warnings.
RDLogger.DisableLog('rdApp.*')

def read_protein_ecif_map():
    pemap = {}
    with open('../input/PDB_Atom_Keys.csv', 'r') as f:
        for line in f:
            nline = line.strip().split(',')
            amino_acid_atom_type = nline[0]
            ecif_val = nline[1]
            pemap[amino_acid_atom_type] = ecif_val
    
    return pemap

def read_molecule_ids(subset):
    filepath = '../input/ligands/' + subset + '/'
    filenames = [name.split('.')[0] for name in os.listdir(filepath)]
    filenames.sort()

    return filenames

def read_ligand(subset, molecule_id):
    filepath = '../input/ligands/' + subset + '/' + molecule_id + '.mol'
    return Chem.MolFromMolFile(filepath)

def make_property(atom_details):
    return ';'.join([str(val) for val in atom_details])

def process_ligand(mol):
    lproperties = {}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H':
            details = (atom.GetSymbol(),
                       atom.GetExplicitValence(),
                       len([natom.GetIdx() for natom in atom.GetNeighbors()
                        if natom.GetSymbol() != 'H']),
                       len([natom.GetIdx() for natom in atom.GetNeighbors()
                        if natom.GetSymbol() == 'H']),
                       int(atom.GetIsAromatic()),
                       int(atom.IsInRing()))
            
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            atom_ecif = make_property(details)
            lproperties[atom.GetIdx()] = (atom_ecif, pos.x, pos.y, pos.z)

    return lproperties

def process_protein(subset, molecule_id, pemap):
    protein = {}
    filepath = '../input/proteins/' + subset + '/' + molecule_id + '_pocket.mol2'
    with open(filepath, 'r') as f:
        in_atoms = False
        in_bonds = False
        for line in f:
            line = line.strip()
                        
            if line == '@<TRIPOS>ATOM':
                in_atoms = True
                continue
            if line == '@<TRIPOS>BOND':
                in_atoms = False
                in_bonds = True
                continue

            if in_atoms: 
                atom_line = line.split()
                index = atom_line[0]
                atom_name = atom_line[1]
                x_coord = float(atom_line[2])
                y_coord = float(atom_line[3])
                z_coord = float(atom_line[4])
                amino_acid = atom_line[7][:3]
                aapair = amino_acid + '-' + atom_name

                try:
                    atom_ecif = pemap[aapair]
                    protein[index] = (atom_ecif, x_coord, y_coord, z_coord)
                except KeyError:
                    pass

    return protein

def get_molecule_feature(protein_properties, ligand_properties, distance):        
    protein_ecif, px, py, pz = protein_properties
    ligand_ecif, lx, ly, lz = ligand_properties
    edist = math.sqrt((px-lx)**2 + (py-ly)**2 + (pz-lz)**2)

    dist = ''
    if edist <= distance:
        dist = 'l'
    elif edist > distance:
        dist = 'h'
    else:
        return None

    return protein_ecif + '-' + ligand_ecif + '-' + dist

def get_molecule_features(protein, ligand, distance):
    molecule_features = {}
    for pidx, protein_properties in protein.items():
        for lidx, ligand_properties in ligand.items():
            molfeature = get_molecule_feature(protein_properties, 
                                              ligand_properties, 
                                              distance)
                        
            if molfeature:
                try:
                    molecule_features[molfeature] += 1
                except KeyError:
                    molecule_features[molfeature] = 1

    return molecule_features

def write_dataset(subset, molecule_ids, all_features, all_molecule_features, distance):
    all_features = sorted(all_features)
    filepath = '../output/pecif_' + subset.lower() + '_' + str(distance) + '.csv'

    header = ',' + ','.join(all_features)
    with open(filepath, 'w') as f:
        f.write(header + '\n')
        for molecule_id in molecule_ids:
            molecule_features = all_molecule_features[molecule_id]
            
            vals = []
            for molfeature in all_features:
                try:
                    vals.append(str(molecule_features[molfeature]))
                except KeyError:
                    vals.append(str(0))
            mol_line = molecule_id + ',' + ','.join(vals)
            f.write(mol_line + '\n')

def generate_ecif_dataset(subset, pemap, distance):
    molecule_ids = read_molecule_ids(subset)

    all_features = set()
    all_molecule_features = {}
    for molecule_id in molecule_ids:
        ligand = read_ligand(subset, molecule_id)
        ligand = process_ligand(ligand)
        protein = process_protein(subset, molecule_id, pemap)
        molecule_features = get_molecule_features(protein, ligand, distance)
        all_molecule_features[molecule_id] = molecule_features
        all_features = all_features.union(molecule_features.keys())

    write_dataset(subset, molecule_ids, all_features, all_molecule_features, distance)

def main():
    pemap = read_protein_ecif_map()
    distances = (4,6,8,10)
    subsets = ['CASF-07', 'CASF-13', 'CASF-16', 'CASF-19']

    for subset in subsets:
        print(subset)
        for distance in distances:
            print(distance)
            generate_ecif_dataset(subset, pemap, distance)

if __name__ == '__main__':
    main()