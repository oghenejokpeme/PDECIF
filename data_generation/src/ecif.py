import os
import shutil
import math
from rdkit import Chem
from rdkit import RDLogger

# Dictionaries extracted from https://github.com/DIFACQUIM/ECIF
# Paper: https://doi.org/10.1093/bioinformatics/btaa982
ECIF_ProteinAtoms = {'C;4;1;3;0;0', 'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1',
                     'C;4;3;0;0;0', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                     'C;5;3;0;0;0', 'C;6;3;0;0;0', 'N;3;1;2;0;0', 'N;3;2;0;1;1',
                     'N;3;2;1;0;0', 'N;3;2;1;1;1', 'N;3;3;0;0;1', 'N;4;1;2;0;0',
                     'N;4;1;3;0;0', 'N;4;2;1;0;0', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                     'S;2;1;1;0;0', 'S;2;2;0;0;0'}

# Possible ligand atoms according to the PDBbind 2016 "refined set"
ECIF_LigandAtoms = {'Br;1;1;0;0;0', 'C;3;3;0;1;1', 'C;4;1;1;0;0', 'C;4;1;2;0;0',
                     'C;4;1;3;0;0', 'C;4;2;0;0;0', 'C;4;2;1;0;0', 'C;4;2;1;0;1',
                     'C;4;2;1;1;1', 'C;4;2;2;0;0', 'C;4;2;2;0;1', 'C;4;3;0;0;0',
                     'C;4;3;0;0;1', 'C;4;3;0;1;1', 'C;4;3;1;0;0', 'C;4;3;1;0;1',
                     'C;4;4;0;0;0', 'C;4;4;0;0;1', 'C;5;3;0;0;0', 'C;5;3;0;1;1',
                     'C;6;3;0;0;0', 'Cl;1;1;0;0;0', 'F;1;1;0;0;0', 'I;1;1;0;0;0',
                     'N;3;1;0;0;0', 'N;3;1;1;0;0', 'N;3;1;2;0;0', 'N;3;2;0;0;0',
                     'N;3;2;0;0;1', 'N;3;2;0;1;1', 'N;3;2;1;0;0', 'N;3;2;1;0;1',
                     'N;3;2;1;1;1', 'N;3;3;0;0;0', 'N;3;3;0;0;1', 'N;3;3;0;1;1',
                     'N;4;1;2;0;0', 'N;4;1;3;0;0', 'N;4;2;1;0;0', 'N;4;2;2;0;0',
                     'N;4;2;2;0;1', 'N;4;3;0;0;0', 'N;4;3;0;0;1', 'N;4;3;1;0;0',
                     'N;4;3;1;0;1', 'N;4;4;0;0;0', 'N;4;4;0;0;1', 'N;5;2;0;0;0',
                     'N;5;3;0;0;0', 'N;5;3;0;1;1', 'O;2;1;0;0;0', 'O;2;1;1;0;0',
                     'O;2;2;0;0;0', 'O;2;2;0;0;1', 'O;2;2;0;1;1', 'P;5;4;0;0;0',
                     'P;6;4;0;0;0', 'P;6;4;0;0;1', 'P;7;4;0;0;0', 'S;2;1;0;0;0',
                     'S;2;1;1;0;0', 'S;2;2;0;0;0', 'S;2;2;0;0;1', 'S;2;2;0;1;1',
                     'S;3;3;0;0;0', 'S;3;3;0;0;1', 'S;4;3;0;0;0', 'S;6;4;0;0;0',
                     'S;6;4;0;0;1', 'S;7;4;0;0;0'}

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
    mol = Chem.AddHs(mol)
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
            if atom_ecif in ECIF_LigandAtoms:
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
    
    if edist <= distance:
        return f'{protein_ecif}-{ligand_ecif}'
    else:
        return None

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
    filepath = '../output/ecif_' + subset.lower() + '_' + str(distance) + '.csv'
    print(f'  {len(all_features)}')
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
        if not ligand:
            print(f'{ligand} Empty!')
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
            print()

if __name__ == '__main__':
    main()