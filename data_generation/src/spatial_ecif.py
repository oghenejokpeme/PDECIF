import os
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

def read_molecule_filenames(subset):
    filepath = '../input/' + subset + '/'
    filenames = os.listdir(filepath)

    seen = set()
    molecule_ids = []
    for filename in filenames:
        molecule_id = filename.split('_')[0]
        if molecule_id not in seen:
            molecule_ids.append(molecule_id)
            seen.add(molecule_id)
    molecule_ids.sort()

    return molecule_ids

def read_ligand_molfile(subset, molecule_id):
    filepath = '../input/' + subset + '/' + molecule_id + '_lig.mol'

    if os.path.isfile(filepath):
        return Chem.MolFromMolFile(filepath)
    else:
        return None

def read_ligand_mol2file(subset, molecule_id):
    filepath = '../input/' + subset + '/' + molecule_id + '_lig.mol2'
    
    if os.path.isfile(filepath):
        return Chem.MolFromMol2File(filepath)
    else:
        return None

def read_protein_mol2file(subset, molecule_id, pemap):
    protein = {}
    filepath = '../input/' + subset + '/' + molecule_id + '_pocket.mol2'
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

def make_property(atom_details):
    return ';'.join([str(val) for val in atom_details])

def generate_ligand_properties(mol, protonate):
    if protonate: mol = Chem.AddHs(mol)
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

def get_molecule_feature(protein_properties, ligand_properties, rdist = (6, 12)):        
    protein_ecif, px, py, pz = protein_properties
    ligand_ecif, lx, ly, lz = ligand_properties
    distance = math.sqrt((px-lx)**2 + (py-ly)**2 + (pz-lz)**2)

    dist = ''
    rmin, rmax = rdist
    if distance <= rmin:
        dist = 'low'
    elif distance > rmin and distance <= rmax:
        dist = 'mid'
    else:
        return None

    return protein_ecif + '-' + ligand_ecif + '-' + dist

def get_molecule_features(protein, ligand):
    molecule_features = {}
    for pidx, protein_properties in protein.items():
        for lidx, ligand_properties in ligand.items():
            molfeature = get_molecule_feature(protein_properties, ligand_properties)
                        
            if molfeature:
                try:
                    molecule_features[molfeature] += 1
                except KeyError:
                    molecule_features[molfeature] = 1

    return molecule_features

def write_dataset(subset, molecule_ids, all_features, all_molecule_features, protonate):
    all_features = sorted(all_features)
    filepath = ''
    if protonate:
        filepath = '../output/x_spatial_ecif_with_hs_' + subset.lower() + '.csv'
    else:
        filepath = '../output/x_spatial_ecif_without_hs_' + subset.lower() + '.csv'

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

def generate_dataset(subset, pemap, molecule_ids, protonate):
    total_molecules = len(molecule_ids)
    print(subset)
    all_features = set()
    all_molecule_features = {}
    for midx, molecule_id in enumerate(molecule_ids, 1):
        message = f'{midx}/{total_molecules} {molecule_id}'
        print(message)
        ligand = None
        # If mol2 ligand file cannot be read by RDKit, mol file
        mol_lig = read_ligand_mol2file(subset, molecule_id)
        if not mol_lig:
            mol_lig = read_ligand_molfile(subset, molecule_id)
            if not mol_lig:
                print('Incorrect mol and mol2 ligand file')
                continue
            ligand = generate_ligand_properties(mol_lig, protonate)
        else:
            ligand = generate_ligand_properties(mol_lig, protonate)
        protein = read_protein_mol2file(subset, molecule_id, pemap)
        molecule_features = get_molecule_features(protein, ligand)
        all_molecule_features[molecule_id] = molecule_features
        all_features = all_features.union(molecule_features.keys())

    write_dataset(subset, molecule_ids, all_features, all_molecule_features, protonate)

def main():
    pemap = read_protein_ecif_map()
    subsets = ['CASF-07', 'CASF-13', 'CASF-16']
    for subset in subsets:
        message = f'Processing subset: {subset}'
        print(message)    
        molecule_ids = read_molecule_filenames(subset)
        generate_dataset(subset, pemap, molecule_ids, True)
        
if __name__ == '__main__':
    main()