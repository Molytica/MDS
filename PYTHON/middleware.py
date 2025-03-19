from rdkit import Chem
from rdkit.Chem import AllChem
import gzip
import os
import plotly.graph_objects as go
from tqdm import tqdm

def smiles_to_structures(smiles):
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    # Generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)

    # Extract atom point cloud and symbols
    atom_coords = [[mol.GetConformer().GetAtomPosition(i).x,
                    mol.GetConformer().GetAtomPosition(i).y,
                    mol.GetConformer().GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())]
    atom_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # Create connection graph
    connection_graph = []
    for bond in mol.GetBonds():
        connection_graph.append({
            "start": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            'weight': bond.GetBondTypeAsDouble()
        })

    return atom_coords, atom_symbols, connection_graph

def get_molecule_data(smiles):
    atom_coords, atom_symbols, connection_graph = smiles_to_structures(smiles)

    return {
        "atom_coords": atom_coords,
        "atom_symbols": atom_symbols,
        "covalent_edge_list": connection_graph
    }


amino_acid_edge_list = {
    "GLY": [(0, 1, 1), (1, 2, 1), (2, 3, 2)],
    "LYS": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1), (6, 7, 1), (7, 8, 1)],
    "ASP": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1.5), (5, 7, 1.5)],
    "VAL": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (3, 6, 1)],
    "ASN": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1), (5, 7, 2)],
    "GLU": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1), (6, 7, 1.5), (6, 8, 1.5)],
    "PHE": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1.5), (5, 7, 1.5), (6, 8, 1.5), (7, 9, 1.5), (9, 10, 1.5), (8, 10, 1.5)],
    "THR": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (3, 6, 1)],
    "GLN": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1), (6, 7, 1), (6, 8, 2)],
    "ILE": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 7, 1), (3, 6, 1)],
    "HIS": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 7, 1.5), (5, 6, 1.5), (6, 9, 1.5), (9, 8, 1.5), (7, 8, 1.5)],
    "LEU": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 7, 1), (5, 6, 1)],
    "SER": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1)],
    "CYS": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1)],
    "TRP": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 7, 1.5), (5, 6, 1.5), (6, 10, 1.5), (10, 8, 1.5), (8, 12, 1.5), (7, 8, 1.5), (7, 9, 1.5), (9, 13, 1.5), (13, 11, 1.5), (12, 11, 1.5)],
    "PRO": [(0, 1, 1), (0, 6, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1)],
    "ARG": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1), (6, 7, 1), (7, 10, 1), (10, 8, 1), (10, 9, 2)],
    "ALA": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1)],
    "MET": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 6, 1), (6, 7, 1)],
    "TYR": [(0, 1, 1), (1, 2, 1), (2, 4, 2), (1, 3, 1), (3, 5, 1), (5, 7, 1.5), (5, 6, 1.5), (6, 8, 1.5), (8, 11, 1.5), (7, 9, 1.5), (9, 11, 1.5), (11, 10, 1)]
}

def get_max_min_edge_list(edge_list):
    max_, min_ = -100, 100
    for edge in edge_list:
        max_ = max(max_, edge[0], edge[1])
        min_ = min(min_, edge[0], edge[1])
    return max_, min_

def get_protein_data(uniprot_id, path=None):
    if path == None:
        pdbgz_file = os.path.join("AFDATA", f"AF-{uniprot_id}-F1-model_v4.pdb.gz")
    else:
        pdbgz_file = path + f"AF-{uniprot_id}-F1-model_v4.pdb.gz"

    atom_coords = []
    atom_types = []

    amino_acids = {}
    should_print = False
    print_count = 0

    last_amino_acid_number = 1
    last_amino_acid_type = None
    last_amino_acid_structure = []

    protein_edge_list = []

    amino_acid_atom_index_dict = {}
    amino_acid_type_dict = {}

    amino_acid_type_atom_list = []

    with gzip.open(pdbgz_file, 'rt') as f:
        atom_idx = -1
        for line in f:
            if line.startswith("ATOM"):
                atom_idx += 1
                # Use regular expression to extract atom type and coordinates
                cols = line.replace("-", " -").split()
                if len(cols[4]) > 2:
                    cols = cols[:4] + [cols[4][:1], cols[4][1:]] + cols[5:]
                    
                    if len(cols) != 12:
                        print(cols)
                                
                atom_type = cols[-1]
                x_coord = float(cols[-6])
                y_coord = float(cols[-5])
                z_coord = float(cols[-4])

                amino_acid_number = int(cols[5])
                amino_acid_label = cols[3]

                # This is done for each atom
                amino_acid_type_atom_list.append(amino_acid_label)

                if amino_acid_number not in amino_acid_atom_index_dict:
                    amino_acid_atom_index_dict[amino_acid_number] = [atom_idx]
                else:
                    amino_acid_atom_index_dict[amino_acid_number].append(atom_idx)

                if amino_acid_number not in amino_acid_type_dict:
                    amino_acid_type_dict[amino_acid_number] = amino_acid_label
                
                last_amino_acid_structure.append([atom_type, x_coord, y_coord, z_coord])
                last_amino_acid_number = amino_acid_number
                last_amino_acid_type = amino_acid_label

                atom_coords.append([x_coord, y_coord, z_coord])
                atom_types.append(atom_type)
                try:
                    pass
                except Exception as e:
                    print(e)
                    print(line)

    # Remove last oxygen
    #print(amino_acid_atom_index_dict[list(amino_acid_atom_index_dict.keys())[-1]][-1])
    #amino_acid_atom_index_dict[list(amino_acid_atom_index_dict.keys())[-1]].pop(-1)
    #print(amino_acid_atom_index_dict[list(amino_acid_atom_index_dict.keys())[-1]][-1])


    for amino_idx in amino_acid_atom_index_dict.keys():
        amino_acid_ref = amino_acid_edge_list[amino_acid_type_dict[amino_idx]]

        if amino_idx > list(amino_acid_atom_index_dict.keys())[0]:
            # Add edges between amino acids
            protein_edge_list.append((amino_acid_atom_index_dict[amino_idx - 1][2], amino_acid_atom_index_dict[amino_idx][0], 1))

        # Add edges inside each amino acid
        for edge in amino_acid_ref:
            protein_edge_list.append((amino_acid_atom_index_dict[amino_idx][edge[0]], amino_acid_atom_index_dict[amino_idx][edge[1]], edge[2]))


    # Add the last oxygen
    protein_edge_list.append((amino_acid_atom_index_dict[list(amino_acid_atom_index_dict.keys())[-1]][2], len(atom_types) - 1, 1.5))

    for atom_idx in range(len(atom_types)):
        print(f"Atom {atom_idx}: Type={atom_types[atom_idx]}, Amino Acid={amino_acid_type_atom_list[atom_idx]}")


    exit() # TODO: Add hydrogen to the amino acids

    return {
        "atom_coords": atom_coords,
        "atom_symbols": atom_types,
        "covalent_edge_list": protein_edge_list,
        "amino_acid_type_belongings": amino_acid_type_atom_list,
    }

def plot_structure(coords, label):
    # Unpack coordinates and atom symbols, and also prepare text labels that include atom index
    x, y, z, symbols = zip(*[(atom[1], atom[2], atom[3], atom[0]) for atom in coords])
    text_labels = [f'{symbol} {idx}' for idx, symbol in enumerate(symbols)]  # Include atom index in label

    # Create a Plotly 3D scatter plot with both symbols and indices
    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers+text',  # Use both markers and text
        text=text_labels,  # Use text labels that include atom indices
        marker=dict(size=5, color='blue')
    )])
    fig.update_layout(
        title=f"Amino Acid Structure: {label}",
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z')
    )
    fig.show()

"""
GLY:
(0, 1, 1)
(1, 2, 1)
(2, 3, 2)

LYS:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(6, 7, 1)
(7, 8, 1)

ASP:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1.5)
(5, 7, 1.5)

VAL:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(3, 6, 1)

ASN:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(5, 7, 2)

GLU:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(6, 7, 1.5)
(6, 8, 1.5)

PHE:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1.5)
(5, 7, 1.5)
(6, 8, 1.5)
(7, 9, 1.5)
(9, 10, 1.5)
(8, 10, 1.5)

THR:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(3, 6, 1)

GLN:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(6, 7, 1)
(6, 8, 2)

ILE:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 7, 1)
(3, 6, 1)

HIS:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 7, 1.5)
(5, 6, 1.5)
(6, 9, 1.5)
(9, 8, 1.5)
(7, 8, 1.5)

LEU:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 7, 1)
(5, 6, 1)

SER:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)

CYS:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)

TRP:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 7, 1.5)
(5, 6, 1.5)
(6, 10, 1.5)
(10, 8, 1.5)
(8, 12, 1.5)
(7, 8, 1.5)
(7, 9, 1.5)
(9, 13, 1.5)
(13, 11, 1.5)
(12, 11, 1.5)

PRO:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(0, 6, 1)


ARG:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(6, 7, 1)
(7, 10, 1)
(10, 8, 1)
(10, 9, 2)

ALA:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)

MET:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 6, 1)
(6, 7, 1)

TYR:
(0, 1, 1)
(1, 2, 1)
(2, 4, 2)

(1, 3, 1)
(3, 5, 1)
(5, 7, 1.5)
(5, 6, 1.5)
(6, 8, 1.5)
(8, 11, 1.5)
(7, 9, 1.5)
(9, 11, 1.5)
(11, 10, 1)

"""


def plot_structure_electronegativity(atom_coords, atom_symbols, edge_list, amino_acid_type_belongings, title):
    x, y, z = zip(*atom_coords)  # Unpack the coordinates
    text_labels = [f'{symbol} {i} ({aa})' for i, (symbol, aa) in enumerate(zip(atom_symbols, amino_acid_type_belongings))]  # Include amino acid type in label

    # Electronegativity values for different atoms (not directly used for colors anymore)
    electronegativities = {
        "C": 2.55,
        "N": 3.04,
        "O": 3.44,
        "S": 2.58,
        "P": 2.19,
    }

    # Define a function to map atom symbols to colors
    def get_color(symbol):
        color_map = {
            "C": 'rgb(128, 128, 128)',  # Gray
            "N": 'rgb(255, 0, 0)',      # Red
            "O": 'rgb(0, 0, 255)',      # Blue
            "S": 'rgb(255, 165, 0)',    # Orange
            "P": 'rgb(255, 255, 0)'     # Yellow for phosphorus
        }
        return color_map.get(symbol, 'rgb(255, 255, 255)')  # Default white if not specified

    # Map atom symbols to colors
    colors = [get_color(symbol) for symbol in atom_symbols]

    # Create scatter plot for atoms
    atom_trace = go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers+text',
        text=text_labels,
        marker=dict(size=5, color=colors),
        name='Atoms'
    )

    # Create line plot for bonds
    bond_traces = []
    for start, end, bond_type in edge_list:
        # Assign colors based on bond type
        bond_color_map = {
            1: 'green',    # Single bond
            1.5: 'red',    # Aromatic bond
            2: 'blue',     # Double bond
            3: 'black'     # Triple bond (added for completeness)
        }
        color = bond_color_map.get(bond_type, 'grey')  # Default grey for other unspecified bond types

        bond_trace = go.Scatter3d(
            x=[atom_coords[start][0], atom_coords[end][0]],
            y=[atom_coords[start][1], atom_coords[end][1]],
            z=[atom_coords[start][2], atom_coords[end][2]],
            mode='lines',
            line=dict(color=color, width=2),
            name=f'Bond {bond_type}'
        )
        bond_traces.append(bond_trace)

    # Combine atom and bond traces
    fig = go.Figure(data=[atom_trace, *bond_traces])
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        legend_title_text='Element'
    )
    fig.show()

def plot_structure_standard(atom_coords, atom_symbols, edge_list, amino_acid_type_belongings, title):
    x, y, z = zip(*atom_coords)  # Unpack the coordinates
    text_labels = [f'{symbol} {i} ({aa})' for i, (symbol, aa) in enumerate(zip(atom_symbols, amino_acid_type_belongings))]  # Include amino acid type in label

    # Define a color map for amino acid types
    color_map = {
        'ALA': 'blue', 'CYS': 'yellow', 'ASP': 'red', 'GLU': 'orange',
        'PHE': 'green', 'GLY': 'pink', 'HIS': 'gray', 'ILE': 'cyan',
        'LYS': 'magenta', 'LEU': 'brown', 'MET': 'olive', 'ASN': 'lime',
        'PRO': 'teal', 'GLN': 'purple', 'ARG': 'maroon', 'SER': 'navy',
        'THR': 'chocolate', 'VAL': 'gold', 'TRP': 'violet', 'TYR': 'black'
    }
    # Default color for unknown or undefined amino acids
    default_color = 'white'

    # Map amino acid types to colors
    colors = [color_map.get(aa, default_color) for aa in amino_acid_type_belongings]

    # Create scatter plot for atoms
    atom_trace = go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers+text',
        text=text_labels,
        marker=dict(size=5, color=colors),
        name='Atoms'
    )

    # Create line plot for bonds
    bond_traces = []
    for start, end, bond_type in edge_list:
        bond_trace = go.Scatter3d(
            x=[atom_coords[start][0], atom_coords[end][0]],
            y=[atom_coords[start][1], atom_coords[end][1]],
            z=[atom_coords[start][2], atom_coords[end][2]],
            mode='lines',
            line=dict(color='green' if bond_type == 1 else 'blue' if bond_type == 2 else 'red', width=2),
            name=f'Bond {bond_type}'
        )
        bond_traces.append(bond_trace)

    # Combine atom and bond traces
    fig = go.Figure(data=[atom_trace, *bond_traces])
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        legend_title_text='Element'
    )
    fig.show()

if __name__ == "__main__":
    af_uniprots = set()
    for file in os.listdir("AFDATA"):
        if file.endswith(".pdb.gz"):
            uniprot_id = file.split("-")[1]
            af_uniprots.add(uniprot_id)
            
    uniprots = sorted(list(af_uniprots))

    for uniprot_id in tqdm(uniprots):
        uniprot_id = "Q15109"
        protein_data = get_protein_data(uniprot_id)

        atom_coords = protein_data["atom_coords"]
        atom_symbols = protein_data["atom_symbols"]
        edge_list = protein_data["covalent_edge_list"]
        amino_acid_type_belongings = protein_data["amino_acid_type_belongings"]

        # Plot the protein structure
        plot_structure_standard(atom_coords, atom_symbols, edge_list, amino_acid_type_belongings, f"Protein Structure for {uniprot_id}")
        plot_structure_electronegativity(atom_coords, atom_symbols, edge_list, amino_acid_type_belongings, f"Protein Structure for {uniprot_id}")
        print(uniprot_id)
        # Show only once
        break



