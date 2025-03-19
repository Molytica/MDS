import middleware as m
import numpy as np
import python_tower as pt

def dock(smiles, uniprot_id):

    prot_voxel = get_protein_voxel(uniprot_id)

    filter_voxels = get_filter_voxels(smiles)

    filter_voxel_outputs = []
    for filter in filter_voxels:
        filter_voxel_outputs.append(convolute(prot_voxel, filter))


def get_protein_voxel(uniprot_id):
    atom_coords = m.get_protein_data(uniprot_id)["atom_coords"]

    offset =np.array([20, 20, 20])
    
    np_atom_coords = np.array(atom_coords)
    min_coords = np.min(np_atom_coords, axis=0)
    max_coords = np.max(np_atom_coords, axis=0)

    # Shift the coordinates by subtracting the min_coords
    np_atom_coords = np_atom_coords - min_coords + offset / 2

    min_coords = np.min(np_atom_coords, axis=0)
    max_coords = np.max(np_atom_coords, axis=0)

    print(min_coords, max_coords)

    # Create a voxel with the shape of max_coords
    voxel = np.zeros(np.round(max_coords).astype(int) + np.round(offset / 2).astype(int))

    # Fill the voxel with the atom coordinates
    for atom_coord in np_atom_coords:
        voxel[np.round(atom_coord).astype(int)] += 1

    print(voxel)
    print(np.average(voxel))
    print(f"Voxel shape: {voxel.shape}")
    print(f"Voxel non-zero elements: {np.count_nonzero(voxel)}")
    #visualize_voxel(voxel)

def visualize_voxel(voxel):
    pt.spin({
        "voxel": voxel.tolist()
    })

def get_filter_voxels(smiles):
    pass

def convolute(prot_voxel, filter_voxel):
    pass

if __name__ == "__main__":
    dock("CCO", "A0A024R1R8")

