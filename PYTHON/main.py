import python_tower as pt
import middleware as m

if __name__ == "__main__":

    target_model = lambda id: {
        "id": id,
        "dynamic": "true",
        "point_cloud": "here_goes_a_serialized_numpy_array_with_atom_types_and_atom_cords",
        "molecule_graph": "here_goes_a_cso_graph_or_something",
    }

    molecule_model = lambda id: {
        "id": id,
        "dynamic": "true",
        "point_cloud": "here_goes_a_serialized_numpy_array_with_atom_types_and_atom_cords",
        "molecule_graph": "here_goes_a_cso_graph_or_something",
    }

    # Define the data you want to send
    data = {
        "target": m.get_protein_data("P51659", path="PYTHON/"),
        "molecules": [
            m.get_molecule_data(smiles="CCO"),
            m.get_molecule_data(smiles="C1=CC=CC=C1")
        ]
    }

    # Send the data and wait for the response
    r = pt.spin(data)

    print(r)