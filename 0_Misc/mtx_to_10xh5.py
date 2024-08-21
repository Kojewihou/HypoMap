from scipy.sparse import csc_array
import h5py
import numpy as np
import sys
import argparse

def main (args):
    matrix_name = args.matrix_name if args.matrix_name else "matrix.mtx"
    matrix_file = f"{args.matrix_path}/{matrix_name}"
    barcodes_file = f"{args.matrix_path}/barcodes.tsv"
    features_file = f"{args.matrix_path}/features.tsv"


    # Unpack matrix into numpy arrays
    with open(matrix_file, 'r') as file:
        third_line = file.readlines()[2]
    shape = np.array([int(value) for value in third_line.split()])

    # Load coordinates and data from matrix file
    rows, cols, values = np.loadtxt(matrix_file, unpack=True, skiprows=3, dtype='i')

    # Convert from 1-base coordinates to 0-base indices
    rows = rows - 1
    cols = cols -1

    # Unpack barcodes into numpy array
    barcodes=np.loadtxt(barcodes_file, unpack=True, dtype='S256')

    # Unpack features into numpy array
    ids, names, feature_type = np.loadtxt(features_file, unpack=True, dtype='S256', delimiter='\t')

    # Generate Genome array
    genome = np.full_like(ids, "GRCm39", dtype='S256')

    sparse_matrix=csc_array((values,(rows, cols)),(shape[0],shape[1]))

    with h5py.File(args.output_file, 'w') as f:

        # Create matrix group
        matrix_group = f.create_group('matrix')

        # Create datasets for rows, cols, and data
        matrix_group.create_dataset('barcodes', data=barcodes, dtype="S256", maxshape=(None,))
        matrix_group.create_dataset('data', data=sparse_matrix.data, dtype='i', maxshape=(None,))
        matrix_group.create_dataset('indices', data=sparse_matrix.indices, dtype='i', maxshape=(None,))
        matrix_group.create_dataset('indptr', data=sparse_matrix.indptr, dtype='i', maxshape=(None,))

        # Create shape dataset
        matrix_group.create_dataset('shape', data=sparse_matrix.shape, dtype='I', maxshape=(None,))

        # Create features group
        features_group = matrix_group.create_group('features')

        # Add relevant datasets and attributes
        features_group.create_dataset('feature_type', data=feature_type, dtype='S256', maxshape=(None,))
        features_group.create_dataset('genome', data=genome, dtype='S256', maxshape=(None,))
        features_group.create_dataset('id', data=ids, dtype='S256', maxshape=(None,))
        features_group.create_dataset('name', data=names, dtype='S256', maxshape=(None,))

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Your script description here.")

    parser.add_argument("matrix_path", help="Path to the matrix files")
    parser.add_argument("output_file", help="File where the output h5 be stored")
    parser.add_argument("-m", "--matrix_name", help="Name of the matrix file (default: matrix.mtx)")

    args = parser.parse_args()
    
    main(args)

    sys.exit(0)