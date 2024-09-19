import argparse
import scanpy as sc
import pandas as pd

def main(args: argparse.Namespace) -> None:
    # Load the 10X Genomics data from the provided matrix directory
    adata = sc.read_10x_mtx(args.mtx)
    print(f'Loaded {args.mtx}')
    # If barcodes file is provided, filter the AnnData object
    if args.barcodes:
        selected_barcodes = pd.read_csv(args.barcodes, header=None)[0].tolist()
        adata = adata[adata.obs_names.isin(selected_barcodes)]
    
    # Write the AnnData object to the specified output directory
    adata.write_h5ad(args.outdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Setup raw data for integration model setup')
    
    parser.add_argument('mtx', metavar='MTX_DIRECTORY', type=str, 
                        help='Path to the directory containing the 10X Genomics matrix (.mtx) files')
    parser.add_argument('-b', '--barcodes', type=str, default=None, 
                        help='Path to a CSV file containing barcodes to filter the data')
    parser.add_argument('-O', '--outdir', type=str, default='./output.h5ad', 
                        help='Path to the output .h5ad file')

    args = parser.parse_args()
    main(args)