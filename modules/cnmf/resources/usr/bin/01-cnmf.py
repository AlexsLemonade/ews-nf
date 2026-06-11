#!/usr/bin/env python3

import argparse
import multiprocessing
from pathlib import Path

import numpy
import anndata
import pandas 
import cnmf

import scipy.sparse


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for cNMF processing."""
    parser = argparse.ArgumentParser(description="Run cNMF on a given dataset.")
    parser.add_argument(
        "--h5ad_file",
        type=Path,
        required=True,
        help="Path to the H5AD file containing a processed AnnData object",
    )
    parser.add_argument(
        "--annotations_file",
        type=Path,
        required=True,
        help="Path to the TSV file containing the tumor cell annotations",
    )
    parser.add_argument(
        "--annotations_column",
        type=str,
        required=True,
        help="Column name in the annotations file containing the cell type annotations to use for cNMF.",
    )
    parser.add_argument(
        "--celltype_value",
        type=str,
        default="tumor",
        help="Value in the annotations column that corresponds to the cell type to subset for cNMF",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=Path,
        required=True,
        help="Path to the directory for output and intermediate results.",
    )
    parser.add_argument(
        "--unique_id",
        type=str,
        required=True,
        help="Unique ID of the sample to process. Will be created as a subdirectory of the output directory.",
    )
    parser.add_argument(
        "--k_components_lower",
        type=int,
        default=5,
        help="Lower boundary for number of components (k) for cNMF. (default: %(default)d)",
    )
    parser.add_argument(
        "--k_components_upper",
        type=int,
        default=20,
        help="Upper boundary for number of components (k) for cNMF. (default: %(default)d)",
    )
    parser.add_argument(
        "--k_step_size",
        type = int,
        default = 1,
        help="Step size to use for creating the range of k components to use when running cNMF. (default: %(default)d)",
    )
    parser.add_argument(
        "--n-iter",
        "-n",
        type=int,
        default=100,
        help="Number of iterations for cNMF. (default: %(default)d)",
    )
    parser.add_argument(
        "--max-nmf-iter",
        "-m",
        type=int,
        default=2000,
        help="Maximum number of iterations for NMF. (default: %(default)d)",
    )
    parser.add_argument(
        "--density-threshold",
        type=float,
        default=0.1,
        help="Density threshold for consensus. (default: %(default)f)",
    )
    parser.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=multiprocessing.cpu_count(),
        help="Number of parallel jobs to run. (default: %(default)d)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="Random seed for reproducibility. (default: %(default)d)",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    # Validate input files
    if not args.h5ad_file.is_file() and args.h5ad_file.suffix == ".h5ad":
        raise ValueError(
            f"H5AD file not found at {args.h5ad_file}. Please ensure the file exists."
        )
    
    if not args.annotations_file.is_file() and args.annotations_file.suffix == ".tsv":
        raise ValueError(
            f"Annotations file not found at {args.annotations_file}. Please ensure the file exists."
        )
    
    # check that k_upper and k_lower are numbers and k_lower is < k_upper
    if args.k_components_lower >= args.k_components_upper:
        raise ValueError(
            "k_components_lower must be less than k_components_upper"
        )

    # assign k range to use for cNMF 
    k_range = numpy.arange(args.k_components_lower, args.k_components_upper + 1, args.k_step_size)

    # read in h5ad
    adata = anndata.read_h5ad(args.h5ad_file)

    # read in annotations tsv
    annotations_df = pandas.read_csv(args.annotations_file, sep= "\t")
    # make sure the specified annotations column is provided
    if not args.annotations_column in annotations_df.columns:
        raise KeyError(
            f"Provided annotations file does not contain the required column, `{args.annotations_column}`"
        )
    
    # get vector of tumor cells to keep 
    tumor_barcodes = annotations_df.loc[
        annotations_df[args.annotations_column] == args.celltype_value, "barcodes"
    ].tolist()

    # make sure .X contains the raw counts 
    # if adata.raw is present, then .X is normalized so need to replace it with the raw counts from adata.raw.X
    if adata.raw is not None:
        adata.X = adata.raw.X

    # check that adata.X has rounded counts and not normalized counts
    x_vals = adata.X.data if scipy.sparse.issparse(adata.X) else adata.X.ravel()
    if not numpy.all(x_vals == numpy.floor(x_vals)):
        raise ValueError(f"adata.X does not contain integer counts")
    
    # create a new anndata object with just tumor cells 
    subset_adata = adata[tumor_barcodes, ].copy()

    # get the mean of each gene across all cells and only keep those with expression > 0
    detected_genes = subset_adata.X.mean(axis=0) > 0
    subset_adata = adata[tumor_barcodes, detected_genes].copy()

    # make a simplified version of the object for easy exporting
    new_adata = anndata.AnnData(X=subset_adata.X)
    new_adata.obs_names = subset_adata.obs_names
    new_adata.var_names = subset_adata.var_names

    # save counts to h5 file to use as input to cnmf
    anndata_file = "anndata.h5ad"
    anndata.AnnData.write_h5ad(new_adata, anndata_file)

    # Prepare cNMF object
    cnmf_obj = cnmf.cNMF(output_dir=args.output_dir, name=args.unique_id)
    cnmf_obj.prepare(
        counts_fn=anndata_file,
        components=k_range,
        n_iter=args.n_iter,
        max_NMF_iter=args.max_nmf_iter,
        seed=args.seed,
    )

    # Run cNMF factorization
    cnmf_obj.factorize_multi_process(args.jobs)
    cnmf_obj.combine()  # combine results from all runs

    # get the consensus for each value of k
    for k in k_range:
        try:
            cnmf_obj.consensus(k=k, density_threshold=args.density_threshold)
        except RuntimeError as e:
            print(f"WARNING: consensus failed for k={k}: {e}")


if __name__ == "__main__":
    main()
    