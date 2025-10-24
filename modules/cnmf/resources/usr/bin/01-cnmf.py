#!/usr/bin/env python3

import argparse
import multiprocessing
from pathlib import Path

import anndata
import numpy
import cnmf


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
        default=15,
        help="Upper boundary for number of components (k) for cNMF. (default: %(default)d)",
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

    # Validate input file
    if not args.h5ad_file.is_file() and args.h5ad_file.suffix == ".h5ad":
        raise ValueError(
            f"H5AD file not found at {args.h5ad_file}. Please ensure the file exists."
        )
    
    # check that k_upper and k_lower are numbers and k_lower is < k_upper
    if args.k_components_lower >= args.k_components_upper:
        raise ValueError(
            "k_components_lower must be less than k_components_upper"
        )

    # assign k range to use for cNMF 
    k_range = numpy.arange(args.k_components_lower, args.k_components_upper + 1)

    # read in h5ad
    adata = anndata.read_h5ad(args.h5ad_file)

    # create a new anndata object with just the raw counts for hvgs
    hv_genes = adata.var['highly_variable']
    adata.X = adata.raw.X # make sure X has the raw counts for subsetting
    subset_adata = adata[:, hv_genes].copy()
    new_adata = anndata.AnnData(X=subset_adata.X)

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
        cnmf_obj.consensus(k=k, density_threshold=args.density_threshold)


if __name__ == "__main__":
    main()