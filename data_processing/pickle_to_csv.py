"""
pickle_to_csv.py

Reads all .pkl files in a directory, extracts 'weights' (2D numpy array),
labels rows with 'sources' and columns with 'targets', and writes each to a
.csv file in a 'csv_output' subdirectory.

Usage:
    python pickle_to_csv.py <input_directory>
    python pickle_to_csv.py <input_directory> --output <output_subdirectory>
"""

import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd


def load_pickle(filepath: str) -> dict:
    with open(filepath, "rb") as f:
        return pickle.load(f)


def extract_and_label(data: dict, filepath: str) -> pd.DataFrame:
    """Extract weights from a dict and label rows/cols with sources/targets."""
    for key in ("weights", "sources", "targets"):
        if key not in data:
            raise KeyError(f"Missing required key '{key}' in {filepath}")

    weights = data["weights"]
    sources = data["sources"]
    targets = data["targets"]

    if not isinstance(weights, np.ndarray) or weights.ndim != 2:
        raise ValueError(
            f"'weights' in {filepath} must be a 2D numpy array, "
            f"got {type(weights).__name__} with shape {getattr(weights, 'shape', '?')}"
        )

    n_rows, n_cols = weights.shape

    if len(sources) != n_rows:
        raise ValueError(
            f"'sources' length ({len(sources)}) does not match "
            f"weights row count ({n_rows}) in {filepath}"
        )
    if len(targets) != n_cols:
        raise ValueError(
            f"'targets' length ({len(targets)}) does not match "
            f"weights column count ({n_cols}) in {filepath}"
        )

    return pd.DataFrame(weights, index=sources, columns=targets)


def process_directory(input_dir: str, output_subdir: str = "csv_output") -> None:
    if not os.path.isdir(input_dir):
        print(f"Error: '{input_dir}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    output_dir = os.path.join(input_dir, output_subdir)
    os.makedirs(output_dir, exist_ok=True)

    pickle_files = sorted(
        f for f in os.listdir(input_dir)
        if f.lower().endswith(".pkl") or f.lower().endswith(".pickle")
    )

    if not pickle_files:
        print(f"No .pkl / .pickle files found in '{input_dir}'.")
        return

    print(f"Found {len(pickle_files)} pickle file(s). Writing CSVs to '{output_dir}/'\n")

    success, failed = 0, 0

    for filename in pickle_files:
        src_path = os.path.join(input_dir, filename)
        stem = os.path.splitext(filename)[0]
        dst_path = os.path.join(output_dir, stem + ".csv")

        try:
            data = load_pickle(src_path)
            df = extract_and_label(data, src_path)
            df.to_csv(dst_path)
            print(f"  [OK]  {filename}  →  {stem}.csv  (shape {df.shape})")
            success += 1
        except Exception as exc:
            print(f"  [FAIL] {filename}: {exc}", file=sys.stderr)
            failed += 1

    print(f"\nDone. {success} succeeded, {failed} failed.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert pickle weight matrices to labelled CSV files."
    )
    parser.add_argument(
        "input_dir",
        help="Directory containing the .pkl files to process.",
    )
    parser.add_argument(
        "--output",
        default="csv_output",
        metavar="SUBDIR",
        help="Name of the output subdirectory (default: csv_output).",
    )
    args = parser.parse_args()
    process_directory(args.input_dir, args.output)


if __name__ == "__main__":
    main()