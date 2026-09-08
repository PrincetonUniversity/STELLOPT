#!/usr/bin/env python3
from argparse import ArgumentParser

import h5py
import numpy as np


DATASETS = ("S_ARR", "U_ARR", "B_R", "B_PHI", "B_Z", "S_lines", "B_lines")


def compare(reference_path, candidate_path):
    failed = False
    with h5py.File(reference_path) as reference, h5py.File(candidate_path) as candidate:
        reference_mask = reference["S_ARR"][:] == 4.0
        candidate_mask = candidate["S_ARR"][:] == 4.0
        mask_changes = np.count_nonzero(reference_mask != candidate_mask)
        print(f"default-mask changes: {mask_changes}")
        failed = mask_changes != 0
        for name in DATASETS:
            reference_values = reference[name][:]
            candidate_values = candidate[name][:]
            changes = np.count_nonzero(reference_values != candidate_values)
            max_difference = np.max(np.abs(reference_values - candidate_values))
            print(f"{name}: changed={changes} max_abs={max_difference:.17e}")
            failed = failed or changes != 0
    return failed


parser = ArgumentParser()
parser.add_argument("reference")
parser.add_argument("candidates", nargs="+")
args = parser.parse_args()
failed = False
for candidate in args.candidates:
    failed = compare(args.reference, candidate) or failed
raise SystemExit(failed)
