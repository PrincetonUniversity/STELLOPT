#!/usr/bin/env python3
from argparse import ArgumentParser

import h5py
import numpy as np


def support_layers(flux):
    interior = (flux >= 0.0) & (flux <= 1.0)
    exterior = flux > 1.0
    layer1 = np.zeros_like(interior)
    layer2 = np.zeros_like(interior)

    layer1[2:] |= interior[1:-1] & interior[:-2]
    layer1[:-2] |= interior[1:-1] & interior[2:]
    layer1[:, :, 2:] |= interior[:, :, 1:-1] & interior[:, :, :-2]
    layer1[:, :, :-2] |= interior[:, :, 1:-1] & interior[:, :, 2:]

    layer2[3:] |= exterior[2:-1] & interior[1:-2] & interior[:-3]
    layer2[:-3] |= exterior[1:-2] & interior[2:-1] & interior[3:]
    layer2[:, :, 3:] |= (
        exterior[:, :, 2:-1] & interior[:, :, 1:-2] & interior[:, :, :-3]
    )
    layer2[:, :, :-3] |= (
        exterior[:, :, 1:-2] & interior[:, :, 2:-1] & interior[:, :, 3:]
    )
    return layer1 & exterior, layer2 & exterior


parser = ArgumentParser()
parser.add_argument("reference")
parser.add_argument("candidate")
args = parser.parse_args()

failed = False
with h5py.File(args.reference) as reference, h5py.File(args.candidate) as candidate:
    for name in ("S_ARR", "U_ARR", "wall_vertex", "wall_faces"):
        changed = np.count_nonzero(reference[name][:] != candidate[name][:])
        print(f"{name}: changed={changed}")
        failed |= changed != 0

    flux = reference["S_ARR"][:]
    interior = (flux >= 0.0) & (flux <= 1.0)
    layer1, layer2 = support_layers(flux)
    changed = np.zeros_like(interior)
    for name in ("B_R", "B_PHI", "B_Z"):
        reference_values = reference[name][:]
        candidate_values = candidate[name][:]
        component_changes = reference_values != candidate_values
        interior_changes = np.count_nonzero(component_changes & interior)
        nonfinite = np.count_nonzero(~np.isfinite(candidate_values))
        print(f"{name}: interior_changed={interior_changes} nonfinite={nonfinite}")
        changed |= component_changes
        failed |= interior_changes != 0 or nonfinite != 0

    layer1_changes = np.count_nonzero(changed & layer1)
    layer2_only_changes = np.count_nonzero(changed & layer2 & ~layer1)
    unsupported_changes = np.count_nonzero(changed & ~(layer1 | layer2))
    print(
        f"support_changed={np.count_nonzero(changed)} "
        f"layer1={layer1_changes} layer2_only={layer2_only_changes} "
        f"unsupported={unsupported_changes}"
    )
    failed |= layer1_changes == 0 or layer2_only_changes == 0
    failed |= unsupported_changes != 0

raise SystemExit(failed)
