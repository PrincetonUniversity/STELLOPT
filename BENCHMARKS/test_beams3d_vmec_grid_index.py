#!/usr/bin/env python3
from pathlib import Path

import numpy as np


source = (
    Path(__file__).resolve().parent.parent
    / "BEAMS3D/Sources/beams3d_init_vmec.f90"
).read_text()
helper = source[
    source.index("SUBROUTINE beams3d_vmec_grid_index"):source.index(
        "END SUBROUTINE beams3d_vmec_grid_index"
    )
]
compact = "".join(helper.lower().split())
assert "i=mod(s-1,nr)+1" in compact
assert "j=mod(s-1,nr*nphi)/nr+1" in compact
assert "k=(s-1)/(nr*nphi)+1" in compact
assert source.count("CALL beams3d_vmec_grid_index(s,i,j,k)") == 7
assert "CEILING(REAL" not in source
assert "FLOOR(REAL" not in source


def integer_index(s, nr, nphi):
    return (
        (s - 1) % nr + 1,
        ((s - 1) % (nr * nphi)) // nr + 1,
        (s - 1) // (nr * nphi) + 1,
    )


def old_plane(s, nr, nphi):
    return int(np.ceil(np.float32(s) / np.float32(nr * nphi)))


for nr, nphi, nz in ((64, 32, 64), (128, 64, 128), (256, 128, 256)):
    plane = nr * nphi
    samples = {1, nr, nr + 1, plane, plane + 1, nr * nphi * nz}
    samples.update((k - 1) * plane + 1 for k in range(1, nz + 1))
    for s in samples:
        i, j, k = integer_index(s, nr, nphi)
        assert s == i + (j - 1) * nr + (k - 1) * plane
        assert k == old_plane(s, nr, nphi)

nr, nphi, nz = 384, 192, 384
plane = nr * nphi
defects = []
for k in range(1, nz + 1):
    s = (k - 1) * plane + 1
    i_new, j_new, k_new = integer_index(s, nr, nphi)
    assert (i_new, j_new, k_new) == (1, 1, k)
    if old_plane(s, nr, nphi) != k:
        defects.append(s)

assert len(defects) == 156
print(f"384-grid plane-boundary defects in float decoder: {len(defects)}")
