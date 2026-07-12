#!/usr/bin/env python3

import argparse
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np


def event_data(path):
    with h5py.File(path, "r") as data:
        return {name: data[name][:] for name in (
            "B_lines", "PHI_lines", "R_lines", "S_lines", "U_lines",
            "Z_lines", "end_state", "lplasma_only", "lwall_from_vmec", "mass",
            "moment_lines", "t_end",
            "time_lines", "vll_lines", "vr_lines", "vphi_lines", "vz_lines",
            "wall_hit_b", "wall_hit_energy", "wall_hit_field_valid",
            "wall_hit_fraction", "wall_hit_model", "wall_hit_moment",
            "wall_hit_phi", "wall_hit_r", "wall_hit_s", "wall_hit_time",
            "wall_hit_u", "wall_hit_z",
            "wall_hit_valid", "wall_hit_vll", "wall_hit_vr",
            "wall_hit_vphi", "wall_hit_vz",
        )}


def run_case(source, executable, launcher, name, transform, hitonly=False):
    with tempfile.TemporaryDirectory(prefix=f"beams3d-{name}-") as temporary:
        run_dir = Path(temporary)
        text = transform((source / "input.ORBITS_loss").read_text())
        (run_dir / "input.ORBITS_loss").write_text(text)
        shutil.copy2(source / "wout_ORBITS_loss.nc", run_dir / "wout_ORBITS_loss.nc")
        command = [*launcher, str(executable), "-vmec", "ORBITS_loss", "-plasma"]
        if hitonly:
            command.append("-hitonly")
        result = subprocess.run(command, cwd=run_dir, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            output = "\n".join((result.stdout + result.stderr).splitlines()[-80:])
            raise RuntimeError(f"{name} failed with status {result.returncode}\n{output}")
        return event_data(run_dir / "beams3d_ORBITS_loss.h5")


def add_full_orbit(text):
    anchor = "  FOLLOW_TOL =  1.00000000000000E-08\n"
    if text.count(anchor) != 1:
        raise ValueError("unexpected FOLLOW_TOL assignment")
    settings = (
        "  RHO_FULLORBIT = 0.0\n"
        "  VR_START_IN = 40*5.0E7\n"
        "  VPHI_START_IN = 40*0.0\n"
        "  VZ_START_IN = 40*0.0\n"
    )
    text = text.replace(anchor, anchor + settings)
    start = "  R_START_IN = 40*10.85"
    if text.count(start) != 1:
        raise ValueError("unexpected R_START_IN assignment")
    text = text.replace(start, "  R_START_IN = 40*10.98")
    rmax = "  RMAX =  11.00"
    if text.count(rmax) != 1:
        raise ValueError("unexpected RMAX assignment")
    text = text.replace(rmax, "  RMAX =  11.10")
    return set_trace_time(text, "1.0E-7")


def add_full_orbit_outside(text):
    text = add_full_orbit(text).replace("  RMAX =  11.10", "  RMAX =  11.00")
    return text.replace("  RMIN =   9.00", "  RMIN =  10.90")


def use_rkh68(text):
    old = "  INT_TYPE = 'LSODE'"
    if text.count(old) != 1:
        raise ValueError("unexpected INT_TYPE assignment")
    return set_trace_time(text.replace(old, "  INT_TYPE = 'RKH68'"), "2.0E-4")


def set_trace_time(text, value):
    old = "  T_END_IN = 40*1.0E-3"
    if text.count(old) != 1:
        raise ValueError("unexpected T_END_IN assignment")
    return text.replace(old, f"  T_END_IN = 40*{value}")


def write_outside_plane(path):
    near = 11.0000000003
    far = 11.0000000008
    vertices = np.array(((near, -1.0, -1.0), (near, 1.0, -1.0),
                         (near, 1.0, 1.0), (near, -1.0, 1.0),
                         (far, -1.0, -1.0), (far, 1.0, -1.0),
                         (far, 1.0, 1.0), (far, -1.0, 1.0)))
    with path.open("w") as handle:
        handle.write("MACHINE: ORBITS plane boundary\nDATE: regression fixture\n")
        handle.write("8 4\n")
        np.savetxt(handle, vertices, fmt="%.17e")
        handle.write("1 2 3\n1 3 4\n5 6 7\n5 7 8\n")


def run_external_wall(source, executable, launcher):
    with tempfile.TemporaryDirectory(prefix="beams3d-full-orbit-outside-") as temporary:
        run_dir = Path(temporary)
        (run_dir / "input.ORBITS_loss").write_text(
            add_full_orbit_outside((source / "input.ORBITS_loss").read_text())
        )
        shutil.copy2(source / "wout_ORBITS_loss.nc", run_dir / "wout_ORBITS_loss.nc")
        write_outside_plane(run_dir / "wall.dat")
        command = [*launcher, str(executable), "-vmec", "ORBITS_loss", "-plasma",
                   "-vessel", "wall.dat"]
        result = subprocess.run(command, cwd=run_dir, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            output = "\n".join((result.stdout + result.stderr).splitlines()[-80:])
            raise RuntimeError(f"full-orbit-outside failed with status {result.returncode}\n{output}")
        return event_data(run_dir / "beams3d_ORBITS_loss.h5")


def check_full_orbit(data, wall_from_vmec=True):
    hit = data["wall_hit_valid"] == 1
    if not np.any(hit):
        raise AssertionError(
            f"full-orbit case recorded no hits; max R={np.max(data['R_lines']):.16e}, "
            f"max s={np.max(data['S_lines']):.16e}, states={np.unique(data['end_state'])}"
        )
    assert (data["lwall_from_vmec"][0] == 1) == wall_from_vmec
    assert data["lplasma_only"][0] == 1
    assert np.all(data["wall_hit_model"][hit] == 2)
    field_hit = hit & (data["wall_hit_field_valid"] == 1)
    for name in ("wall_hit_vr", "wall_hit_vphi", "wall_hit_vz"):
        assert np.all(data[name][~hit] < -1.0e300)
    energy = 0.5*data["mass"][hit]*(
        data["wall_hit_vr"][hit]**2 + data["wall_hit_vphi"][hit]**2 +
        data["wall_hit_vz"][hit]**2
    )
    assert np.array_equal(data["wall_hit_energy"][hit], energy)
    populated = (data["R_lines"][hit] != 0.0) | (data["Z_lines"][hit] != 0.0)
    last = populated.shape[1] - 1 - np.argmax(populated[:, ::-1], axis=1)
    rows = np.flatnonzero(hit)
    for event_name, line_name in (
        ("wall_hit_vr", "vr_lines"),
        ("wall_hit_vphi", "vphi_lines"),
        ("wall_hit_vz", "vz_lines"),
    ):
        assert np.array_equal(data[event_name][hit], data[line_name][rows, last])
    if np.any(field_hit):
        field_rows = np.flatnonzero(field_hit)
        field_last = (populated.shape[1] - 1 -
                      np.argmax(((data["R_lines"][field_hit] != 0.0) |
                                 (data["Z_lines"][field_hit] != 0.0))[:, ::-1], axis=1))
        for event_name, line_name in (
            ("wall_hit_vll", "vll_lines"),
            ("wall_hit_moment", "moment_lines"),
            ("wall_hit_b", "B_lines"),
        ):
            assert np.array_equal(data[event_name][field_hit],
                                  data[line_name][field_rows, field_last])
        projected = (
            data["wall_hit_vll"][field_hit]**2 +
            2.0*data["wall_hit_moment"][field_hit]*data["wall_hit_b"][field_hit] /
            data["mass"][field_hit]
        )
        speed = (data["wall_hit_vr"][field_hit]**2 +
                 data["wall_hit_vphi"][field_hit]**2 +
                 data["wall_hit_vz"][field_hit]**2)
        assert np.allclose(projected, speed, rtol=1.0e-12, atol=0.0)
    for name in ("wall_hit_vll", "wall_hit_moment", "wall_hit_b",
                 "wall_hit_s", "wall_hit_u"):
        assert np.all(data[name][hit & ~field_hit] < -1.0e300)


def check_times(data):
    hit = data["wall_hit_valid"] == 1
    assert np.any(hit)
    assert np.all(data["wall_hit_time"][hit] >= 0.0)
    assert np.all(data["wall_hit_time"][hit] <= data["t_end"][hit])
    populated = (data["R_lines"] != 0.0) | (data["Z_lines"] != 0.0)
    assert np.all(data["time_lines"][populated] >= 0.0)
    for particle in range(data["time_lines"].shape[0]):
        times = data["time_lines"][particle, populated[particle]]
        assert np.all(np.diff(times) >= 0.0)


def check_hitonly(data):
    populated = (data["R_lines"] != 0.0) | (data["Z_lines"] != 0.0)
    assert np.all(data["time_lines"][populated] >= 0.0)
    for particle in range(data["time_lines"].shape[0]):
        times = data["time_lines"][particle, populated[particle]]
        assert np.all(np.diff(times) >= 0.0)
    for particle in np.flatnonzero(data["wall_hit_valid"] == 1):
        slots = np.flatnonzero(populated[particle])
        assert np.array_equal(slots, np.array([0, 1, 2]))
        fraction = data["wall_hit_fraction"][particle]
        times = data["time_lines"][particle, slots]
        assert times[1] == times[0] + fraction*(times[2]-times[0])
        xyz = np.column_stack((
            data["R_lines"][particle, slots]*np.cos(data["PHI_lines"][particle, slots]),
            data["R_lines"][particle, slots]*np.sin(data["PHI_lines"][particle, slots]),
            data["Z_lines"][particle, slots],
        ))
        assert np.allclose(xyz[1], xyz[0] + fraction*(xyz[2]-xyz[0]),
                           rtol=0.0, atol=32*np.finfo(float).eps)
        assert data["wall_hit_r"][particle] == data["R_lines"][particle, 1]
        assert data["wall_hit_phi"][particle] == data["PHI_lines"][particle, 1]
        assert data["wall_hit_z"][particle] == data["Z_lines"][particle, 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--launcher", default="")
    parser.add_argument("--case", choices=("all", "full-orbit", "rkh68", "hitonly"),
                        default="all")
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    launcher = shlex.split(args.launcher)
    if args.case in ("all", "rkh68"):
        rkh68 = run_case(source, args.executable.resolve(), launcher,
                         "rkh68", use_rkh68)
        check_times(rkh68)
    if args.case in ("all", "hitonly"):
        hitonly = run_case(source, args.executable.resolve(), launcher,
                           "hitonly", lambda text: set_trace_time(text, "2.0E-4"), hitonly=True)
        check_times(hitonly)
        check_hitonly(hitonly)
    if args.case in ("all", "full-orbit"):
        full_orbit = run_case(source, args.executable.resolve(), launcher,
                              "full-orbit", add_full_orbit)
        check_full_orbit(full_orbit)
        outside = run_external_wall(source, args.executable.resolve(), launcher)
        check_full_orbit(outside, wall_from_vmec=False)
        outside_hit = outside["wall_hit_valid"] == 1
        assert np.any(outside_hit & (outside["wall_hit_field_valid"] == 0))
    print("BEAMS3D event-mode checks: PASS")


if __name__ == "__main__":
    main()
