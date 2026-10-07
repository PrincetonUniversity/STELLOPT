#!/usr/bin/env python3
"""Compile actual native axis-storage behavior with no full VMEC dependencies."""
import argparse
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import tempfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_after(text, anchor):
    if text.count(anchor) != 1:
        raise ValueError(f'Expected one production anchor: {anchor!r}')
    return text.index(anchor)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran',
                        help='Compiler command, optionally with flags')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    native = root/'VMEC2000/Sources/General/symforce.f'
    production = root/'VMEC2000/Sources/General/funct3d.f'
    fixture = here/'serial_axis_storage_oracle.f90'
    source = native.read_text()
    routine_start = unique_after(source, '      SUBROUTINE symforce(')
    routine_end = source.index('      END SUBROUTINE symforce', routine_start)
    routine = source[routine_start:routine_end + len('      END SUBROUTINE symforce')] + '\n'
    source = production.read_text().split('      END SUBROUTINE funct3d_par')[1]
    save_start = unique_after(source, '         IF (lmove_axis .AND. iter2.EQ.1) THEN')
    save_end = unique_after(source, '         CALL symforce (')
    restore_start = unique_after(source, '         IF (ALLOCATED(axis_r_save)) THEN')
    restore_end = source.index('         END IF', restore_start) + len('         END IF')
    save = source[save_start:save_end]
    restore = source[restore_start:restore_end] + '\n'
    with tempfile.TemporaryDirectory(prefix='vmec-axis-storage-') as tmp:
        build = Path(tmp)
        (build/'symforce_serial.f').write_text(routine)
        (build/'serial_axis_save.inc').write_text(save)
        include = build/'serial_axis_restore.inc'
        command = shlex.split(args.compiler) + [
            '-O0', '-g', '-fcheck=all', '-ffixed-line-length-none',
            '-I'+str(build), str(fixture), str(build/'symforce_serial.f'),
            '-o', str(build/'oracle')]
        records = []
        for label, block in [('repaired', restore), ('without_restore', ''), ('later_iteration', restore)]:
            include.write_text(block)
            compiled = subprocess.run(command, cwd=build, capture_output=True,
                                      text=True, timeout=60)
            if compiled.returncode != 0:
                raise RuntimeError(compiled.stdout + compiled.stderr)
            oracle_command = [str(build/'oracle')]
            if label == 'later_iteration':
                oracle_command.append('later_iteration')
            result = subprocess.run(oracle_command, capture_output=True,
                                    text=True, timeout=10)
            if label != 'without_restore':
                if result.returncode != 0:
                    raise AssertionError(result.stdout + result.stderr)
            elif (result.returncode == 0 or
                  'Axis scan would receive corrupted geometry' not in result.stderr):
                raise AssertionError('Missing restoration did not fail the geometry oracle')
            records.append(dict(control=label, process_exit=result.returncode,
                                stdout=result.stdout.strip(),
                                expected_geometry_rejection=label == 'without_restore'))
        print(json.dumps(dict(production_sha256=sha(production),
            symforce_source_sha256=sha(native), fixture_sha256=sha(fixture),
            extraction_driver_sha256=sha(Path(__file__)), controls=records,
            scope='Actual native scratch overwrite and source-extracted R/Z preservation; no PDE solve'),
            indent=2))


if __name__ == '__main__':
    main()
