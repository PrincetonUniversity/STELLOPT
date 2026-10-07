#!/usr/bin/env python3
"""Actual writer/loader component round-trip with independent physical lambda."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def run(source, require_fixed, baseline):
    tests=source/'tests/test_lambda_restart'
    writer=source/'VMEC2000/Sources/Input_Output/wrout.f'
    loader=source/'VMEC2000/Sources/Initialization_Cleanup/load_xc_from_wout.f'
    text=writer.read_text()
    export=re.search(r'^\s*lmns\(:,js\) = \(lmns1\(:\)/phipf\(js\)\) \* lamscale$',text,re.M)
    start=text.index('      WHERE (NINT(xm) .le. 1) lmns(:,1) = lmns(:,2)')
    end=text.index('      lmns(:,1) = 0',start)+len('      lmns(:,1) = 0')
    if not export:
        raise ValueError('Production export convention changed; review the oracle')
    controls=[]
    invalid_controls=[]
    with tempfile.TemporaryDirectory(prefix='vmec-lambda-native-roundtrip-') as tmp:
        tmp=Path(tmp)
        (tmp/'native_export.inc').write_text(export.group().strip()+'\n')
        (tmp/'native_halfmesh_export.inc').write_text(text[start:end]+'\n')
        for name in ['stubs.f90','oracle.f90']:
            shutil.copy2(tests/name,tmp/name)
        binary=tmp/'oracle'
        if baseline:
            raw=subprocess.check_output(['git','-C',str(source),'show',
                '8060f5e5b1bfe11b2f8809c8dfa90459ee72f9ba:VMEC2000/Sources/Initialization_Cleanup/load_xc_from_wout.f'])
            loader=tmp/'baseline_loader.f';loader.write_bytes(raw)
        loader_digest=hashlib.sha256(loader.read_bytes()).hexdigest()
        compiled=subprocess.run(['gfortran','-cpp','-O0','-g','-fcheck=all',
            '-ffpe-trap=invalid,zero,overflow','-ffixed-line-length-72','-I',str(tmp),str(tmp/'stubs.f90'),
            str(loader),str(tmp/'oracle.f90'),'-o',str(binary)],cwd=tmp,
            capture_output=True,text=True,timeout=30)
        if compiled.returncode:
            raise RuntimeError(compiled.stderr)
        controls_to_run=[(scale, sign*scale) for scale in [1.,2.,.5] for sign in [1.,-1.]]
        for scale, phi in controls_to_run:
            result=subprocess.run([str(binary),str(scale),str(phi)],capture_output=True,
                                  text=True,timeout=5)
            values={k:float(v) for k,v in (line.split('=') for line in result.stdout.splitlines())}
            expected_exit=0 if require_fixed or scale==1. else 6
            if result.returncode!=expected_exit:
                raise AssertionError((scale,result.returncode,result.stdout,result.stderr))
            controls.append(dict(lamscale=scale,phipf=phi,process_exit=result.returncode,
                                 independent_errors=values,expected_exit=expected_exit))
        if require_fixed:
            for scale in ['0', '-1', 'NaN', 'Infinity', '-Infinity']:
                result=subprocess.run([str(binary),scale,'1'],capture_output=True,
                                      text=True,timeout=5)
                if result.returncode!=1 or 'Invalid lambda scale in load_xc' not in result.stdout:
                    raise AssertionError((scale,result.returncode,result.stdout,result.stderr))
                invalid_controls.append(dict(scale_argument=scale,process_exit=result.returncode,
                    diagnostic=result.stdout.strip(),scope='Rejected before native reader and division'))
    return dict(source_base=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD']).decode().strip(),
        writer_sha256=hashlib.sha256(writer.read_bytes()).hexdigest(),
        actual_loader_sha256=loader_digest,
        tests_sha256={name:hashlib.sha256((tests/name).read_bytes()).hexdigest()
                      for name in ['stubs.f90','oracle.f90','run.py']},
        mode='fixed allscales mustpass' if require_fixed else 'baseline nonunit scales mustfail independentlambda oracle',
        controls=controls,invalid_scale_controls=invalid_controls,scope='Production wrout normalization/halfmesh block and actual native loader compiled; synthetic read module supplies exactm1/m2 full coefficients. Independent physicallambda/fielddensity, geometry and profileinvariance gates; no PDE or privateinputs.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[2])
    mode=p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--require-fixed',action='store_true')
    mode.add_argument('--baseline',action='store_true')
    p.add_argument('--output',type=Path)
    a=p.parse_args();result=run(a.source,a.require_fixed,a.baseline)
    if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
