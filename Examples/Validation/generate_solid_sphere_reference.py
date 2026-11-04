"""Regenerate the R=1 sphere form-factor fixtures without importing pySEB.

Uses 80-digit mpmath evaluation of the continuous-sphere amplitude squared.
Preserves the existing q grid and both historical filename aliases.
"""
from pathlib import Path
import hashlib
import json
import mpmath as mp

ROOT=Path(__file__).resolve().parent
FOLDER=ROOT/'SolidSphere_R1'
OLD_MISLABELED_SHA256='a9598e7703294f6ebfdeb9d6093bbe327b463aa0934c4c45cccfccc45a1821ac'


def main():
    qstrings=[line.split()[0] for line in (FOLDER/'F.dat').read_text().splitlines() if line.strip()]
    with mp.workdps(80):
        q=[mp.mpf(value) for value in qstrings]
        if len(q)!=201 or any(value<=0 for value in q) or any(a>=b for a,b in zip(q,q[1:])):
            raise ValueError('Unexpected historical q grid')
        values=[(3*(mp.sin(x)-x*mp.cos(x))/x**3)**2 for x in q]
        text=''.join(f'{qvalue}\t{float(value):.17g}\n' for qvalue,value in zip(qstrings,values))
    for name in ('F.dat','FF.dat'):
        (FOLDER/name).write_text(text)
    sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    provenance=dict(formula='[3 (sin(q R) - q R cos(q R)) / (q R)^3]^2',
                    radius=1,precision_decimal_digits=80,mpmath_version=mp.__version__,
                    q_grid='201 unchanged q values from the historical SolidSphere_R1 fixture',
                    generator_sha256=sha(Path(__file__)),
                    replaces_mislabeled_sha256=OLD_MISLABELED_SHA256,
                    reason='The previous F.dat and FF.dat were exact copies of the thick spherical-shell form-factor fixture.',
                    output_sha256={name:sha(FOLDER/name) for name in ('F.dat','FF.dat')})
    (FOLDER/'reference-provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print('Regenerated both solid-sphere aliases on the unchanged 201-point q grid.')


if __name__=='__main__':main()
