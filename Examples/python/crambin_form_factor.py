"""Download crambin (1CRN) and calculate an illustrative Debye form factor.

The coordinates and radii are in Angstrom, so q is reported in inverse
Angstrom. The beta values below are approximate electron counts and are meant
to demonstrate the PDB-to-cloud pipeline. They are not solvent-corrected SAXS
contrasts or q-dependent atomic form factors.
"""

import hashlib
import math
from pathlib import Path
import tempfile
import urllib.request

import pyseb


PDB_URL = "https://files.rcsb.org/download/1CRN.pdb"
PDB_SHA256 = (
    "42199a30a0701864a2a5cc76cd7f35cc"
    "544cd0e65fbcf63e03c166543249b811"
)


def logarithmic_q(start, stop, count):
    if count < 2:
        return [start]
    ratio = (stop / start) ** (1.0 / (count - 1))
    return [start * ratio**index for index in range(count)]


def crambin_profile():
    profile = pyseb.AtomParameterProfile()
    # van der Waals radius [Angstrom], approximate electron count.
    profile.set_element("C", radius=1.70, beta=6.0)
    profile.set_element("N", radius=1.55, beta=7.0)
    profile.set_element("O", radius=1.52, beta=8.0)
    profile.set_element("S", radius=1.80, beta=16.0)
    return profile


def download_crambin(directory):
    pdb_path = Path(directory) / "1CRN.pdb"
    with urllib.request.urlopen(PDB_URL, timeout=30) as response:
        contents = response.read()

    digest = hashlib.sha256(contents).hexdigest()
    if digest != PDB_SHA256:
        raise RuntimeError(f"Unexpected SHA-256 for 1CRN.pdb: {digest}")
    pdb_path.write_bytes(contents)
    return pdb_path


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = download_crambin(tmpdir)

        cloud = pyseb.StructureCloudLoader.load_pdb(
            str(pdb_path),
            crambin_profile(),
        )
        world = pyseb.World()
        world.Add(cloud, "crambin")

        q_values = [0.0] + logarithmic_q(0.005, 1.0, 40)
        form_factor = world.EvaluateFormFactor(
            "crambin",
            {},
            q_values,
        )
        radius_of_gyration2 = world.EvaluateRadiusOfGyration2(
            "crambin",
            {},
        )
        radius_of_gyration = math.sqrt(radius_of_gyration2)

        print(f"# scatterers: {cloud.scattererCount()}")
        print(f"# Rg [Angstrom]: {radius_of_gyration:.10g}")
        print(f"# Rg^2 [Angstrom^2]: {radius_of_gyration2:.10g}")
        print("# q [1/Angstrom]    F(q)")
        for q, value in zip(q_values, form_factor):
            if not math.isfinite(value):
                raise RuntimeError(f"Non-finite F(q) at q={q}")
            print(f"{q:16.8g} {value:16.10g}")


if __name__ == "__main__":
    main()
