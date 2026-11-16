"""Attach an analytic Gaussian polymer to a numerical Crambin atom cloud.

The polymer end is connected to the backbone C atom in Crambin's final
residue. The combined structure is evaluated through World's mixed numerical
traversal: the PDB cloud supplies numerical Debye terms, while GaussianPolymer
supplies its existing analytic scattering expressions.
"""

import math
import tempfile

import pyseb

from crambin_form_factor import (
    crambin_profile,
    download_crambin,
    logarithmic_q,
)


def c_terminal_carbon(atoms):
    """Return the backbone C atom belonging to the final parsed residue."""
    final_atom = atoms[-1]
    final_residue = (
        final_atom.chain_id,
        final_atom.residue_number,
        final_atom.insertion_code,
    )
    candidates = [
        atom
        for atom in atoms
        if (
            atom.chain_id,
            atom.residue_number,
            atom.insertion_code,
        )
        == final_residue
        and atom.atom_name == "C"
    ]
    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one backbone C atom in the final residue"
        )
    return candidates[0]


def total_beta(atoms, profile):
    return sum(
        profile.resolve(atom).beta * atom.occupancy
        for atom in atoms
    )


def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = download_crambin(tmpdir)
        atoms = pyseb.PDBParser().parse_file(str(pdb_path))
        profile = crambin_profile()

        terminal_atom = c_terminal_carbon(atoms)
        build_options = pyseb.AtomCloudBuildOptions()
        build_options.reference_atom_serials = {
            "cterminus": terminal_atom.serial,
        }
        cloud = pyseb.AtomCloudBuilder.build(
            atoms,
            profile,
            build_options,
        )

        protein_beta = total_beta(atoms, profile)
        parameters = {
            "beta_polymer": 0.25 * protein_beta,
            "Rg_polymer": 15.0,
        }

        world = pyseb.World()
        graph = world.Add(cloud, "crambin")
        world.Link(
            "GaussianPolymer",
            "polymer.end1",
            "crambin.cterminus",
        )
        world.Add(graph, "crambinpolymer")

        q_values = [0.0] + logarithmic_q(0.005, 1.0, 40)
        protein_form_factor = world.EvaluateFormFactor(
            "crambin",
            parameters,
            q_values,
        )
        conjugate_form_factor = world.EvaluateFormFactor(
            "crambinpolymer",
            parameters,
            q_values,
        )
        protein_rg = math.sqrt(
            world.EvaluateRadiusOfGyration2("crambin", parameters)
        )
        conjugate_rg = math.sqrt(
            world.EvaluateRadiusOfGyration2(
                "crambinpolymer",
                parameters,
            )
        )

        print(f"# scatterers in Crambin: {cloud.scattererCount()}")
        print(
            "# polymer attached to "
            f"{terminal_atom.residue_name}{terminal_atom.residue_number} "
            f"atom {terminal_atom.atom_name} (serial {terminal_atom.serial})"
        )
        print(f"# protein beta: {protein_beta:.10g}")
        print(f"# polymer beta: {parameters['beta_polymer']:.10g}")
        print(f"# polymer Rg [Angstrom]: {parameters['Rg_polymer']:.10g}")
        print(f"# protein Rg [Angstrom]: {protein_rg:.10g}")
        print(f"# conjugate Rg [Angstrom]: {conjugate_rg:.10g}")
        print("# q [1/Angstrom]    F_protein(q)    F_conjugate(q)")
        for q, protein_value, conjugate_value in zip(
            q_values,
            protein_form_factor,
            conjugate_form_factor,
        ):
            if not (
                math.isfinite(protein_value)
                and math.isfinite(conjugate_value)
            ):
                raise RuntimeError(f"Non-finite form factor at q={q}")
            print(
                f"{q:16.8g} "
                f"{protein_value:16.10g} "
                f"{conjugate_value:16.10g}"
            )


if __name__ == "__main__":
    main()
