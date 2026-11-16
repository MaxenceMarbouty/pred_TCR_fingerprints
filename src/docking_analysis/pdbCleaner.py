from Bio.PDB import PDBParser, PDBIO, Select
import sys
import os

class CleanSelect(Select):
    # Residues to always exclude even if listed as ATOM
    EXCLUDE_RESNAMES = {
        "HOH", "WAT", "TIP", "TIP3",  # waters
        "MG", "ZN", "CA", "NA", "CL", "K",   # ions
    }

    def accept_residue(self, residue):
        # Skip HETATM and insertion codes
        if residue.id[0].strip() and residue.id[0] != ' ':
            return False
        if residue.id[2] != ' ':
            return False
        # Skip residues missing backbone atoms
        atom_names = {a.get_name() for a in residue.get_atoms()}
        if not {' N  ', ' CA ', ' C  '}.issubset(atom_names):
            # Try without spaces in case of naming variation
            atom_names_stripped = {a.get_name().strip() for a in residue.get_atoms()}
            if not {'N', 'CA', 'C'}.issubset(atom_names_stripped):
                print(f"  Dropping residue {residue.id} in chain "
                      f"{residue.get_parent().id} — missing backbone atoms")
                return False
        return True

    def accept_atom(self, atom):
        # Keep only blank altloc or the first one (A)
        return atom.get_altloc() in (' ', 'A')


def renumber_and_clean(input_path, output_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('s', input_path)

    for model in structure:
        for chain in model:
            residues = [r for r in chain.get_residues()
                        if (r.id[0].strip() == '' or r.id[0] == ' ')
                        and r.id[2] == ' ']
            # Filter to only backbone-complete residues before renumbering
            complete = []
            for r in residues:
                names = {a.get_name().strip() for a in r.get_atoms()}
                if {'N', 'CA', 'C'}.issubset(names):
                    complete.append(r)
                else:
                    print(f"  Skipping {r.id} chain "
                          f"{r.get_parent().id} — incomplete backbone")
            for i, residue in enumerate(complete, start=1):
                residue.id = (' ', i, ' ')

        io = PDBIO()
        io.set_structure(structure)
        io.save(output_path, CleanSelect())


if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    os.makedirs(output_path, exist_ok=True)
    for m in os.listdir(input_path):
        if m.endswith(".pdb") or m.endswith(".cif"):
            renumber_and_clean(os.path.join(input_path, m), os.path.join(output_path, m))

