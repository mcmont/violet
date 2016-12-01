# Third-party module imports
import indigo as indigo_module


def rotatable_bonds(smiles):
    """
    Find the number of rotatable bonds in a molecule.
    Rotatable bonds are defined as any single bond, not in a ring, bound
    to a nonterminal heavy (i.e., non-hydrogen) atom.
    Excluded from the count are amide C-N bonds because of their high
    rotational energy barrier.

    The SMARTS definition is that according to
    J. Med. Chem., 2002, 45 (12), pp 2615-2623
    """
    rotatable = 0
    try:
        # Initialise the Indigo library
        indigo = indigo_module.Indigo()
        # Load the molecule into Indigo
        try:
            mol = indigo.loadMolecule(smiles)
        except IndigoException as e:
            # If loading the molecule fails, try turning off chirality error checking and try again.
            indigo.setOption("ignore-stereochemistry-errors", True)
            mol = indigo.loadMolecule(smiles)
        # Aromatize the molecule
        mol.aromatize()
        matcher = indigo.substructureMatcher(mol)
        q = indigo.loadSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]')
        rotatable = matcher.countMatches(q)

    except IndigoException as e:
        print("Indigo Exception: %s" % (e))
    finally:
        return rotatable

