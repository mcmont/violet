# PLEASE NOTE
# ===========
# 
# The code in this module is a slightly modified version of the code from
# the Chemistry Toolkit Rosetta Wiki.
# http://ctr.wikia.com/wiki/Calculate_TPSA
#
# The algorithm follows the approach of Ertl et al., which is to sum partial
# surface contributions based on fragments defined in a SMARTS string.
# Ertl, Rohde, and Selzer (J. Med. Chem., 43:3714-3717, 2000)
# The SMARTS string is from TJ O'Donnell's CHORD chemistry extension for
# PostgreSQL.

# Core module imports
import collections
import logging
import os

# Third-party module imports
import indigo as indigo_module


def tpsa_count_matches(indigo, subsearch, mol_obj):
    """
    Helper function for tpsa()
    """
    matcher = indigo.substructureMatcher(mol_obj)
    return matcher.countMatches(subsearch)


def tpsa(smiles):
    """
    Compute the topological polar surface area of a molecule specified as a
    SMILES string.
    """

    return_value = False

    # Variables to store the pattern defintions
    Pattern = collections.namedtuple("Pattern", ["value", "subsearch"])
    patterns = []

    try:
        # Initialise the Indigo library
        indigo = indigo_module.Indigo()
        # Build the path to the tpsa data file, relative to this file.
        fn = os.path.join(os.path.dirname(__file__), 'data/tpsa.tab')
        # Get the patterns from the tpsa.tab file, ignoring the header line
        for line in open(fn).readlines()[1:]:
            # Extract the fields
            value, smarts, comment = line.split("\t")
            subsearch = indigo.loadSmarts(smarts)
            # Store for later use
            patterns.append(Pattern(float(value), subsearch))

        # Load the molecule
        mol = indigo.loadMolecule(smiles)
        # Molecules MUST be dearomatized for this TPSA calculation to work correctly.
        mol.dearomatize()
        return_value = sum(tpsa_count_matches(indigo, pattern.subsearch, mol)*pattern.value for pattern in patterns)

    except IndigoException as e:
        logging.error("Indigo exception: %s" % (e))
    except Exception as e:
        logging.error("Exception: %s" % (e))
    finally:
        return return_value


