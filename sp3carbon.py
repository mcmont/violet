# Core imports
import logging

# Third-party module imports
import indigo as indigo_module


def sp3carbon(smiles):
    """
    Calculate the number of carbon atoms, sp3 hybridised carbon atoms and sp3
    fraction of a SMILES string.
    """

    return_value = False
    num_carbon = float(0)
    num_carbon_sp3 = float(0)
    sp3_fraction = float(0)

    try:
        indigo = indigo_module.Indigo()
        try:
            mol = indigo.loadMolecule(smiles)
        except IndigoException as e:
            # If loading the molecule fails, try turning off chirality error checking and try again.
            indigo.setOption("ignore-stereochemistry-errors", True)
            mol = indigo.loadMolecule(smiles)

        mol.aromatize()

        mol.unfoldHydrogens() # Must add hydrogens to get the correct answer
        for atom in mol.iterateAtoms():
            if (atom.symbol() == 'C'):
                num_carbon = num_carbon + 1
                if (atom.degree() == 4):
                    num_carbon_sp3 = num_carbon_sp3 + 1

        if num_carbon == 0:
            sp3_fraction = 0
        else:
            sp3_fraction = num_carbon_sp3 / num_carbon
        return_value = (int(num_carbon), int(num_carbon_sp3), sp3_fraction)

    except IndigoException as e:
        logging.error("sp3carbon(): Indigo exception: %s" % (e))
    except Exception as e:
        logging.error("sp3carbon(): Exception: %s" % (e))
    finally:
        return return_value


