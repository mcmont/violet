# Core imports
import logging

# Third-party module imports
import indigo as indigo_module


def murcko(indigo_obj):
    """
    Generate the Bemis-Murcko scaffold for a molecule.
    """
    indigo = indigo_module.Indigo()
    if not indigo:
        raise Exception('The Indigo library was not initialised.')
    mol = indigo_obj.clone()

    try:
        hanging_atoms = []
        while True:
            # Exit the while loop if we've removed all of the atoms
            if mol == None:
                break
            natoms = mol.countAtoms()
            nei = False
            found = False
            present = False

            # Atoms with 0 or 1 connections are appended to the 
            # hanging_atoms[] list
            for atom in mol.iterateAtoms():
                if (atom.degree() <= 1):
                    hanging_atoms.append(atom.index())

            while (len(hanging_atoms) > 0):
                # Clear down the to_remove and hanging_next lists
                to_remove = []
                hanging_next = []

                for atom in hanging_atoms:
                    if (mol.getAtom(atom).degree() == 0):
                        to_remove.append(atom)
                    else:
                        nei = mol.getAtom(atom).iterateNeighbors().next()

                    if (nei.degree() <= 2 or nei.bond().bondOrder() == 1):
                        to_remove.append(atom)
                    else:
                        present = False
                        for a in hanging_next:
                            if (a == atom):
                                present = True
                        if (present):
                            hanging_next.append(atom)

                for atom in to_remove:
                    if (mol.getAtom(atom).degree() > 0):
                        nei = mol.getAtom(atom).iterateNeighbors().next()
                        if (nei.degree() == 2):
                            found = False
                            for a in to_remove:
                                if (a == nei.index()):
                                    found = True
                                    break

                            if not found:
                                for a in hanging_next:
                                    if (mol.getAtom(a).index() == nei.index()):
                                        found = True
                                        break
                                if not found:
                                    hanging_next.append(nei.index())

                # If there are no atoms left in the to_remove list, break out
                # of the while (len(hanging_atoms) > 0) loop.
                if (len(to_remove) == 0):
                    break

                # Remove the atoms whose indices are in the to_remove array
                mol.removeAtoms(to_remove)
                hanging_atoms = hanging_next

            if mol == None:
                break
            # If we haven't removed any atoms during the last cycle, we've finished.
            if (natoms == mol.countAtoms()):
                break

            # If the number of remaining atoms is less than 3, there can't be a ring, therefore there can't be a scaffold.
            # Return a 'None' result
            if (mol.countAtoms() < 3):
                mol = None

    except (Exception, TypeError, ValueError) as e:
        logging.error("Exception: %s" % (e))
        mol = None

    except IndigoException as e:
        logging.error("IndigoException: %s" % (e))
        mol = None

    finally:
        return mol

