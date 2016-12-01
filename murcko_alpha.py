# Core imports
import logging

# Third-party module imports
import indigo as indigo_module

# Violet library imports
from violet import murcko


def murcko_alpha(indigo_obj, murcko_fwk=False):
    """
    Generate the Bemis-Murcko scaffold for a molecule, but leave the 
    alpha atoms attached.

    The general strategy is:
      1) Obtain the molecule's Murcko framework, which is either passed as an
         argument or calculated if necessary.
      2) Find the atom ID numbers of the Murcko framework ring atoms.
      3) Subtract the ring atom IDs from all atom IDs in the Murcko 
         scaffold to recover the remaining (linker) atoms.
      4) Highlight the linker atoms.
      5) Loop until the atom count hasn't changed between two cycles:
         a) If an atom has one or zero bonds, add it to the hanging_atoms list.
            For all hanging atoms:
            i) If the atom isn't connected to any others, flag it for deletion.
           ii) Check whether the neighbouring atom is part of a ring. If it is, 
               we don't want to delete the current atom because it is a ring alpha
               atom.
          iii) If the neighbouring atom isn't in a ring, and 
               (the neighbouring atom is bonded to 2 or fewer atoms and the bond order is 1) and
               the neighbouring atom isn't highlighted,
               flag the atom for deletion.
         b) Remove all atoms that have been flagged for deletion.
      6) Un-highlight all atoms.
    """

    try:
        indigo = indigo_module.Indigo()

        if not indigo:
            raise Exception('The Indigo library instance was passed uninitialised')

        mol = indigo_obj.clone()

        if (murcko_fwk == False):
            # Get the molecule's Murcko framework.
            murcko_fwk = murcko(mol)

        if (murcko_fwk is None):
            raise Exception("This molecule doesn't have a Murcko framework. Stopping.") 

        # Initialise a set to hold the atom indices of the Murcko framework.
        murcko_atoms = set()
        # Populate the murcko_atoms set with the atom indices.
        for atom in murcko_fwk.iterateAtoms():
            murcko_atoms.add(atom.index())

        # Initialise a set to hold the atom indices of the ring atoms in the Murcko framework.
        murcko_ring_atoms = set()
        for r in murcko_fwk.iterateRings(3, 1000):
            for atom in r.iterateAtoms():
                #print "Atom index %s is in a ring." % (atom.index())
                murcko_ring_atoms.add(atom.index())

        # Subtract the set of Murcko ring atoms from the whole Murcko set to get the linker atom indices.
        linker_atom_indices = list(murcko_atoms - murcko_ring_atoms)
        # Highlight the linker atoms.
        for atom in mol.iterateAtoms():    
            if atom.index() in linker_atom_indices:
                atom.highlight()

        # Hanging atoms are those with one or fewer bonds to other atoms.
        hanging_atoms = []

        # Loop until done
        while True:
            if mol == None:
                break
            # Record the number of atoms in the molecule. This will be 
            # compared to the number of atoms at the end of this cycle.
            # If the numbers are the same, we'll have finished. 
            natoms = mol.countAtoms()

            # Reset some variables
            nei = False
            found = False
            present = False

            # If an atom is connected to one or fewer other atoms, consider it for deletion.
            for atom in mol.iterateAtoms():
                if (atom.degree() <= 1):
                    hanging_atoms.append(atom.index())

            while (len(hanging_atoms) > 0):

                # Clear down the to_remove and hanging_next lists
                to_remove = []
                hanging_next = []

                for atom in hanging_atoms:
                    if (mol.getAtom(atom).degree() == 0):
                        # The hanging atom is not connected to any other atoms. Add it to to_remove array.
                        to_remove.append(atom)
                    else:
                        nei = mol.getAtom(atom).iterateNeighbors().next()

                        # Check whether the neighbouring atom is part of a ring. If it is, we don't want to delete the current atom.
                        nei_not_in_ring = True
                        # Loop over all bonds in the molecule
                        for bond in mol.iterateBonds():
                            # If the bond is somehow connected to the neighbouring atom, check whether it's a chain or a ring bond.
                            if (bond.source().index() == nei.index()) or (bond.destination().index() == nei.index()):
                                # The Indigo.RING constant is 10
                                if bond.topology() == 10:
                                    # The neighbouring atom is in a ring. The atom will not be deleted.
                                    nei_not_in_ring = False

                        # If the neighbouring atom isn't in a ring, 
                        # and (the neighbouring atom is bonded to 2 
                        # or fewer atoms and the bond order is 1) and 
                        # the neighbouring atom isn't highlighted, 
                        # flag the atom for deletion.
                        if nei_not_in_ring and (nei.degree() <= 2 \
                          or nei.bond().bondOrder() == 1) \
                          and (not(nei.isHighlighted())):
                            # If the atom meets the deletion criteria, add it to the list of atoms to be deleted.
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
                            # Is the neighbour's index() already in the to_remove array?
                            for a in to_remove:
                                if (a == nei.index()):
#                                    print "Atom is already in the to_remove array"
                                    found = True
                                    break

                            if not found:
                                for a in hanging_next:
                                    if (mol.getAtom(a).index() == nei.index()):
                                        found = True
                                        break
                                if not found:
                                    hanging_next.append(nei.index())

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

        # Remove any atom highlights to prevent them from appearing 
        # in the canonical SMILES string when the molecule object is 
        # passed back to the calling function.
        if mol is not None:
            for atom in mol.iterateAtoms():
                if atom.index() in linker_atom_indices:
                    atom.unhighlight()


    except (Exception, TypeError, ValueError) as e:
        logging.error("murcko_alpha(): Exception: %s" % (e))
        mol = None

    except IndigoException as e:
        logging.error("murcko_alpha(): IndigoException: %s" % (e))
        mol = None

    finally:
        return mol

