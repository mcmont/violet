# Third-party module imports
import indigo as indigo_module


def reaction(smarts_string, mol1, mol2=False):
    """
    This is a convenience function that quickly runs a virtual reaction
    and return the products.

    """
    products = []
    try: 
    
        # Setup the Indigo library.
        # Enable generation of multiple products from a reaction, but don't
        # allow product molecules to be fed back in to be reacted a second
        # time around.
        indigo = indigo_module.Indigo()
        indigo.setOption("rpe-multistep-reactions", "true")
        indigo.setOption("rpe-max-depth", "1")
        indigo.setOption("rpe-self-reaction", "false")
        indigo.setOption("rpe-mode", "grid")

        # Load the molecules into Indigo
        molecules = indigo.createArray()
        m1 = indigo.loadMolecule(mol1)
        molecules.arrayAdd(m1)

        if (mol2):
            m2 = indigo.loadMolecule(mol2)
            molecules.arrayAdd(m2)

        # Initialise a 2D array to hold the reactant table
        reactant_table = indigo.createArray()

        # Create the reaction Indigo variable
        rxn = indigo.loadReactionSmarts(smarts_string)
        # Keep mapped atoms
        rxn.automap("keep")

        # Build the Indigo array of reactants
        for i in range(0, rxn.countReactants()):
            array = indigo.createArray()
            if (i < molecules.count()):
                array.arrayAdd(molecules.at(i))
            reactant_table.arrayAdd(array)

        # Enumerate the products
        output_reactions = indigo.reactionProductEnumerate(rxn, reactant_table)
        num_products = output_reactions.count()

        if (num_products > 0):
            # The reaction was successful. Add the products to products[] 
            for i in range(num_products):
                product = output_reactions.at(i)
                for p in product.iterateProducts():
                    products.append(p.canonicalSmiles())

    except IndigoException as e:
        print("Indigo Exception: %s" % (e))

    except Exception as e:
        print("Exception: %s" % (e))

    finally:
        return products

