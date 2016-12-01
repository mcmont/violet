# violet
A library of useful helper functions to the Indigo computational chemistry toolkit, written in Python.

## Installation
Make sure you have a working version of the Indigo Python bindings installed. Please follow the [installation guide on EPAM's website](http://lifescience.opensource.epam.com/indigo/index.html). 

Clone this repo and move it into your PYTHONPATH.

To check that everything's working correctly unit tests are included. To run them:

    python tests.py
    
## Modules
murcko.py - Generate the Bemis-Murcko framework of an Indigo object.
murcko_alpha.py - Generate the Bemis-Murcko framework of an Indigo object with the alpha connection atoms intact.
reaction.py - Run a chemical reaction by specifying a SMARTS string and one or two SMILES strings.
rotatable_bonds.py - Find the number of rotatable bonds in a molecule.
sp3carbon.py - Calculate the number of carbon atoms, sp3 hybridised carbon atoms and sp3 fraction of a SMILES string.
tpsa.py - Compute the topological polar surface area of a molecule specified as a SMILES string. Derived from code from the [Chemistry Toolkit Rosetta Wiki](http://ctr.wikia.com/wiki/Calculate_TPSA).

## License
Violet is made available under the Creative Commons Attribution (CC-BY) license v4.0.
[https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/ "View license")

## Acknowledgements
This code was developed as part of the [LLAMA (Lead-Likeness and Molecular Analysis)](https://llama.leeds.ac.uk) software by Dr. Chris Empson at the School of Chemistry, University of Leeds, UK. 

We thank [EPSRC](https://www.epsrc.ac.uk/) and [GSK](http://www.gsk.com/) for funding, award reference [EP/J00894X/1](http://gow.epsrc.ac.uk/NGBOViewGrant.aspx?GrantRef=EP/J00894X/1).
