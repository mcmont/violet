# Core module imports
import unittest

# Third-party module imports
import indigo as indigo_module

# Violet library imports
from violet import murcko
from violet import murcko_alpha
from violet import reaction
from violet import rotatable_bonds
from violet import sp3carbon
from violet import tpsa


class violetTests(unittest.TestCase):

    def test_murcko(self):
        indigo = indigo_module.Indigo()
        icb = "N#CC1C=CC=C(C=1)N=C=O" # 3-isocyanatobenzonitrile
        mol = indigo.loadMolecule(icb)
        self.assertEqual(murcko(mol).canonicalSmiles(),'C1C=CC=CC=1')

    def test_murcko_alpha(self):
        indigo = indigo_module.Indigo()
        icb = "N#CC1C=CC=C(C=1)N=C=O" # 3-isocyanatobenzonitrile
        mol = indigo.loadMolecule(icb)
        self.assertEqual(murcko_alpha(mol).canonicalSmiles(),'CC1=CC(N)=CC=C1')

    # Bimolecular reaction: Amide formation
    def test_reaction_bimolecular(self):
        amide_formation_smarts = "([NX3;H1,H2;!$(NC=O);!$(NS(=O)=O)]).([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-]),Cl])>>[NX3;H1;!$(NC=O);!$(NS(=O)=O)][CX3;$([R0][#6]),$([H1R0])](=[OX1])"
        mol1 = "C(C1NCCC1C)"
        mol2 = "O=C(Cl)C1CCCCN1"
        self.assertEqual(reaction(amide_formation_smarts, mol1, mol2), ['CC1CCN(C1C)C(=O)C1CCCCN1'])

    # Bimolecular reaction: Amide formation with an invalid reactant
    def test_reaction_bimolecular_expectedfailure(self):
        amide_formation_smarts = "([NX3;H1,H2;!$(NC=O);!$(NS(=O)=O)]).([CX3;$([R0][#6]),$([H1R0])](=[OX1])[$([OX2H]),$([OX1-]),Cl])>>[NX3;H1;!$(NC=O);!$(NS(=O)=O)][CX3;$([R0][#6]),$([H1R0])](=[OX1])"
        mol1 = "C(C1NCCC1C)"
        mol2 = "O=C(C)C1CCCCN1"
        self.assertEqual(reaction(amide_formation_smarts, mol1, mol2), [])

    # Unimolecular reaction: BOC deprotection
    def test_reaction_unimolecular(self):
        boc_deprotection_smarts = "[*;$([NX3]([C,N,O])([CX4])C(=O)OC(C)(C)C):1]C(=O)OC(C)(C)C>>[*:1]"
        mol = "C1CCCCC1N(C)C(=O)OC(C)(C)C"
        self.assertEqual(reaction(boc_deprotection_smarts, mol), ['CNC1CCCCC1'])

    # Unimolecular reaction: BOC deprotection where no BOC group is present
    def test_reaction_unimolecular_expectedfailure(self):
        boc_deprotection_smarts = "[*;$([NX3]([C,N,O])([CX4])C(=O)OC(C)(C)C):1]C(=O)OC(C)(C)C>>[*:1]"
        mol = "C1CCCCC1N"
        self.assertEqual(reaction(boc_deprotection_smarts, mol), [])

    # Calculate the number of rotatable bonds for benzene. The answer should be zero.
    def test_rotatable_bonds_0(self):
        self.assertEqual(rotatable_bonds("c1ccccc1"), 0)
        
    # Calculate the number of rotatable bonds for 3-methylpentane. There are two.
    def test_rotatable_bonds_2(self):
        self.assertEqual(rotatable_bonds("C(C)C(C)CC"), 2)

    # There are no carbon atoms in sulfamoyl chloride, so the answer should contain zeros.
    def test_sp3carbon_nocarbons(self):
        self.assertEqual(sp3carbon("NS(=O)(=O)Cl"), (0, 0, 0))

    # There are no sp3 hybridised carbon atoms in benzene.
    def test_sp3carbon_0(self):
        self.assertEqual(sp3carbon("c1ccccc1"), (6, 0, 0.0))

    # A naked carbon atom is sp3 hybridised.
    def test_sp3carbon_1(self):
        self.assertEqual(sp3carbon("C"), (1, 1, 1.0))

    # benzyl(methyl)amine contains 8 carbon atoms, 2 of which are sp3 hybridised.
    def test_sp3carbon_2(self):
        self.assertEqual(sp3carbon("CNCc1ccccc1"), (8, 2, 0.25))

    # Benzene has a polar surface area of zero.
    def test_tpsa_0(self):
        self.assertEqual(tpsa("c1ccccc1"), 0)

    # 1-hydroxycyclopropane-1-carboxylic acid has a polar surface area of
    # 57.53 angstroms squared.
    def test_tpsa_57(self):
        self.assertEqual(round(tpsa("OC(=O)C1(O)CC1"), 2), 57.53)


if __name__ == '__main__':
    unittest.main()
