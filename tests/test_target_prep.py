"""Scientific regression checks. Run after prep: python3 -I -m unittest discover -s tests -v."""
import math
from pathlib import Path
import runpy
import unittest

M = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/prep_target.py"))


class TargetPreparationTests(unittest.TestCase):
    def test_all_optimal_alignment_mappings_are_exposed(self):
        _, _, choices = M["align"]("AA", "A", with_options=True)
        self.assertEqual(choices[1], {(1, "A"), (None, "-")})
        self.assertEqual(choices[2], {(1, "A"), (None, "-")})

    def test_sasa_of_single_carbon_has_analytic_area(self):
        atom = {("A", 1): {"xyz": [(0, 0, 0)], "names": ["CB"], "elements": ["C"]}}
        area = M["sidechain_sasa"](atom, "A", (1, 1))[1]
        self.assertAlmostEqual(area, 4 * math.pi * (1.7 + 1.4)**2)

    def test_sasa_accounts_for_other_chain_occlusion(self):
        atom = {("A", 1): {"xyz": [(0, 0, 0)], "names": ["CB"], "elements": ["C"]}}
        free = M["sidechain_sasa"](atom, "A", (1, 1))[1]
        atom[("B", 1)] = {"xyz": [(3, 0, 0)], "names": ["CB"], "elements": ["C"]}
        bound = M["sidechain_sasa"](atom, "A", (1, 1))[1]
        self.assertGreater(bound, 0)
        self.assertLess(bound, free)

    def test_3s9l_assembly_uses_symmetry_mate_not_chain_b(self):
        path = M["RAW"] / "3S9L.pdb"
        original = M["read_atoms"](path)
        assembly = M["biological_assembly"](path, 1)
        self.assertEqual({c for c, _ in assembly}, {"A:1", "A:2", "C:1", "C:2"})
        x, y, z = original[("A", 240)]["ca"]
        self.assertEqual(assembly[("A:2", 240)]["ca"], (231.676 - x, y, -z))

    def test_second_3s9l_assembly_has_its_own_transform(self):
        path = M["RAW"] / "3S9L.pdb"
        original = M["read_atoms"](path)
        assembly = M["biological_assembly"](path, 2)
        x, y, z = original[("B", 240)]["ca"]
        self.assertEqual(assembly[("B:2", 240)]["ca"], (y, x, -126.71625 - z))

    def test_old_crystal_contact_is_not_biological_dimer_contact(self):
        path = M["RAW"] / "3S9L.pdb"
        raw = M["read_atoms"](path)
        assembly = M["biological_assembly"](path, 1)
        self.assertIn(211, M["touching"](raw, "A", (190, 382), {"B"}))
        self.assertNotIn(211, M["touching"](assembly, "A:1", (190, 382), {"A:2"}))

    def test_recommended_patch_matches_both_uniprot_sequences(self):
        human, mouse = M["uniprot"]("P02786"), M["uniprot"]("Q62351")
        _, _, choices = M["align"](human, mouse, with_options=True)
        for h, m, aa in ((237, 239, "F"), (240, 242, "K"), (241, 243, "K")):
            self.assertEqual(human[h - 1], aa)
            self.assertEqual(mouse[m - 1], aa)
            self.assertEqual(choices[h], {(m, aa)})


if __name__ == "__main__":
    unittest.main()
