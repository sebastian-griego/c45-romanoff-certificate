from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import verify_c45_romanoff_upper_bound as verifier


class VerifierTests(unittest.TestCase):
    def test_cluster_polynomial_matches_bruteforce_small_cluster(self) -> None:
        self.assertEqual(
            verifier.cluster_polynomial(11, [23, 89]),
            verifier.brute_force_cluster_polynomial(11, [23, 89]),
        )

    def test_multiplier_distribution_has_expected_total_weight(self) -> None:
        dist = verifier.multiplier_distribution()
        expected = 1
        for primes in verifier.ADDED_CLUSTERS.values():
            for prime in primes:
                expected *= prime

        self.assertEqual(sum(dist.values()), expected)

    def test_log2_lower_bound_is_monotone(self) -> None:
        self.assertLess(verifier.log2_lower_bound(10), verifier.log2_lower_bound(20))
        self.assertLess(verifier.log2_lower_bound(20), verifier.log2_lower_bound(40))

    def test_load_seed_histogram_rejects_duplicate_nu(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "hist.txt"
            path.write_text("0 1\n0 2\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "duplicate nu"):
                verifier.load_seed_histogram(path)


if __name__ == "__main__":
    unittest.main()
