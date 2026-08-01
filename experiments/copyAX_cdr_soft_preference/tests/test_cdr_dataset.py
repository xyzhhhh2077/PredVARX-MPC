import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from collect_cdr_data import (  # noqa: E402
    INPUT_NAMES,
    OUTPUT_NAMES,
    PRBS_AMPLITUDE,
    PRBS_HOLD,
    SEGMENT_BASE_SEED,
    STEPS_PER_SEGMENT,
    TRAIN_SEGMENTS,
    VALIDATION_SEGMENTS,
    build_dataset,
    prbs_sequence,
)


class CdrDatasetTest(unittest.TestCase):
    def test_fixed_30_output_8_input_contract(self):
        dataset = build_dataset()
        total = len(TRAIN_SEGMENTS + VALIDATION_SEGMENTS) * STEPS_PER_SEGMENT
        self.assertEqual(dataset["y"].shape, (30, total))
        self.assertEqual(dataset["u"].shape, (8, total))
        self.assertEqual(
            dataset["segment_id"].tolist(),
            [s for s in TRAIN_SEGMENTS + VALIDATION_SEGMENTS for _ in range(STEPS_PER_SEGMENT)],
        )
        self.assertEqual(dataset["output_names"].tolist(), list(OUTPUT_NAMES))
        self.assertEqual(dataset["input_names"].tolist(), list(INPUT_NAMES))

    def test_prbs_is_deterministic_and_bounded(self):
        u1 = prbs_sequence(1, 100)
        u2 = prbs_sequence(1, 100)
        np.testing.assert_array_equal(u1, u2)
        self.assertTrue(np.all(np.abs(u1) == PRBS_AMPLITUDE))
        # Hold pattern: within each block of PRBS_HOLD time samples the level
        # is constant (u is time-major here: rows = time, columns = channels).
        self.assertTrue(np.all(u1[1::PRBS_HOLD] == u1[0::PRBS_HOLD]))
        self.assertTrue(np.all(u1[2::PRBS_HOLD] == u1[0::PRBS_HOLD]))
        self.assertTrue(np.all(u1[3::PRBS_HOLD] == u1[0::PRBS_HOLD]))

    def test_prbs_differs_across_segments(self):
        self.assertFalse(
            np.array_equal(prbs_sequence(1, 200), prbs_sequence(2, 200))
        )

    def test_dataset_is_finite(self):
        dataset = build_dataset()
        self.assertTrue(np.all(np.isfinite(dataset["y"])))
        self.assertTrue(np.all(np.isfinite(dataset["u"])))

    def test_segment_boundaries_are_exact(self):
        dataset = build_dataset()
        ids = dataset["segment_id"]
        boundaries = np.flatnonzero(np.diff(ids) != 0) + 1
        self.assertEqual(
            boundaries.tolist(),
            [k * STEPS_PER_SEGMENT for k in range(1, len(TRAIN_SEGMENTS + VALIDATION_SEGMENTS))],
        )

    def test_seed_formula_is_stable(self):
        self.assertEqual(SEGMENT_BASE_SEED + 1 * 101, 20260801 + 101)


if __name__ == "__main__":
    unittest.main()
