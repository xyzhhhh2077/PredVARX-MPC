import sys
import unittest
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from prepare_pelican_data import (  # noqa: E402
    ALL_SEGMENTS,
    INPUT_NAMES,
    OUTPUT_NAMES,
    TRAIN_SEGMENTS,
    VALIDATION_SEGMENTS,
    build_dataset,
    flight_to_arrays,
    load_flights,
)

SOURCE_MAT = HERE.parent / "data" / "AscTec_Pelican_Flight_Dataset.mat"
DERIVED_MAT = HERE.parent / "data" / "copyAY_pelican_dataset.mat"


class PelicanDatasetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Derived dataset contract is mandatory (it is what MATLAB consumes).
        if not DERIVED_MAT.exists():
            raise unittest.SkipTest("derived dataset not built; run prepare_pelican_data.py first")
        import scipy.io as sio

        cls.dataset = sio.loadmat(DERIVED_MAT, squeeze_me=True, struct_as_record=False)

    def test_fixed_10_output_4_input_contract(self):
        self.assertEqual(self.dataset["y"].shape, (10, 1_388_410))
        self.assertEqual(self.dataset["u"].shape, (4, 1_388_410))
        self.assertEqual(self.dataset["segment_id"].size, 1_388_410)
        self.assertEqual(self.dataset["output_names"].tolist(), list(OUTPUT_NAMES))
        self.assertEqual(self.dataset["input_names"].tolist(), list(INPUT_NAMES))

    def test_segment_ids_cover_exactly_54_flights_in_order(self):
        ids = np.unique(self.dataset["segment_id"])
        self.assertEqual(ids.tolist(), list(ALL_SEGMENTS))

    def test_segment_boundaries_are_exact_flight_boundaries(self):
        ids = self.dataset["segment_id"]
        # each flight must be contiguous: no sample of flight k inside flight j
        for seg in ALL_SEGMENTS:
            idx = np.flatnonzero(ids == seg)
            self.assertEqual(idx.tolist(), list(range(idx[0], idx[-1] + 1)))
        # first flight starts at 0, flights in increasing order
        starts = []
        for seg in ALL_SEGMENTS:
            idx = np.flatnonzero(ids == seg)
            starts.append(idx[0])
        self.assertEqual(starts, sorted(starts))
        self.assertEqual(starts[0], 0)

    def test_dataset_is_finite(self):
        self.assertTrue(np.all(np.isfinite(self.dataset["y"])))
        self.assertTrue(np.all(np.isfinite(self.dataset["u"])))

    def test_vel_excluded_because_linear_dependence(self):
        # The reason Vel/pqr are excluded: Vel == 100*diff(Pos) up to smoothing,
        # i.e. C_y would be singular. Verify the residual of Vel ~ [dPos, 1].
        if not SOURCE_MAT.exists():
            self.skipTest("source .mat not present")
        flights = load_flights(SOURCE_MAT)
        f = flights[0]
        P = np.asarray(f.Pos, dtype=float)
        V = np.asarray(f.Vel, dtype=float)
        dP = np.diff(P, axis=0)
        n = min(len(V), len(dP))
        A = np.hstack([dP[:n], np.ones((n, 1))])
        coef, *_ = np.linalg.lstsq(A, V[:n], rcond=None)
        pred = A @ coef
        rel = np.linalg.norm(V[:n] - pred) / np.linalg.norm(V[:n])
        self.assertLess(rel, 1e-6)

    def test_train_validation_split_is_strict(self):
        train_ids = np.unique(self.dataset["segment_id"][
            np.isin(self.dataset["segment_id"], list(TRAIN_SEGMENTS))])
        val_ids = np.unique(self.dataset["segment_id"][
            np.isin(self.dataset["segment_id"], list(VALIDATION_SEGMENTS))])
        self.assertEqual(train_ids.tolist(), list(TRAIN_SEGMENTS))
        self.assertEqual(val_ids.tolist(), list(VALIDATION_SEGMENTS))
        self.assertFalse(set(train_ids.tolist()) & set(val_ids.tolist()))

    def test_flight_arrays_shapes(self):
        if not SOURCE_MAT.exists():
            self.skipTest("source .mat not present")
        flights = load_flights(SOURCE_MAT)
        for f in flights[:5]:
            y, u = flight_to_arrays(f)
            n = int(f.len)
            self.assertEqual(y.shape, (10, n))
            self.assertEqual(u.shape, (4, n))


if __name__ == "__main__":
    unittest.main()
