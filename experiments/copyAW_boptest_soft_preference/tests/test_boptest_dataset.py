import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from collect_boptest_data import (  # noqa: E402
    CONTROL_POINTS,
    OUTPUT_POINTS,
    build_dataset,
    read_checkpoint,
    write_checkpoint,
)


class BoptestDatasetTest(unittest.TestCase):
    def make_records(self):
        records = []
        for segment in (1, 2):
            for step in range(4):
                measurements = {
                    name: 1000 * segment + 10 * i + step
                    for i, name in enumerate(OUTPUT_POINTS)
                }
                controls = {
                    name: 100 * segment + 10 * i + step
                    for i, name in enumerate(CONTROL_POINTS)
                }
                records.append(
                    {
                        "segment_id": segment,
                        "step": step,
                        "time": 900 * step,
                        "measurements": measurements,
                        "controls": controls,
                    }
                )
        return records

    def test_builds_fixed_30_output_6_input_contract(self):
        dataset = build_dataset(self.make_records())

        self.assertEqual(dataset["y"].shape, (30, 8))
        self.assertEqual(dataset["u"].shape, (6, 8))
        self.assertEqual(dataset["segment_id"].tolist(), [1] * 4 + [2] * 4)
        self.assertEqual(dataset["output_names"].tolist(), list(OUTPUT_POINTS))
        self.assertEqual(dataset["input_names"].tolist(), list(CONTROL_POINTS))
        self.assertEqual(dataset["y"][0, 0], 1000)
        self.assertEqual(dataset["u"][5, -1], 253)

    def test_rejects_missing_measurement(self):
        records = self.make_records()
        del records[0]["measurements"][OUTPUT_POINTS[-1]]

        with self.assertRaisesRegex(ValueError, OUTPUT_POINTS[-1]):
            build_dataset(records)

    def test_rejects_nonincreasing_time_inside_segment(self):
        records = self.make_records()
        records[2]["time"] = records[1]["time"]

        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            build_dataset(records)

    def test_checkpoint_round_trip_is_validated(self):
        records = self.make_records()[4:7]
        metadata = {"segment": 2, "testid": "transient-test-id"}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "segment_2.json"
            write_checkpoint(path, records, metadata)
            loaded_records, loaded_metadata = read_checkpoint(path, 2, 3)
        self.assertEqual(loaded_records, records)
        self.assertEqual(loaded_metadata, metadata)

    def test_checkpoint_rejects_wrong_sample_count(self):
        records = self.make_records()[4:5]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "segment_2.json"
            write_checkpoint(path, records, {"segment": 2})
            with self.assertRaisesRegex(ValueError, "does not match"):
                read_checkpoint(path, 2, 3)


if __name__ == "__main__":
    unittest.main()