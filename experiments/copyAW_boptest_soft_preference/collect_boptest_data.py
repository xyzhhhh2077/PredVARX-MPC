"""Collect an auditable BOPTEST trajectory for the copyAW MATLAB experiment."""

from __future__ import annotations

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from scipy.io import savemat


API_BASE = "https://api.boptest.net"
TESTCASE = "multizone_office_complex_air"
FLOORS = (1, 2, 3)
ZONES = ("Cor", "Eas", "Nor", "Sou", "Wes")
TEMPERATURE_POINTS = tuple(
    f"floor{floor}_reaZon{zone}_TZon_y" for floor in FLOORS for zone in ZONES
)
CO2_POINTS = tuple(
    f"floor{floor}_reaZon{zone}_CO2Zon_y" for floor in FLOORS for zone in ZONES
)
OUTPUT_POINTS = TEMPERATURE_POINTS + CO2_POINTS
CONTROL_POINTS = tuple(
    point
    for floor in FLOORS
    for point in (f"floor{floor}_TSupAirSet_u", f"floor{floor}_dpSet_u")
)
ACTIVATE_POINTS = tuple(name.replace("_u", "_activate") for name in CONTROL_POINTS)


def build_dataset(records: list[dict]) -> dict[str, np.ndarray]:
    """Convert API records to the fixed 30-output/6-input column contract."""
    if not records:
        raise ValueError("At least one BOPTEST record is required")

    segment_id = np.asarray([record["segment_id"] for record in records], dtype=np.int32)
    step = np.asarray([record["step"] for record in records], dtype=np.int32)
    timestamp = np.asarray([record["time"] for record in records], dtype=float)

    for segment in np.unique(segment_id):
        segment_time = timestamp[segment_id == segment]
        if segment_time.size > 1 and np.any(np.diff(segment_time) <= 0):
            raise ValueError(f"Time must be strictly increasing inside segment {segment}")

    def matrix(group: str, names: tuple[str, ...]) -> np.ndarray:
        columns = []
        for sample_index, record in enumerate(records):
            values = record[group]
            missing = [name for name in names if name not in values]
            if missing:
                raise ValueError(
                    f"Record {sample_index} is missing {group} field {missing[0]}"
                )
            columns.append([float(values[name]) for name in names])
        return np.asarray(columns, dtype=float).T

    return {
        "y": matrix("measurements", OUTPUT_POINTS),
        "u": matrix("controls", CONTROL_POINTS),
        "segment_id": segment_id,
        "step_index": step,
        "time": timestamp,
        "output_names": np.asarray(OUTPUT_POINTS, dtype=object),
        "input_names": np.asarray(CONTROL_POINTS, dtype=object),
    }


def prbs_controls(segment: int, sample: int, rng: np.random.Generator) -> dict[str, float]:
    """Return bounded piecewise-constant excitation, changing every four samples."""
    block = sample // 4
    block_rng = np.random.default_rng(20260801 + 1000 * segment + block)
    supply = block_rng.choice((289.15, 293.15), size=3)
    pressure = block_rng.choice((120.0, 320.0), size=3)
    # Add a small deterministic dither so the six channels are not identical PRBS copies.
    supply += rng.uniform(-0.20, 0.20, size=3)
    pressure += rng.uniform(-5.0, 5.0, size=3)
    values = np.empty(6)
    values[0::2] = supply
    values[1::2] = pressure
    return {name: float(value) for name, value in zip(CONTROL_POINTS, values)}


def request_json(session: requests.Session, method: str, url: str, **kwargs) -> dict:
    last_error = None
    for attempt in range(4):
        try:
            response = session.request(method, url, timeout=180, **kwargs)
            response.raise_for_status()
            payload = response.json()
            if isinstance(payload, dict) and payload.get("status", 200) >= 400:
                raise RuntimeError(payload.get("message", "BOPTEST API error"))
            return payload
        except (requests.RequestException, ValueError, RuntimeError) as exc:
            last_error = exc
            if attempt == 3:
                break
            time.sleep(2**attempt)
    raise RuntimeError(f"BOPTEST request failed: {method} {url}: {last_error}")


def collect_segment(
    session: requests.Session,
    segment: int,
    start_time: int,
    samples: int,
    step_seconds: int,
) -> tuple[list[dict], dict]:
    select = request_json(
        session,
        "POST",
        f"{API_BASE}/testcases/{TESTCASE}/select",
    )
    testid = select["testid"]
    rng = np.random.default_rng(20260801 + segment)
    records = []
    metadata = {"testid": testid, "segment": segment, "start_time": start_time}
    try:
        inputs = request_json(session, "GET", f"{API_BASE}/inputs/{testid}")["payload"]
        measurements = request_json(
            session, "GET", f"{API_BASE}/measurements/{testid}"
        )["payload"]
        missing_inputs = [name for name in CONTROL_POINTS + ACTIVATE_POINTS if name not in inputs]
        missing_outputs = [name for name in OUTPUT_POINTS if name not in measurements]
        if missing_inputs or missing_outputs:
            raise RuntimeError(
                f"BOPTEST contract changed: inputs={missing_inputs}, outputs={missing_outputs}"
            )
        metadata["input_metadata"] = {name: inputs[name] for name in CONTROL_POINTS}
        metadata["output_metadata"] = {name: measurements[name] for name in OUTPUT_POINTS}

        request_json(
            session,
            "PUT",
            f"{API_BASE}/step/{testid}",
            json={"step": step_seconds},
        )
        request_json(
            session,
            "PUT",
            f"{API_BASE}/initialize/{testid}",
            json={"start_time": start_time, "warmup_period": 0},
        )

        for sample in range(samples):
            controls = prbs_controls(segment, sample, rng)
            body = dict(controls)
            body.update({name: 1 for name in ACTIVATE_POINTS})
            response = request_json(
                session,
                "POST",
                f"{API_BASE}/advance/{testid}",
                json=body,
            )["payload"]
            records.append(
                {
                    "segment_id": segment,
                    "step": sample,
                    "time": float(response["time"]),
                    "measurements": {name: response[name] for name in OUTPUT_POINTS},
                    "controls": controls,
                }
            )
            if (sample + 1) % 12 == 0 or sample + 1 == samples:
                print(f"segment {segment}: {sample + 1}/{samples}", flush=True)
    finally:
        try:
            session.put(f"{API_BASE}/stop/{testid}", timeout=180)
        except requests.RequestException:
            pass
    return records, metadata


def write_checkpoint(path: Path, records: list[dict], metadata: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(
        json.dumps({"records": records, "metadata": metadata}, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)


def read_checkpoint(path: Path, segment: int, samples: int) -> tuple[list[dict], dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload["records"]
    if len(records) != samples or any(record["segment_id"] != segment for record in records):
        raise ValueError(f"Checkpoint {path} does not match segment {segment} / {samples} samples")
    build_dataset(records)
    return records, payload["metadata"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-per-segment", type=int, default=24)
    parser.add_argument("--step-seconds", type=int, default=900)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "data")
    parser.add_argument(
        "--segments",
        type=int,
        nargs="+",
        choices=(1, 2, 3),
        default=(1, 2, 3),
        help="Segments to collect; valid existing checkpoints are reused.",
    )
    args = parser.parse_args()
    if args.samples_per_segment < 24:
        parser.error("samples-per-segment must be at least 24")

    # Three seasonal pilot segments; segments 1-2 train and segment 3 validates.
    starts = (14 * 86400, 104 * 86400, 195 * 86400)
    session = requests.Session()
    args.output.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = args.output / "checkpoints"
    all_records: list[dict] = []
    segment_metadata: list[dict] = []
    for segment, start_time in enumerate(starts, start=1):
        checkpoint = checkpoint_dir / f"segment_{segment}.json"
        if checkpoint.exists():
            records, metadata = read_checkpoint(
                checkpoint, segment, args.samples_per_segment
            )
            print(f"segment {segment}: reused checkpoint", flush=True)
        elif segment in args.segments:
            records, metadata = collect_segment(
                session, segment, start_time, args.samples_per_segment, args.step_seconds
            )
            write_checkpoint(checkpoint, records, metadata)
            print(f"segment {segment}: checkpoint saved", flush=True)
        else:
            raise RuntimeError(
                f"Segment {segment} is missing; include it in --segments or provide {checkpoint}"
            )
        all_records.extend(records)
        segment_metadata.append(metadata)

    dataset = build_dataset(all_records)
    provenance = {
        "api_base": API_BASE,
        "testcase": TESTCASE,
        "collected_utc": datetime.now(timezone.utc).isoformat(),
        "step_seconds": args.step_seconds,
        "samples_per_segment": args.samples_per_segment,
        "train_segments": [1, 2],
        "validation_segments": [3],
        "study_scope": "three-season pilot after public instance allocation blocked segment 4",
        "control_policy": "bounded two-level PRBS base held four samples with samplewise small deterministic dither; seed 20260801",
        "segments": segment_metadata,
    }
    savemat(
        args.output / "boptest_multizone_office_complex_air.mat",
        {**dataset, "provenance_json": json.dumps(provenance, sort_keys=True)},
        do_compression=True,
    )
    (args.output / "boptest_multizone_office_complex_air_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(f"saved {dataset['y'].shape[1]} samples to {args.output}")


if __name__ == "__main__":
    main()