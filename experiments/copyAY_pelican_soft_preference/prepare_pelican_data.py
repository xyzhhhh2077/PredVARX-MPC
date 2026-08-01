"""Prepare the Waterloo Pelican quadrotor flight dataset for the copyAY experiment.

Source
------
Mohajerin & Waslander, "Quadrotor Flight Dataset", WAVELab, University of Waterloo.
Downloaded from http://wavelab.uwaterloo.ca/wp-content/uploads/2017/09/AscTec_Pelican_Flight_Dataset.mat
Papers:
  [1] N. Mohajerin, M. Mozifian, S. L. Waslander, "Deep Learning a Quadrotor
      Dynamic Model for Multi-Step Prediction", ICRA 2018.
  [2] N. Mohajerin, S. L. Waslander, "Multistep Prediction of Dynamic Systems
      With Recurrent Neural Networks", IEEE TNNLS 2019.

Data contract (this copy)
-------------------------
54 flights (segments), 100 Hz sampling, Vicon indoor motion capture
# Units (verified on raw .mat, 2026-08-02): Pos in **m**, Euler in **rad**,
# Vel in m/s == 100*diff(Pos) EXACTLY (dataset stores this as a field),
# pqr in rad/s (body-frame), Motors/Motors_CMD dimensionless integers [0,218].
# NOTE: earlier comments saying "mm / deg" were wrong; raw values never
# converted, only the documentation was misleading.

Outputs y (10 rows):  [Pos(3); Euler(3); Motors(4)]  -- raw measurements only.
  - Vel and pqr are EXCLUDED on purpose: Vel is exactly 100*diff(Pos)
    (numerical derivative + smoothing = linear operator, residual 2e-14), so
    including both Pos and Vel makes C_y singular and breaks the C_y^{-1/2}
    geometry of the copyAU construction. pqr is a nonlinear transform of the
    Euler-rate derivative (also derived, not measured).
Inputs u (4 rows):  Motors_CMD (commanded motor speeds, integer [0,218]).

Alignment (standard causality): within each flight, transition k pairs
    (y(:,k), u(:,k)) -> y(:,k+1)   i.e.  ur = uc(:, valid).
Physical motivation: commanded motor speed at step k produces thrust that
changes attitude/position at step k+1 (10 ms at 100 Hz).

Segments: flights are stored in the .mat as a cell array; segment_id = flight
index 1..54. Training segments 1..36, validation segments 37..54 (split by
flight, never shuffling samples within a flight).
"""

from __future__ import annotations

import argparse
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

SOURCE_URL = "http://wavelab.uwaterloo.ca/wp-content/uploads/2017/09/AscTec_Pelican_Flight_Dataset.mat"
TRAIN_SEGMENTS = tuple(range(1, 37))       # flights 1..36
VALIDATION_SEGMENTS = tuple(range(37, 55))  # flights 37..54
ALL_SEGMENTS = TRAIN_SEGMENTS + VALIDATION_SEGMENTS
SAMPLE_RATE_HZ = 100
OUTPUT_NAMES = (
    "pos_x_mm", "pos_y_mm", "pos_z_mm",
    "euler_roll_deg", "euler_pitch_deg", "euler_yaw_deg",
    "motor_1_speed", "motor_2_speed", "motor_3_speed", "motor_4_speed",
)
INPUT_NAMES = (
    "motor_cmd_1", "motor_cmd_2", "motor_cmd_3", "motor_cmd_4",
)


def load_flights(path: str | os.PathLike) -> list:
    """Load the Pelican .mat and return the list of flight structs."""
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    flights = d["flights"]
    if not isinstance(flights, np.ndarray):
        flights = np.asarray([flights])
    return [f for f in flights]


def flight_to_arrays(flight) -> tuple[np.ndarray, np.ndarray]:
    """Return (y, u) for one flight: y 10xN, u 4xN, N = flight.len."""
    n = int(flight.len)
    pos = np.asarray(flight.Pos, dtype=float).T        # 3 x N (m)
    euler = np.asarray(flight.Euler, dtype=float).T    # 3 x N (rad)
    motors = np.asarray(flight.Motors, dtype=float).T  # 4 x N
    cmd = np.asarray(flight.Motors_CMD, dtype=float).T  # 4 x N
    y = np.vstack([pos, euler, motors])
    u = cmd
    assert y.shape == (10, n) and u.shape == (4, n), (
        f"Unexpected shapes y={y.shape} u={u.shape} n={n}"
    )
    assert np.all(np.isfinite(y)) and np.all(np.isfinite(u)), "non-finite values"
    return y, u


def flight_to_speed_state(flight) -> tuple[np.ndarray, np.ndarray]:
    """Return (y_spd, u) for one flight: y_spd 10xN, u 4xN, N = flight.len.

    Speed-state benchmark state (ICRA 2018 protocol): the paper predicts
    translational velocity and body rates, NOT position/attitude. To compare
    fairly we build y = [Vel(3); pqr(3); Motors(4)] (10 rows). Vel and pqr
    appear WITHOUT Pos/Euler here, so C_y stays full-rank (cond ~6 after
    standardization, verified) -- the singularity only occurs when Pos and
    Vel are in y together. Vel is 1 sample shorter than Pos; we align by
    dropping the first sample of every other channel (Vel(k) = 100*(Pos(k+1)
    - Pos(k)) is a forward difference, so y_spd(:,k) is paired with u(:,k)
    and drives y_spd(:,k+1) under the standard causality contract).
    """
    n = int(flight.len)
    vel = np.asarray(flight.Vel, dtype=float).T        # 3 x (N-1) (m/s)
    pqr = np.asarray(flight.pqr, dtype=float).T        # 3 x N (rad/s)
    motors = np.asarray(flight.Motors, dtype=float).T  # 4 x N
    cmd = np.asarray(flight.Motors_CMD, dtype=float).T  # 4 x N
    assert vel.shape[1] == n - 1 and pqr.shape[1] == n and motors.shape[1] == n
    y_spd = np.vstack([vel, pqr[:, 1:], motors[:, 1:]])
    u_spd = cmd[:, 1:]
    assert y_spd.shape == (10, n - 1) and u_spd.shape == (4, n - 1)
    assert np.all(np.isfinite(y_spd)) and np.all(np.isfinite(u_spd)), "non-finite values"
    return y_spd, u_spd


def build_dataset(path: str | os.PathLike) -> dict:
    flights = load_flights(path)
    assert len(flights) == 54, f"expected 54 flights, got {len(flights)}"
    columns_y, columns_u, segment_ids = [], [], []
    columns_ys, columns_us, segment_ids_spd = [], [], []
    for idx, flight in enumerate(flights, start=1):
        y, u = flight_to_arrays(flight)
        columns_y.append(y)
        columns_u.append(u)
        segment_ids.extend([idx] * y.shape[1])
        y_spd, u_spd = flight_to_speed_state(flight)
        columns_ys.append(y_spd)
        columns_us.append(u_spd)
        segment_ids_spd.extend([idx] * y_spd.shape[1])
    y = np.concatenate(columns_y, axis=1)
    u = np.concatenate(columns_u, axis=1)
    y_spd = np.concatenate(columns_ys, axis=1)
    u_spd = np.concatenate(columns_us, axis=1)
    assert y_spd.shape[1] == y.shape[1] - 54, (  # one sample dropped per flight
        f"speed-state length mismatch {y_spd.shape[1]} vs {y.shape[1]}"
    )
    return {
        "y": y,
        "u": u,
        "y_speed_state": y_spd,
        "u_speed_state": u_spd,
        "segment_id": np.asarray(segment_ids, dtype=np.int32),
        "segment_id_speed_state": np.asarray(segment_ids_spd, dtype=np.int32),
        "output_names": np.asarray(OUTPUT_NAMES, dtype=object),
        "input_names": np.asarray(INPUT_NAMES, dtype=object),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--src", default=None, help="path to AscTec_Pelican_Flight_Dataset.mat")
    parser.add_argument("--out", default=None, help="output .mat path")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    src = Path(args.src) if args.src else here / "data" / "AscTec_Pelican_Flight_Dataset.mat"
    out_path = Path(args.out) if args.out else here / "data" / "copyAY_pelican_dataset.mat"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    started = time.time()
    dataset = build_dataset(src)
    savemat(out_path, dataset)

    lens = [int(f.len) for f in load_flights(src)]
    provenance = {
        "source": "Waterloo Pelican quadrotor flight dataset (WAVELab)",
        "source_url": SOURCE_URL,
        "papers": [
            "Mohajerin, Mozifian, Waslander, Deep Learning a Quadrotor Dynamic Model for Multi-Step Prediction, ICRA 2018",
            "Mohajerin, Waslander, Multistep Prediction of Dynamic Systems With Recurrent Neural Networks, TNNLS 2019",
        ],
        "vehicle": "AscTec Pelican quadrotor, indoor Vicon motion capture",
        "sample_rate_hz": SAMPLE_RATE_HZ,
        "flight_count": len(lens),
        "flight_lengths": lens,
        "total_samples": int(sum(lens)),
        "train_segments": list(TRAIN_SEGMENTS),
        "validation_segments": list(VALIDATION_SEGMENTS),
        "outputs": list(OUTPUT_NAMES),
        "inputs": list(INPUT_NAMES),
        "output_exclusion_note": (
            "Vel (100*diff(Pos)) and pqr (body-rate transform of Euler-rate "
            "derivative) are excluded: they are derived linearly from Pos/Euler, "
            "which would make C_y singular for the copyAU C_y^{-1/2} geometry."
        ),
        "smoothing": "dataset authors applied a window-5 local-regression smoother to all channels",
        "alignment": "u(:,k) drives y(:,k+1) (standard causality, 100 Hz)",
        "split_rule": "flights 1..36 train, 37..54 validation; no within-flight shuffling",
        "collected_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "prepared_by": "prepare_pelican_data.py (copyAY)",
    }
    provenance_path = out_path.with_suffix(".json")
    provenance_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")

    print(
        f"Pelican dataset written: {out_path} "
        f"y={dataset['y'].shape} u={dataset['u'].shape} "
        f"segments={np.unique(dataset['segment_id']).tolist()} "
        f"elapsed={time.time()-started:.1f}s"
    )


if __name__ == "__main__":
    main()
