"""Collect an auditable ControlGym CDR dataset for the copyAX MATLAB experiment.

Contract
--------
The CDR environment (controlgym, convection_diffusion_reaction) evolves
    x_{k+1} = A x_k + B2 u_k + w_k,
and returns, at step k, the observation computed from the state BEFORE the
update (controlgym pde.py step order):
    y_k = C x_k + v_k.
Therefore the action u_k sent at step k drives y_{k+1}, returned by step k+1.
This collector stores the aligned arrays
    y(:,k) -- observation y_k (reset obs is y_0; obs returned by step k is y_k)
    u(:,k) -- action u_k sent at step k (u_{T-1} is a sentinel, not used by the
              regression: within a segment, transition k pairs (y_k, u_k) -> y_{k+1})
so that the MATLAB pipeline uses the STANDARD causal alignment
    z_{k+1} = A z_k + B u_k,   i.e.  ur = uc(:, valid).

Reproducibility
---------------
Everything is seeded per segment: the environment seed, the PRBS generator seed
and the initial condition draw all derive from SEGMENT_BASE_SEED + segment*101.
Process noise is disabled (process_noise_cov=0): CDR has no exogenous
disturbance channel, matching the "d_k may be absent" acceptance criterion.
Sensor noise is small (std 0.01) so the linear dynamics dominate the fit.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.io import savemat

CONTROLMGYM_PATH = Path(
    os.environ.get("CDR_CONTROLMGYM_PATH", r"C:/tmp/controlgym")
).resolve()
sys.path.insert(0, str(CONTROLMGYM_PATH))

import controlgym  # noqa: E402

# Environment configuration (documented in provenance.json).
ENV_KWARGS = dict(
    n_steps=1001,  # env will terminate after 1001 steps; we only run 1000
    sample_time=0.01,  # 10 Hz sampling; keeps the DC reaction mode bounded
    process_noise_cov=0.0,  # no exogenous disturbance channel
    sensor_noise_cov=1e-4,  # std 0.01 per sensor, field is O(1)
    n_state=200,
    n_observation=30,
    n_action=8,
    control_sup_width=0.1,
    action_limit=None,
    observation_limit=None,
    reward_limit=None,
)
STEPS_PER_SEGMENT = 1001  # stored columns per segment (u_{T-1} is sentinel)
TRAIN_SEGMENTS = (1, 2, 3, 4, 5, 6)
VALIDATION_SEGMENTS = (7, 8)
PRBS_AMPLITUDE = 1.0
PRBS_HOLD = 4
SEGMENT_BASE_SEED = 20260801

OUTPUT_NAMES = tuple(f"sensor_{idx:02d}" for idx in range(ENV_KWARGS["n_observation"]))
INPUT_NAMES = tuple(f"u_{idx:02d}" for idx in range(ENV_KWARGS["n_action"]))


def prbs_sequence(segment: int, length: int) -> np.ndarray:
    """Deterministic two-level PRBS, held PRBS_HOLD samples, amplitude +/-1."""
    rng = np.random.default_rng(SEGMENT_BASE_SEED + segment * 101)
    levels = rng.integers(0, 2, size=(int(np.ceil(length / PRBS_HOLD)) + 1, ENV_KWARGS["n_action"]))
    held = np.repeat(levels, PRBS_HOLD, axis=0)[:length]
    return PRBS_AMPLITUDE * (2.0 * held - 1.0)


def collect_segment(segment: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (y, u) with STEPS_PER_SEGMENT columns; u[:,k] drives y[:,k+1]."""
    seed = SEGMENT_BASE_SEED + segment * 101
    env = controlgym.make(
        "convection_diffusion_reaction", **ENV_KWARGS, seed=seed
    )
    T = STEPS_PER_SEGMENT
    u = prbs_sequence(segment, T).T  # (n_action, T); u[:,k] drives y[:,k+1]
    y = np.zeros((ENV_KWARGS["n_observation"], T))
    obs, _ = env.reset()  # y_0
    y[:, 0] = obs
    for k in range(1, T):  # steps 1..T-1 apply u_0..u_{T-2}
        obs, _, _, _, _ = env.step(u[:, k - 1])  # returns y_k
        y[:, k] = obs
    return y, u


def build_dataset() -> dict:
    segments = TRAIN_SEGMENTS + VALIDATION_SEGMENTS
    columns = []
    segment_ids = []
    for segment in segments:
        y, u = collect_segment(segment)
        columns.append((y, u))
        segment_ids.extend([segment] * y.shape[1])
    y = np.concatenate([c[0] for c in columns], axis=1)
    u = np.concatenate([c[1] for c in columns], axis=1)
    segment_id = np.asarray(segment_ids, dtype=np.int32)
    return {
        "y": y,
        "u": u,
        "segment_id": segment_id,
        "output_names": np.asarray(OUTPUT_NAMES, dtype=object),
        "input_names": np.asarray(INPUT_NAMES, dtype=object),
    }


def controlgym_revision() -> str:
    try:
        import subprocess

        return subprocess.check_output(
            ["git", "-C", str(CONTROLMGYM_PATH), "rev-parse", "HEAD"],
            text=True,
        ).strip()
    except Exception:  # pragma: no cover - best effort
        return "unknown"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=None, help="output .mat path")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    out_path = Path(args.out) if args.out else here / "data" / "copyAX_cdr_dataset.mat"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    started = time.time()
    dataset = build_dataset()
    savemat(out_path, dataset)

    provenance = {
        "source": "ControlGym convection_diffusion_reaction (linear analytical CDR)",
        "controlgym_revision": controlgym_revision(),
        "controlgym_path": str(CONTROLMGYM_PATH),
        "paper": "Zhang et al., ControlGym: Large-Scale Neuro-Evolutionary PDE Control",
        "env_kwargs": {k: (None if v is None else v) for k, v in ENV_KWARGS.items()},
        "steps_per_segment": STEPS_PER_SEGMENT,
        "train_segments": list(TRAIN_SEGMENTS),
        "validation_segments": list(VALIDATION_SEGMENTS),
        "transitions_per_segment": STEPS_PER_SEGMENT - 1,
        "prbs_amplitude": PRBS_AMPLITUDE,
        "prbs_hold": PRBS_HOLD,
        "segment_seed_formula": "SEGMENT_BASE_SEED + segment*101, base=20260801",
        "alignment": "u[:,k] drives y[:,k+1] (standard causality; "
        "controlgym returns pre-update observation)",
        "sensor_noise_std": float(np.sqrt(ENV_KWARGS["sensor_noise_cov"])),
        "process_noise": "disabled",
        "collected_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "collector": "collect_cdr_data.py (copyAX)",
    }
    provenance_path = out_path.with_suffix(".json")
    provenance_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8")

    print(
        f"CDR dataset written: {out_path} "
        f"y={dataset['y'].shape} u={dataset['u'].shape} "
        f"segments={np.unique(dataset['segment_id']).tolist()} "
        f"elapsed={time.time()-started:.1f}s"
    )


if __name__ == "__main__":
    main()
