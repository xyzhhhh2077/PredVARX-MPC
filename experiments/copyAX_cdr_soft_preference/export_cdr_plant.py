"""Export the exact ControlGym CDR plant used by the copyAX controller."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import savemat

from collect_cdr_data import CONTROLMGYM_PATH, ENV_KWARGS, controlgym_revision

import controlgym

SEED = 20260801


def main() -> None:
    here = Path(__file__).resolve().parent
    env = controlgym.make(
        "convection_diffusion_reaction", **ENV_KWARGS, seed=SEED
    )
    observation, _ = env.reset()
    payload = {
        "A_plant": np.asarray(env.A, dtype=float),
        "B_plant": np.asarray(env.B2, dtype=float),
        "C_plant": np.asarray(env.C, dtype=float),
        "x0": np.asarray(env.state, dtype=float).reshape(-1, 1),
        "y0_noise_free": (np.asarray(env.C) @ np.asarray(env.state)).reshape(-1, 1),
        "y0_observed": np.asarray(observation, dtype=float).reshape(-1, 1),
        "sample_time": np.asarray([[ENV_KWARGS["sample_time"]]], dtype=float),
        "seed": np.asarray([[SEED]], dtype=np.int64),
    }
    out = here / "data" / "copyAX_cdr_plant.mat"
    savemat(out, payload)
    provenance = {
        "source": "ControlGym convection_diffusion_reaction",
        "controlgym_revision": controlgym_revision(),
        "controlgym_path": str(CONTROLMGYM_PATH),
        "seed": SEED,
        "matrices": {key: list(value.shape) for key, value in payload.items() if hasattr(value, "shape")},
        "dynamics": "x_next=A_plant*x+B_plant*u_raw; y=C_plant*x+v",
        "process_noise": "disabled",
    }
    out.with_suffix(".json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(
        f"CDR plant written: {out} A={payload['A_plant'].shape} "
        f"B={payload['B_plant'].shape} C={payload['C_plant'].shape}"
    )


if __name__ == "__main__":
    main()
