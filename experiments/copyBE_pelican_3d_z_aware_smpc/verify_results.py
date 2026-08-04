"""Independent numerical gates for copyBE, including held-out rollout errors."""

import json
from pathlib import Path

import h5py
import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "copyAY_pelican_soft_preference" / "data" / "copyAY_pelican_dataset.mat"
MODEL_BE = HERE / "results" / "copyBE_pelican_3d_z_aware_data.mat"
MODEL_BD = HERE.parent / "copyBD_pelican_3d_fixed_anchor_smpc" / "results" / "copyBD_pelican_3d_fixed_anchor_data.mat"
CLOSED_LOOP = HERE / "results" / "copyBE_3d_z_aware_trajectory.json"
AUDIT_JSON = HERE / "results" / "copyBE_multistep_audit.json"


def _read_h5(path, names):
    with h5py.File(path, "r") as f:
        model = {name: np.asarray(f["model"][name]).T for name in names}
        E = np.asarray(f["Eanchor"]).T
        norm = {name: np.asarray(f[name]).ravel() for name in
                ("y_offset", "y_scale", "u_offset", "u_scale")}
    return model, E, norm


def _valid_origins(segment_id, horizon, second_order=False):
    start = 1 if second_order else 0
    origins = np.arange(start, segment_id.size - horizon)
    return origins[segment_id[origins] == segment_id[origins + horizon]]


def rollout_rmse(model, E, y, u, segment_id, horizons, order):
    yc = y - model["y_mean"]
    uc = u - model["u_mean"]
    z_measured = model["R"].T @ yc
    output = {}
    for horizon in horizons:
        origins = _valid_origins(segment_id, horizon, order == 2)
        if order == 1:
            state = z_measured[:, origins].copy()
        else:
            state = np.vstack([z_measured[:, origins], z_measured[:, origins - 1]])
        for step in range(horizon):
            state = model["A"] @ state + model["B"] @ uc[:, origins + step]
        pred = model["y_mean"] + model["P"] @ state
        task_pred = E.T @ pred
        task_true = E.T @ y[:, origins + horizon]
        rmse = np.sqrt(np.mean((task_pred - task_true) ** 2, axis=1))
        output[str(horizon)] = {"origins": int(origins.size), "rmse": rmse.tolist()}
    return output


def main():
    raw = loadmat(DATA)
    segment_id = raw["segment_id"].ravel().astype(int)
    validation = np.isin(segment_id, np.arange(37, 55))
    sid = segment_id[validation]

    be, E, norm = _read_h5(MODEL_BE, ("A", "B", "P", "R", "y_mean", "u_mean"))
    bd, E_bd, norm_bd = _read_h5(MODEL_BD, ("A", "B", "P", "R", "y_mean", "u_mean"))
    assert np.allclose(E, E_bd)
    assert all(np.allclose(norm[key], norm_bd[key]) for key in norm)

    y = ((raw["y"][:, validation] - norm["y_offset"][:, None]) /
         norm["y_scale"][:, None])
    u = ((raw["u"][:, validation] - norm["u_offset"][:, None]) /
         norm["u_scale"][:, None])
    horizons = (10, 50, 100)
    audit = {
        "definition": "held-out recorded-input rollout; true output only initializes each origin",
        "horizons": list(horizons),
        "copyBD_order1": rollout_rmse(bd, E, y, u, sid, horizons, 1),
        "copyBE_order2": rollout_rmse(be, E, y, u, sid, horizons, 2),
    }
    AUDIT_JSON.write_text(json.dumps(audit, indent=2), encoding="ascii")

    for horizon in horizons:
        key = str(horizon)
        old = np.asarray(audit["copyBD_order1"][key]["rmse"])
        new = np.asarray(audit["copyBE_order2"][key]["rmse"])
        assert np.all(np.isfinite(new))
        if horizon <= 50:
            assert new[2] < old[2], (
                f"z rollout did not improve at h={horizon}: {new[2]} >= {old[2]}"
            )
        else:
            assert new[2] <= 1.01 * old[2], (
                f"z rollout degraded by more than 1% at h={horizon}: "
                f"{new[2]} vs {old[2]}"
            )
        print(f"h={horizon}: copyBD={old} copyBE={new}")

    with h5py.File(MODEL_BE, "r") as f:
        fit = f["fitstats"]
        spectral_radius = float(np.asarray(fit["spectral_radius"]).ravel()[0])
        dual_error = float(np.asarray(fit["dual_error"]).ravel()[0])
        anchor_error = float(np.asarray(fit["anchor_preservation_error"]).ravel()[0])
        reach = np.asarray(fit["all_input_horizon_gain"]).ravel()
    assert spectral_radius < 1.0
    assert dual_error < 1e-10
    assert anchor_error < 1e-10
    assert reach[2] > 2.0 * min(reach[:2])

    closed = json.loads(CLOSED_LOOP.read_text(encoding="ascii"))
    for seed, result in closed["per_seed"].items():
        assert result["smpc_violation_steps"] == 0
        assert result["mpc_violation_steps"] == 0
        assert result["smpc_fallback_count"] == 0
        assert result["mpc_fallback_count"] == 0
        assert result["smpc_rmse"][2] < 0.02, (seed, result["smpc_rmse"])
    print(f"spectral_radius={spectral_radius:.12f} reach_xyz={reach}")
    print("copyBE independent audit PASS")


if __name__ == "__main__":
    main()
