"""Independent audit for copyBB saved identification and closed-loop results."""

import json
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def main():
    model_path = RESULTS / "copyBB_pelican_fixed_anchor_data.mat"
    json_path = RESULTS / "copyBB_fixed_anchor_closed_loop.json"
    npz_path = RESULTS / "copyBB_fixed_anchor_closed_loop.npz"

    with h5py.File(model_path, "r") as f:
        E = np.asarray(f["Eanchor"]).T
        model = f["model"]
        P = np.asarray(model["P"]).T
        R = np.asarray(model["R"]).T
        A = np.asarray(model["A"]).T
        Sigma = np.asarray(model["Sigma_eps"]).T

    expected = np.zeros((10, 2))
    expected[0, 0] = 1.0
    expected[1, 1] = 1.0
    assert np.array_equal(E, expected), "Controlled directions are not fixed x/y axes"
    assert np.linalg.norm(R.T @ P - np.eye(P.shape[1])) < 1e-7
    assert np.linalg.norm(P @ R.T @ E - E) < 1e-10
    assert np.linalg.eigvalsh((Sigma + Sigma.T) / 2.0).min() > -1e-10

    saved = json.loads(json_path.read_text(encoding="ascii"))
    traces = np.load(npz_path, allow_pickle=True)["traces"].item()
    assert saved["training_flights"] == list(range(1, 37))
    assert saved["held_out_reference_flights"] == [37, 38, 39, 40, 41]
    assert saved["fixed_anchor"] == "standardized x/y position channels"

    improvements = []
    for segid, trace in traces.items():
        metric = saved["per_seg"][segid]
        assert trace["S_ref"].shape == (1200, 2)
        assert all(np.isfinite(array).all() for array in trace.values())
        rmse = np.sqrt(np.mean((trace["S_smpc"] - trace["S_ref"]) ** 2, axis=0))
        rmse_open = np.sqrt(np.mean((trace["S_open"] - trace["S_ref"]) ** 2, axis=0))
        improvement = 100.0 * (rmse_open - rmse) / rmse_open
        assert np.allclose(rmse, metric["rmse_smpc"])
        assert np.allclose(rmse_open, metric["rmse_openloop"])
        assert np.allclose(improvement, metric["improvement_pct"])
        assert metric["fallback_count"] == 0
        improvements.extend(improvement.tolist())

    print(f"dual_error {np.linalg.norm(R.T @ P - np.eye(P.shape[1])):.3e}")
    print(f"anchor_error {np.linalg.norm(P @ R.T @ E - E):.3e}")
    print(f"spectral_radius {max(abs(np.linalg.eigvals(A))):.10f}")
    print(f"minimum_open_loop_improvement_pct {min(improvements):.3f}")
    print("copyBB independent audit PASS")


if __name__ == "__main__":
    main()
