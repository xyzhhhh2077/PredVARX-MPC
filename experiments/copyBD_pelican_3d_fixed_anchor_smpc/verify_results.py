"""Independent audit of copyBD identification and 3D closed-loop artifacts."""

import json
from pathlib import Path

import h5py
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def main():
    with h5py.File(RESULTS / "copyBD_pelican_3d_fixed_anchor_data.mat", "r") as f:
        E = np.asarray(f["Eanchor"]).T
        m = f["model"]
        A = np.asarray(m["A"]).T
        B = np.asarray(m["B"]).T
        P = np.asarray(m["P"]).T
        R = np.asarray(m["R"]).T
        Sigma = np.asarray(m["Sigma_eps"]).T

    expected = np.zeros((10, 3))
    expected[:3, :3] = np.eye(3)
    assert np.array_equal(E, expected)
    assert np.linalg.norm(R.T @ P - np.eye(5)) < 1e-7
    assert np.linalg.norm(P @ R.T @ E - E) < 1e-10
    assert np.linalg.eigvalsh((Sigma + Sigma.T) / 2.0).min() > -1e-10
    ctr = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(5)])
    assert np.linalg.matrix_rank(ctr) == 5

    saved = json.loads((RESULTS / "copyBD_3d_boundary_trajectory.json").read_text(encoding="ascii"))
    arrays = np.load(RESULTS / "copyBD_3d_boundary_trajectory.npz")
    reference, smpc, mpc = arrays["reference"], arrays["smpc"], arrays["mpc"]
    assert reference.shape == smpc.shape == mpc.shape == (1200, 3)
    assert np.allclose(np.ptp(reference, axis=0), [0.5, 1.1, 0.64], atol=0.01)
    assert np.linalg.matrix_rank(reference - np.mean(reference, axis=0), tol=1e-8) == 3
    assert np.all(np.ptp(smpc, axis=0) > np.array([0.45, 0.95, 0.35]))
    assert all(item["smpc_fallback_count"] == 0 for item in saved["per_seed"].values())
    assert all(item["mpc_fallback_count"] == 0 for item in saved["per_seed"].values())
    assert all(item["smpc_violation_steps"] == 0 for item in saved["per_seed"].values())
    assert all(np.all(np.asarray(item["smpc_rmse"]) < [0.08, 0.15, 0.23])
               for item in saved["per_seed"].values())

    image = np.asarray(Image.open(RESULTS / "copyBD_3d_boundary_trajectory.png").convert("RGB"))
    assert image.shape[:2] == (1620, 2700)
    expected_colors = {
        "black_reference": np.array([0, 0, 0]),
        "teal_smpc": np.array([8, 126, 139]),
        "red_mpc": np.array([196, 69, 54]),
    }
    color_counts = {}
    for name, rgb in expected_colors.items():
        color_counts[name] = int(np.all(np.abs(image.astype(int) - rgb) <= 3, axis=2).sum())
        assert color_counts[name] > 100
    nonwhite = np.any(image < 245, axis=2)
    assert nonwhite[:, :1350].sum() > 10000
    assert nonwhite[:, 1350:].sum() > 10000

    sat = 100.0 * np.mean(np.max(np.abs(arrays["smpc_u"]), axis=1) >= 5.999)
    assert sat == 0.0
    print(f"dual_error {np.linalg.norm(R.T @ P - np.eye(5)):.3e}")
    print(f"anchor_error {np.linalg.norm(P @ R.T @ E - E):.3e}")
    print(f"spectral_radius {max(abs(np.linalg.eigvals(A))):.10f}")
    print(f"reference_ptp {np.ptp(reference, axis=0)}")
    print(f"smpc_ptp {np.ptp(smpc, axis=0)}")
    print(f"seed7_input_saturation_pct {sat:.3f}")
    print(f"figure_color_pixels {color_counts}")
    print("copyBD independent audit PASS")


if __name__ == "__main__":
    main()
