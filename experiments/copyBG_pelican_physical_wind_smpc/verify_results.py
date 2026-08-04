"""Hard acceptance checks for copyBG saved artifacts."""

import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
STEM = "copyBG_pelican_physical_wind_smpc"


def _scalar(value):
    return float(np.asarray(value).squeeze())


def verify(results_dir=RESULTS):
    results_dir = Path(results_dir)
    json_path = results_dir / f"{STEM}.json"
    npz_path = results_dir / f"{STEM}.npz"
    mat_path = results_dir / f"{STEM}.mat"
    png_paths = [
        results_dir / f"{STEM}_smpc.png",
        results_dir / f"{STEM}_deterministic_mpc.png",
    ]
    for path in [json_path, npz_path, mat_path, *png_paths]:
        assert path.is_file() and path.stat().st_size > 0, path

    summary = json.loads(json_path.read_text(encoding="ascii"))
    arrays = np.load(npz_path)
    mat = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    required_arrays = {
        "reference_standardized", "smpc_trajectory_standardized",
        "deterministic_mpc_trajectory_standardized", "innovation_latent",
        "wind_velocity_mps", "wind_position_delta_xyz_m",
        "wind_latent_increment", "Sigma_eps", "Sigma_obs_proxy",
        "Sigma_eps_aug", "Sigma_wind_aug", "Sigma_total_aug", "G_w",
        "dt_used_seconds", "shared_noise_exact",
    }
    assert required_arrays <= set(arrays.files)

    requested_steps = int(summary["requested_steps"])
    assert arrays["reference_standardized"].shape == (requested_steps, 3)
    assert arrays["innovation_latent"].shape == (requested_steps, 5)
    assert arrays["wind_velocity_mps"].shape == (requested_steps, 3)
    assert arrays["wind_latent_increment"].shape == (requested_steps, 5)
    assert summary["position_unit"] == "m"
    assert abs(float(summary["dt_used_seconds"]) - 0.01) < 1e-15
    assert _scalar(mat["dt_used_seconds"]) == 0.01
    assert int(summary["training_segments"]) == 36
    assert int(summary["uses_true_Sigma_n"]) == 0
    assert int(summary["shared_noise_exact"]) == 1
    assert int(np.asarray(arrays["shared_noise_exact"]).squeeze()) == 1
    assert int(_scalar(mat["shared_noise_exact"])) == 1

    assert np.array_equal(
        arrays["wind_position_delta_xyz_m"],
        arrays["wind_velocity_mps"] * 0.01,
    )
    y_scale = np.asarray(mat["position_scale_m"], dtype=float).reshape(3)
    standardized_delta = arrays["wind_velocity_mps"] * 0.01 / y_scale
    reconstructed = (arrays["G_w"][:3] @ arrays["wind_velocity_mps"].T).T
    assert np.allclose(reconstructed, standardized_delta, atol=1e-12)
    assert np.allclose(
        arrays["Sigma_total_aug"],
        arrays["Sigma_eps_aug"] + arrays["Sigma_wind_aug"],
        atol=1e-15,
    )
    assert np.allclose(arrays["Sigma_eps"], mat["Sigma_eps"], atol=1e-15)
    assert np.allclose(
        arrays["Sigma_obs_proxy"], mat["Sigma_obs_proxy"], atol=1e-15,
    )
    assert np.allclose(
        arrays["wind_velocity_mps"], mat["wind_velocity_mps"], atol=0.0,
    )
    assert np.allclose(
        arrays["innovation_latent"], mat["innovation_latent"], atol=0.0,
    )

    smpc = summary["smpc"]
    deterministic = summary["deterministic_mpc"]
    for result in (smpc, deterministic):
        assert int(result["fallback_count"]) == 0
        assert float(result["maximum_qp_constraint_residual"]) <= 1e-6
    assert int(smpc["completed_steps"]) == requested_steps
    assert int(smpc["qp_failure_count"]) == 0
    # Chance-layer activity is only required when the reference presses near the
    # hard bound. Interior trackable orbits (copyBE-scale) may stay inactive.
    peak_ref = float(np.max(np.abs(arrays["reference_standardized"])))
    hard_bound = float(summary.get("hard_bound_standardized", 3.4))
    if peak_ref >= 0.70 * hard_bound:
        assert int(smpc["active_qp_steps"]) > 0
        assert int(smpc["positive_dual_qp_steps"]) > 0
    assert float(smpc["hard_violation_rate"]) <= 0.05
    # Trackable interior path should not wander far from the reference.
    if peak_ref <= 1.0:
        rmse = np.asarray(smpc["rmse_standardized"], dtype=float).reshape(3)
        assert float(np.max(rmse)) < 0.25, rmse
    assert summary["disturbance_sha256"] == smpc["disturbance_sha256"]
    assert smpc["disturbance_sha256"] == deterministic["disturbance_sha256"]
    assert arrays["smpc_trajectory_standardized"].shape[0] == requested_steps
    assert arrays["deterministic_mpc_trajectory_standardized"].shape[0] == int(
        deterministic["completed_steps"]
    )
    if int(deterministic["completed_steps"]) < requested_steps:
        assert int(deterministic["qp_failure_count"]) == 1
        assert int(deterministic["qp_failure_step"]) == int(
            deterministic["completed_steps"]
        )

    image_stats = {}
    for path in png_paths:
        with Image.open(path) as image:
            image.verify()
        with Image.open(path).convert("RGB") as image:
            pixels = np.asarray(image)
        assert pixels.shape[0] >= 1000 and pixels.shape[1] >= 1000
        assert float(np.std(pixels)) > 10.0
        assert np.unique(pixels.reshape(-1, 3), axis=0).shape[0] > 100
        image_stats[path.name] = {
            "shape": list(pixels.shape),
            "pixel_std": float(np.std(pixels)),
        }

    return {
        "status": "PASS",
        "requested_steps": requested_steps,
        "smpc_completed_steps": int(smpc["completed_steps"]),
        "deterministic_mpc_completed_steps": int(deterministic["completed_steps"]),
        "smpc_active_qp_steps": int(smpc["active_qp_steps"]),
        "smpc_hard_violation_rate": float(smpc["hard_violation_rate"]),
        "png": image_stats,
    }


def main():
    print(json.dumps(verify(), indent=2))


if __name__ == "__main__":
    main()
