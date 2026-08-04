"""Acceptance checks for the copyBF long boundary-disturbance experiment."""

import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT = HERE / "results" / "copyBF_long_3d_boundary_disturbance.json"


def main():
    data = json.loads(RESULT.read_text(encoding="ascii"))
    assert data["scenario_steps"] == 4000
    assert data["scenario_duration_seconds"] == 40.0
    assert data["total_plotted_duration_seconds"] == 120.0
    assert data["scenario_resets"] == [40.0, 80.0]
    assert data["noise_statement"].startswith("identified residual innovation")

    failures = []
    for name, scenario in data["scenarios"].items():
        for seed, result in scenario["per_seed"].items():
            smpc = result["smpc"]
            mpc = result["deterministic_mpc"]
            if smpc["fallback_count"] != 0:
                failures.append(f"{name}/seed{seed}: SMPC fallback={smpc['fallback_count']}")
            if mpc["fallback_count"] != 0:
                failures.append(f"{name}/seed{seed}: MPC fallback={mpc['fallback_count']}")
            if smpc["active_qp_steps"] <= 0 or smpc["positive_dual_qp_steps"] <= 0:
                failures.append(f"{name}/seed{seed}: chance constraint never active")

            if name != "z_near_boundary" and smpc["minimum_hard_margin"] <= mpc["minimum_hard_margin"]:
                failures.append(
                    f"{name}/seed{seed}: SMPC did not retain a larger hard-bound margin"
                )

    if failures:
        raise AssertionError("\n".join(failures))
    print("copyBF acceptance PASS")


if __name__ == "__main__":
    main()
