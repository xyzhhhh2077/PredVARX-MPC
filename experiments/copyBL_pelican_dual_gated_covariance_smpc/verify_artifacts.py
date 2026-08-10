"""Fail-closed validation for copyBL static figures, event log, and GIF."""

import hashlib
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageChops


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
MEDIA = HERE / "results_media"


def main():
    reject_constant = lambda token: (_ for _ in ()).throw(ValueError(f"non-finite JSON: {token}"))
    summary_path = RESULTS / "copyBL_dual_gated_summary.json"
    events_path = RESULTS / "copyBL_dual_gate_events.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"), parse_constant=reject_constant)
    events = json.loads(events_path.read_text(encoding="utf-8"), parse_constant=reject_constant)
    manifest = json.loads((MEDIA / "copyBL_figure_manifest.json").read_text(encoding="utf-8"))
    figures = [MEDIA / name for name in manifest["figures"].values()]
    if len(figures) != 4 or any(not path.exists() or path.stat().st_size < 20_000 for path in figures):
        raise AssertionError("Expected four nontrivial static figures")
    gate = summary["copyBL"]
    if int(summary["steps"]) != 18_000 or gate["completed_steps"] != 18_000:
        raise AssertionError("Artifact source must be the complete 18,000-step run")
    if gate["qp_failure_count"] or gate["fallback_count"] or gate["hard_violation_steps"]:
        raise AssertionError("Complete run contains a QP failure, fallback, or hard violation")
    event_hash = hashlib.sha256(events_path.read_bytes()).hexdigest()
    if event_hash != summary["event_log_sha256"]:
        raise AssertionError("Event log hash does not match the run summary")
    if gate["feasibility_trigger_count"] != len(events):
        raise AssertionError("Event log and summary trigger counts differ")
    if gate["value_trigger_count"] <= 0:
        raise AssertionError("Control-value gate never triggered")
    value_events = [event for event in events if event["value_triggered"]]
    observed = {
        "feasibility_rejection_count": sum(not event["feasibility_accepted"] for event in events),
        "value_trigger_count": len(value_events),
        "value_accept_count": sum(event["value_accepted"] for event in value_events),
        "value_rejection_count": sum(not event["value_accepted"] for event in value_events),
        "committed_update_count": sum(event["update_committed"] for event in events),
    }
    if any(gate[key] != value for key, value in observed.items()):
        raise AssertionError("Event-derived gate counts differ from the run summary")
    if any(event["value_accepted"] != event["update_committed"] for event in events):
        raise AssertionError("Value acceptance and atomic commit are inconsistent")
    data = np.load(RESULTS / "copyBL_dual_gated_run.npz")
    required_arrays = (
        "reference_standardized", "trajectory_standardized", "control",
        "residual_innovation_latent", "wind_velocity_mps", "stage_cost",
        "estimated_sigma_eps_std", "accepted_sigma_eps_std",
    )
    if any(len(data[name]) != 18_000 for name in required_arrays):
        raise AssertionError("NPZ does not contain a complete aligned 18,000-step run")
    gif_path = MEDIA / "copyBL_dual_gate_summary_200frames.gif"
    image = Image.open(gif_path)
    if image.n_frames != 200:
        raise AssertionError("GIF must contain exactly 200 frames")
    image.seek(0); first = image.convert("RGB").copy()
    image.seek(199); last = image.convert("RGB").copy()
    if ImageChops.difference(first, last).getbbox() is None:
        raise AssertionError("GIF first and final frames are pixel-identical")
    report = {
        "static_figure_count": 4,
        "static_figure_bytes": {path.name: path.stat().st_size for path in figures},
        "gif_frames": image.n_frames,
        "gif_size_pixels": list(image.size),
        "gif_bytes": gif_path.stat().st_size,
        "feasibility_trigger_count": gate["feasibility_trigger_count"],
        "feasibility_rejection_count": gate["feasibility_rejection_count"],
        "value_trigger_count": gate["value_trigger_count"],
        "value_accept_count": gate["value_accept_count"],
        "value_rejection_count": gate["value_rejection_count"],
        "committed_update_count": gate["committed_update_count"],
    }
    (MEDIA / "copyBL_artifact_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()