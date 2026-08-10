import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "experiments" / "copyBM_pelican_20seed_fair_comparison" / "run_fair_comparison.py"


def load_module(name="copybm_contract"):
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_fixed_seed_pairs_are_unique_and_not_result_selected():
    bm = load_module()
    pairs = bm.seed_pairs()
    assert len(pairs) == 20
    assert len(set(pairs)) == 20
    assert pairs == [(16 + i, 7016 + i) for i in range(20)]


def test_same_pair_disturbance_sha_contract():
    bm = load_module()
    arrays = bm.make_pair_inputs(8, 16, 7016)
    assert len({bm.disturbance_sha256(arrays["innovation"], arrays["wind_latent_increment"])
                for _ in ("Frozen", "BJ", "BL")}) == 1
    assert arrays["initial_task"].tolist() == [0.0, 0.0, 0.0]


def test_cost_uses_realized_stage_cost_and_requires_complete_run(monkeypatch):
    bm = load_module()
    seen = {}
    def fake_cost(S, ref, U, mean):
        seen["args"] = (S, ref, U, mean)
        return np.array([2., 3.])
    monkeypatch.setattr(bm.BI, "realized_stage_cost", fake_cost)
    result = {"S": np.ones((2, 3)), "U": np.ones((2, 1)), "completed_steps": 2}
    assert bm.controller_metrics(result, np.zeros((2, 3)), np.zeros(1), 2)["cumulative_realized_stage_cost"] == 5.0
    assert seen["args"][1].shape == (2, 3)
    incomplete = {"S": np.ones((1, 3)), "U": np.ones((1, 1)), "completed_steps": 1}
    assert bm.controller_metrics(incomplete, np.zeros((2, 3)), np.zeros(1), 2)["cumulative_realized_stage_cost"] is None


def test_bl_value_acceptance_rate_from_events():
    bm = load_module()
    events = [{"feasibility_triggered": True, "feasibility_accepted": False, "value_triggered": False,
               "value_accepted": False, "update_committed": False},
              {"feasibility_triggered": True, "feasibility_accepted": True, "value_triggered": True,
               "value_accepted": True, "update_committed": True},
              {"feasibility_triggered": True, "feasibility_accepted": True, "value_triggered": True,
               "value_accepted": False, "update_committed": False}]
    m = bm.gate_metrics(events)
    assert m == {"feasibility_trigger_count": 3, "feasibility_rejection_count": 1,
                 "value_trigger_count": 2, "value_accept_count": 1,
                 "value_acceptance_rate": 0.5, "value_rejection_count": 1,
                 "update_committed_count": 1}


def test_strict_json_and_incomplete_metrics_are_null(tmp_path):
    bm = load_module()
    payload = {"x": float("nan"), "y": [float("inf")]}
    text = bm.strict_json_dumps(payload)
    assert json.loads(text) == {"x": None, "y": [None]}
    assert bm.controller_metrics({"S": np.zeros((0, 3)), "U": np.zeros((0, 1)), "completed_steps": 0},
                                 np.zeros((2, 3)), np.zeros(1), 2)["rmse_xyz"] is None


def test_resume_rejects_configuration_mismatch(tmp_path):
    bm = load_module()
    path = tmp_path / "per_seed.json"
    path.write_text(json.dumps({"config_sha256": "wrong", "seed_pairs": bm.seed_pairs(), "seeds": []}))
    with pytest.raises(ValueError, match="configuration"):
        bm.load_resume(path, "actual")
