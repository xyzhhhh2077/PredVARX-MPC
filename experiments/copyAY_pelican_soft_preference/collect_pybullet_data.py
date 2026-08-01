"""Collect excitation data from gym-pybullet-drones (CF2X) for copyAY pipeline.

Uses CtrlAviary (direct RPM command). Observation (20,) = [pos(3), quat(4),
rpy(3), vel(3), ang_v(3), rpm(4)] -- the last 4 entries are the motor RPMs,
so y = [pos(3); rpy(3); rpm(4)] (10 channels) mirrors the Pelican copyAY
structure Pos+Euler+Motors exactly.

Purpose: closed-loop phase-2 validation of the copyAU/copyAY pipeline on a
dynamic simulator. The simulator dynamics differ from Pelican (CF2X, 33 g vs
Pelican 1.3 kg) -- this validates the METHOD end-to-end (identify -> latent
VARX -> SMPC), not the Pelican plant.

Usage:
  python collect_pybullet_data.py --out data/copyAZ_pybullet_dataset.mat \
      --seed 7 --ctrl_freq 100 --n_seconds 60 --n_flights 6
"""
import argparse
import json
import os
import time

import numpy as np

from gym_pybullet_drones.envs.CtrlAviary import CtrlAviary
from gym_pybullet_drones.utils.enums import DroneModel, Physics

HOVER_RPM = 22000.0  # CF2X hover RPM (gym-pybullet-drones constant)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data/copyAZ_pybullet_dataset.mat")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--ctrl_freq", type=int, default=100)
    ap.add_argument("--n_seconds", type=float, default=60.0)
    ap.add_argument("--style", default="chirp",
                    choices=["chirp", "random", "square"])
    ap.add_argument("--n_flights", type=int, default=6)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    n_steps = int(args.n_seconds * args.ctrl_freq)

    if args.style == "chirp":
        # gentle swept-sine around hover: slow, small amplitude, per-motor
        # phase offsets produce pitch/roll/yaw excitation without tumbling.
        # f: 0.1 -> 1.2 Hz over the flight; mod amplitude 5% of hover.
        u_cmds = np.zeros((args.n_flights, n_steps, 4))
        for f in range(args.n_flights):
            for k in range(n_steps):
                t = k / args.ctrl_freq
                fHz = 0.1 + (1.1 * t / args.n_seconds)
                phi = np.array([0.0, np.pi / 2, np.pi, 3 * np.pi / 2])
                mod = 1.0 + 0.05 * np.sin(2 * np.pi * fHz * t + phi)
                u_cmds[f, k] = mod * HOVER_RPM
    elif args.style == "random":
        # persistent AR(1) modulation around hover, SMALL step noise so the
        # open-loop CF2X does not tumble (previous 0.02 std -> 40 m drift).
        mods = np.clip(0.95 + 0.05 * rng.random((args.n_flights, 4)), 0.9, 1.1)
        u_cmds = np.zeros((args.n_flights, n_steps, 4))
        for f in range(args.n_flights):
            m = mods[f]
            for k in range(n_steps):
                m = 0.997 * m + 0.004 * rng.normal(0, 1, 4)
                m = np.clip(m, 0.9, 1.1)
                u_cmds[f, k] = m * HOVER_RPM
    elif args.style == "square":
        phase = (np.arange(n_steps) / args.ctrl_freq // 2.0) % 2.0
        u_cmds = HOVER_RPM * np.where(phase[:, None] < 1.0, 1.05, 0.95)
        u_cmds = np.broadcast_to(u_cmds[None], (args.n_flights, n_steps, 4)).copy()
    else:
        raise ValueError(args.style)

    env = CtrlAviary(drone_model=DroneModel.CF2X,
                     num_drones=1,
                     physics=Physics.PYB,
                     pyb_freq=args.ctrl_freq * 4,
                     ctrl_freq=args.ctrl_freq,
                     gui=False,
                     user_debug_gui=False)
    # run a few settle-in steps at hover so the drone is airborne & stable
    env.reset()

    y_all, u_all, seg = [], [], []
    t0 = time.time()
    for f in range(args.n_flights):
        env.reset()
        # settle at hover for 2 s
        for _ in range(2 * args.ctrl_freq):
            env.step(np.array([[HOVER_RPM] * 4]))
        for k in range(n_steps):
            u = u_cmds[f, k]
            obs, *_ = env.step(u.reshape(1, 4))
            # obs (1,20): pos(3) quat(4) rpy(3) vel(3) ang_v(3) rpm(4)
            o = obs.reshape(-1)
            y = np.hstack([o[0:3], o[7:10], o[16:20]])
            y_all.append(y)
            u_all.append(u)
            seg.append(f + 1)
        print(f"flight {f + 1}/{args.n_flights} done "
              f"({time.time() - t0:.1f}s elapsed)")
    env.close()

    y = np.asarray(y_all).T          # 12 x T
    u = np.asarray(u_all).T          # 4 x T
    seg = np.asarray(seg, dtype=np.int32)

    # verify no NaN, sane ranges, C_y conditioning
    assert np.all(np.isfinite(y)) and np.all(np.isfinite(u))
    print("y range:", np.round(y.min(1), 3), np.round(y.max(1), 3))
    yc = y - y.mean(1, keepdims=True)
    Cy = yc @ yc.T / y.shape[1]
    print("cond(C_y):", np.linalg.cond(Cy))
    print("u range:", u.min(), u.max())

    # provenance
    meta = {
        "simulator": "gym-pybullet-drones 2.1.0 / pybullet 3.2.7",
        "drone_model": "CF2X", "physics": "PYB", "aviary": "CtrlAviary",
        "ctrl_freq": args.ctrl_freq, "pyb_freq": args.ctrl_freq * 4,
        "n_flights": args.n_flights, "n_seconds": args.n_seconds,
        "style": args.style, "seed": args.seed,
        "settle_seconds": 2.0,
        "y_channels": ["pos_x", "pos_y", "pos_z", "rpy_r", "rpy_p", "rpy_y",
                       "rpm_m0", "rpm_m1", "rpm_m2", "rpm_m3"],
        "u_channels": ["rpm_cmd_0", "rpm_cmd_1", "rpm_cmd_2", "rpm_cmd_3"],
        "u_units": "RPM (commanded, 0.7-1.3 x 22000 hover)",
        "y_units": "m, rad, RPM",
        "note": "y = pos+rpy+measured motor RPM (obs[16:20]); mirrors Pelican "
                "copyAY structure Pos+Euler+Motors; C_y cond "
                f"{np.linalg.cond(Cy):.3f}",
        "collected": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    import scipy.io as sio
    sio.savemat(args.out, {"y": y, "u": u, "segment_id": seg})
    meta_path = args.out.replace(".mat", ".json")
    with open(meta_path, "w") as fp:
        json.dump(meta, fp, indent=2)
    print("saved", args.out, y.shape, u.shape, "->", meta_path)


if __name__ == "__main__":
    main()
