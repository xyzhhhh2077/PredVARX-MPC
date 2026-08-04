# copyBA — Waterloo Pelican: SMPC closed loop on the REAL identified model

Variant of copyAY with a **position-dominant task axis**, then a closed-loop
SMPC experiment **in-the-loop**: the plant is the model identified from real
Waterloo Pelican flight data (training segments 1:36), and the SMPC controller
(the centered_smpc_step port) drives it to reproduce **real flight segment 37**
in the task space.

## Why copyBA exists

copyAY's `learn_segmented_preferred_output_directions` with weights
[1.0; 0.9; 0.7] and preference_strength = 0.78 learned a **motor-dominant** E
(physical weight: rpm ±3.8 vs position ±0.03). On the Pelican data the motor
channels have ~6x the variance of position, and the input-reachability Gram is
motor-dominated, so the preference weights could not pull E to position.
Consequence: an s1 reference of 0.8 was physically a motor-speed-pattern
change, not a position task — the closed loop had nothing meaningful to track.

copyBA re-learns E with a strong position preference
(weights pos=1.0, euler=0.3, motors=0.05, preference_strength=0.99) so the
task axis is a **position task** (E physical weight: position 0.51–0.71 m
dominant). This is the "task-needed output vs input-pushed direction"
distinction: the input-reachability direction is motor, but the task is
position.

## Results

### Offline validation (copyBA_pelican_position_task_metrics.txt)
- Task-axis prediction on validation segments 37:54: task RMSE = [0.020, 0.027],
  task R² = [0.9994, 0.999] — position-dominant task axis is still very well
  predicted by the same latent VARX model.

### Closed loop, model in-the-loop (copyBA_multi_segment_Ru006_fig.png)
- Plant: copyBA model (A, B, P, R, Sigma_eps) with process noise per 100 Hz step.
- Reference: each validation segment's own real flight (segments 37-41, first
  1200 samples each), task axis s = E' y_std. Initial state = real segment start.
- Multirate control: QP solved every d = 10 steps (10 Hz decisions),
  N = 18 horizon → 1.8 s, matching the position time scale (A has 4 unit roots).
- Q = diag([1×6, 0.05×4])/6 × 10, Ru = 0.06 (input magnitude matched to the
  real expert: u std within ~15% of expert on every segment), u ∈ [−6, 6] std,
  task bound h = 2.0, alpha_joint = 0.05.
- Lifted noise: Sigma_d = Σ_{i=0}^{d-1} A^i Σ_eps A'^i (d-step disturbance
  variance, so a decision-step prediction carries the full d-step noise).
- Results (5 segments, 6000 steps total):
  | seg | RMSE s1 | RMSE s2 | u sat | violations |
  |---|---|---|---|---|
  | 37 | 0.589 | 0.452 | 0% | 0 |
  | 38 | 0.209 | 0.492 | 0% | 0 |
  | 39 | 0.780 | 0.713 | 0% | 0 |
  | 40 | 0.529 | 0.432 | 0% | 0 |
  | 41 | 0.559 | 0.301 | 0% | 0 |
- Open-loop baseline (u = u_mean hold, same noise stream), seg37:
  RMSE (0.841, 1.113) — SMPC cuts tracking error roughly in half.
- Input magnitude: SMPC u std per channel ≈ expert u std per channel on all
  5 segments (Ru = 0.06); Ru = 0.01 gives ~4× expert input but better tracking
  (seg37 RMSE 0.426/0.243) — input-matching costs some tracking accuracy.

### Diagnostic findings along the way
1. N = 18 at 100 Hz decisions was too short for the Pelican position dynamics
   (4 integrator modes; s1 moves ~0.007/step under full input) — QP saw little
   horizon gain and backed off. Multirate d = 10 fixes this.
2. Starting from a real segment start (s2 = −1.22) while referencing s2 = 0
   splits the controller between recovery and tracking; referencing the real
   trajectory itself removes the conflict and is the fair "reproduce real
   flight" test.
3. Ru = 1e-4 (copyAZ heritage) gives 84% input saturation; Ru = 0.06 matches
   the expert's input magnitude with 0% saturation.
4. The SMPC controller must know the multirate disturbance variance:
   using the single-step Sigma_eps inside the lifted model understates the
   chance-constraint kappa (nominal α = 5% would not be the true violation
   probability). Sigma_d accumulation fixes this; results unchanged because the
   task bound h = 2.0 was never active (viol = 0 in all runs).

## Files
- copyBA_pelican_position_task.m — identify + offline validation
- run_real_model_closed_loop.py — in-the-loop closed-loop experiment
- plot_track_real_seg37.py — figure
- results/copyBA_pelican_position_task_data.mat — model + Etask + scales
- results/copyBA_track_real_seg37_fig.png — closed-loop figure
- results/copyBA_closed_loop_metrics.txt — closed-loop numbers
