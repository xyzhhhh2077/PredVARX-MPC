# copyAN_crte_fixed_surrogate

## What this copy is

Implements the **executable fixed spectral surrogate** described in CRTE draft
Sec 3.3. It is not the profiled teacher objective:

- `x in R^(p-q)` chooses a free latent direction in `range(Nperp)`.
- For each candidate `(mu, alpha, beta)` in a frozen grid:
  - build complete metric-dual basis `(P, R, Pbar, Rbar)` from the same metric
    `C_mu = (1-mu)*Sigma_perp + mu*tau_G*G`;
  - compute `A_CRTE = Ntr(S_yu) + alpha*Ntr(A_T) - beta*Ntr(Sigma_noise_proxy)`
    with each term **Frobenius-normalized** to fix the trace-zero collapse that
    the paper `Ntr` produces;
  - solve the symmetric generalized problem in `C_mu` (no `+gamma*I` shift);
  - refit VARX on the candidate basis;
  - compute candidate-specific FWL task score, readout-noise score, and
    finite-horizon latent authority `c_x'W_c c_x` for the `s=1` companion.
- Apply outer gates:
  - `min(task_axis) >= task_gate * max(task_axis)`;
  - `max(noise_axis) <= noise_gate_factor * median(noise_axis)`;
  - `min(reach_axis) >= reach_gate * max(reach_axis)`;
  - `dual_error < 1e-8` and `rho(Ahat) < 1.05`.
- Pick the smallest validation one-step NRMSE among survivors.
- Rebuild the chosen candidate on **all offline data**, refit the complete
  dual basis, and run the standard centered SMPC.

## Honest claim boundary

- No true `Sigma_n` is read or accepted. The noise term is the one-step
  free-output residual covariance estimated from training data; it is a
  residual proxy, **not** a sensor-noise covariance.
- The fixed spectral surrogate does not solve the profiled teacher objective.
  No consistency theorem is offered.
- `cover` and `QP success rate` are warm-horizon **time-fraction proxies**.
  No recursive feasibility, stability, or chance-constraint probability
  certificate is claimed.
- Single-seed result. Robustness needs >= 30 seeds.

## Files

| File | Role |
|---|---|
| `crte_fixed_surrogate_varx.m` | identifier with grid search + per-candidate gates |
| `copyAN_crte_fixed_surrogate.m` | same-plant runner; produces `.mat`, `.txt`, `.png` |
| `tests/test_crte_fixed_surrogate_varx.m` | focused gates (signature, dual, stability, selection) |
| `centered_smpc_step.m` / `smooth_noise_profile.m` | local helpers copied from copyAM |
| `results/copyAN_crte_fixed_surrogate_metrics.txt` | last run metrics (gitignored) |
| `results/copyAN_crte_fixed_surrogate_data.mat` | last run raw arrays (gitignored) |
| `results/copyAN_crte_fixed_surrogate_fig.png` | last run 5-row figure (gitignored) |

## Run

```matlab
cd('.../experiments/copyAN_crte_fixed_surrogate')
test_crte_fixed_surrogate_varx
copyAN_crte_fixed_surrogate
```

## Focused test result (seed 2107, p=9, q=2, ell=5)

```
PASS CRTE fixed surrogate: dual=2.20e-15 spectral=0.9092
     selected=[mu=0.25 alpha=0.00 beta=0.00]
     candidates=45 valid=27 val NRMSE=0.0863
```

## Closed-loop result (seed 20260710, p=30, q=2, ell=5)

```
selected mu=0.250 alpha=0.500 beta=0.500 validation NRMSE=0.0647
num candidates=45 valid after gates=17
tracked right=0  left=0  columns=0  dual=3.45e-15
4-piece dual completion errors max=3.45e-15
spectral radius max|eig(Ahat)|=0.9417
MAE=[0.0606 0.0574]   RMSE=[0.0834 0.0796]
upper violation=[0 0.0008]   abs violation=[0 0.0008]
QP success=1.0000   fallback=0   max QP residual=5.80e-12
constraint active rate=0.6638   cover=1.0000
avg J=51.67   u sat < 0.2%
```

Same plant seed `20260710` allows direct comparison with `copyAA` and
`copyAM`. Baseline `main/`, `copyX`, `copyAA`, `copyAM` are **not** modified.