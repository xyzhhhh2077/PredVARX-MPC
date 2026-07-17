# Opinion 3 patch report

status: done

## Goal

Clarify that in the free-IVR loop:

- `M` is the **regularized prediction association matrix** (ridge objective-drop Gram), used only to re-estimate free loading directions via `eig(M)`.
- `trace_hist(iter) = τ_i` is **free-latent prediction energy** `tr(Xpred Xpred')/T`, **not** output predictable energy and **not** `tr(M)`.

Do **not** call `M` a “predicted output energy matrix.”

## Math (from pure-math review)

Ridge map:

```text
L_* = Y_+ Xhat_+' (S + λ I)^{-1},   S = Xhat_+ Xhat_+'
M   = Y_+ Xhat_+' (S + λ I)^{-1} Xhat_+ Y_+'
Yhat = L_* Xhat_+
Yhat Yhat' = Y_+ Xhat_+' (S+λI)^{-1} S (S+λI)^{-1} Xhat_+ Y_+'
```

Identity:

```text
M - Yhat Yhat' = λ Y_+ Xhat_+' (S+λI)^{-2} Xhat_+ Y_+'  ⪰ 0
```

- `λ = 0` and `S` invertible (or consistent Moore–Penrose): `M = Yhat Yhat'`.
- Code uses `λ = 1e-8`: numerically close when `S` is well-conditioned.
- 1-D counterexample `Y X' = S = 1, λ = 1`: `M = 1/2`, `Yhat Yhat' = 1/4`.

`τ_i` depends only on `Xpred` (latent coords `z = V'y`), not on free-ambient `Ycur`.

## Code change

File: `split_control_free_ivr_varx.m` — **comments only** around the free-IVR loop:

- Block comment before the loop defining `Xpred`, `τ_i = trace_hist`, and `M`.
- Inline notes on the `trace_hist` and `M` assignment lines.
- No algorithm / formula change.

## Test

`tests/test_opinion03_M_vs_tau.m` (standalone numeric; does not require a full identifier run):

1. Synthetic `Ycur`, `Xpred` with well-conditioned `S`.
2. `λ = 1e-8`: `rel(M, Yhat Yhat') < 1e-6` (observed `~2e-11`).
3. Residual identity `M - Yhat Yhat' = λ · …` holds.
4. `λ = 0`: exact `M = Yhat Yhat'`.
5. 1-D `λ = 1`: `M = 0.5`, `Yhat Yhat' = 0.25`.
6. `τ` independent of `Ycur`; `M` depends on `Ycur`.

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'tests'));
test_opinion03_M_vs_tau;
```

### MATLAB result (real run)

```text
PASS opinion03 M vs tau: lambda=1e-8 rel(M,YhatYhat')=1.977e-11;
lambda=0 exact; 1-D M=0.5 YhatYhat=0.25; tau=3.08425 (latent only)
```

## Boundary

- Comment/doc fix only: IVR still maximizes free loading directions via `eig(M)`; stop criterion still monitors `τ` + subspace `δ` (see Opinion 4).
- `τ` stationarity is **not** a proof of monotone convergence of the ridge objective.
