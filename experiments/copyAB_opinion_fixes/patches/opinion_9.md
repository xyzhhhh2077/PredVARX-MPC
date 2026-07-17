# Opinion 9 patch report

status: done

## Request

增加可选终端代价开关（默认关）并文档化未证明稳定。

## Design

- Object: optional quadratic terminal weight on final latent state \(z_N\).
- Flag: `opt.use_terminal_cost` (default **false** / missing → original behaviour).
- Construction when enabled:
  1. Build stage weight \(Q_f\) from `opt.Qf` if provided, else \(P^\top Q P\) (fallback \(I\)).
  2. Solve discrete Lyapunov \(A^\top P_{\mathrm{term}} A - P_{\mathrm{term}} + Q_f = 0\) via `dlyap(A', Qf)`.
  3. Predict \(z_N = A^N z + G_z(U-U_0)\) and add \(z_N^\top P_{\mathrm{term}} z_N\) into the raw QP objective.
  4. On any failure (unstable \(A\), missing toolbox, non-PSD \(P_{\mathrm{term}}\)): **skip** and keep stage-only cost.
- Documentation: README states **严格凸 ≠ 递归可行 ≠ 稳定**; terminal cost is not a terminal-set certificate.

## Files

| File | Change |
|---|---|
| `centered_smpc_step.m` | optional terminal cost block + `out.use_terminal_cost` / `out.terminal_cost_applied` / `out.Pterm` |
| `README.md` | stability / recursive-feasibility caveat under hard boundaries |
| `tests/test_opinion09_terminal_default_off.m` | default-off shape regression + ON path callable |
| `patches/opinion_9.md` | this report |

## Baseline preservation

- Default path (flag missing or `false`) does not touch Hessian construction beyond the pre-existing stage terms + `1e-9` regularizer.
- `main/`, `copyAA`, `copyX` untouched.

## Verification

Run:

```matlab
cd('.../experiments/copyAB_opinion_fixes');
addpath(pwd); addpath(fullfile(pwd,'tests'));
test_opinion09_terminal_default_off
```

## Limitations (honest)

- Lyapunov terminal weight is a soft regularizer only.
- No terminal invariant set, no dual-mode / tube certificate, no recursive-feasibility proof.
- If `dlyap` fails the flag is effectively a no-op (`terminal_cost_applied=false`).
