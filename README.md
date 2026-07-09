# PredVAR + MPC MATLAB Workspace

This repository stores MATLAB code, notes, and experiment materials for the PredVAR + MPC research workflow.

## Purpose

- Keep MATLAB scripts and functions under version control.
- Record important experiment changes through Git commits.
- Avoid losing working versions when testing new ideas.

## Suggested structure

```text
.
├── README.md          # Project entry note
├── .gitignore         # Files Git should ignore
├── src/               # MATLAB source code
├── experiments/       # Experiment scripts
├── results/           # Small result summaries or figures, if tracked
└── data/              # Data files, usually not tracked if large
```

## Basic Git workflow

```bash
git status
git add .
git commit -m "Describe the change"
```

If this repository is later connected to GitHub:

```bash
git push
```

## Notes

Large generated files such as `.mat` datasets, Simulink build folders, and temporary MATLAB autosave files are ignored by default. If a data file is small and important, add it manually with `git add -f filename.mat`.
