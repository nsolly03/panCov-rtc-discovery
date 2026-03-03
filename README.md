# RTC Pan-Coronavirus Inhibitor Discovery

**Candidate:** Olivier Nsekuye  
**Institution:** University of Liège — GIGA-VIN Lab  
**Supervisor:** Prof. Jean-Claude Twizere  
**Funding:** FRIA-B1 Fellowship  
**Timeline:** 2025-2029  

## Objective
Identify pan-coronavirus small-molecule inhibitors targeting protein-protein
interfaces within the SARS-CoV-2 Replication-Transcription Complex (RTC).

## 8 Binary Interface Targets
| # | Complex | Priority |
|---|---------|----------|
| 1 | NSP10-NSP16 | Critical |
| 2 | NSP12-NSP7 | Critical |
| 3 | NSP12-NSP8 | Critical |
| 4 | NSP9-NSP12 | Critical |
| 5 | NSP10-NSP14 | High |
| 6 | NSP13-Helicase | High |
| 7 | NSP12-NSP13 | Medium |
| 8 | NSP7-NSP8 | Low-Medium |

## How to reproduce this project
Follow WORKLOG.md step by step from Entry 001.

## Project structure
```
00-reference/        — PDB structures + sequences
01-alphafold3/       — AF3 predictions per complex
02-validation/       — Interface validation results
03-virtual-screening/— Docking inputs and results
04-hits/             — Validated hit compounds
05-experimental/     — Wet lab protocols and data
scripts/             — All Python scripts
notebooks/           — Jupyter analysis notebooks
results/             — Figures and summary tables
docs/                — Literature and notes
```
