<!-- esasim: cover -->
![Virtual patients: drug dosing study](assets/cover.svg)
<!-- /esasim: cover -->

<div align="center">

<!-- Banner: upload a picture named banner.png under Data › Readme pictures; it shows here. -->

# Virtual patients: drug dosing study

**A whole-body physiologically based pharmacokinetic (PBPK) model in thousands of virtual patients: who a dosing regimen puts in the target range, who it under- or over-exposes, the best dose per population, and one patient followed live in 3D.**

![type](https://img.shields.io/badge/type-Research-blue) ![language](https://img.shields.io/badge/language-Python-3776AB) ![runs on](https://img.shields.io/badge/runs%20on-CPU%20%7C%20NVIDIA-6f42c1) [![made with](https://img.shields.io/badge/made%20with-EsaSim-0b7285)](https://esacode.com/)

<!-- More badges: ![licence](https://img.shields.io/badge/licence-MIT-green) — label, value and colour, written in the address. -->

</div>

## Flows

- **1 · Virtual trial**: A population on a dosing regimen: the plasma concentration over time (5th to 95th percentile), who ends up in the target range, below it, above it, and who reaches toxic peaks, overall and by group.
- **Maintenance · Make the 3D models**: The patient (a mockup until the anatomy models are in), the platform, the heartbeat and the scene, into the project's assets (models/).

## What you need

- A computer with the EsaSim worker (Project > Compute); it runs on any worker, faster with: NVIDIA GPU.
- Python packages, installed by the worker the first time: `numpy`, `torch`, `scikit-image`.

## How to run it

1. Under **Compute**, check that a worker is available for this project.
2. Press **Run**, choose the worker, and follow the run on the flow's **View** page.
3. The files it made are under **Results**; they stay on the worker until you save or download them.

<!-- esasim: files of the project -->

---

An [EsaSim](https://esacode.com/) project of the workspace Gianluca De Novi's workspace (Research).

### Blueprints in this repository

- [1 · Virtual trial](blueprints/1-virtual-trial.json) · 0 blocks
- [Maintenance · Make the 3D models](blueprints/maintenance-make-the-3d-models.json) · 1 block · version 1

### Workspace blocks it uses

- [Virtual patients (PBPK dosing)](blocks/ws-virtual-patients-pbpk-dosing) · version 3

Keys pasted in blocks (secret settings) are not saved here. Datasets, results and the team stay in EsaSim.
