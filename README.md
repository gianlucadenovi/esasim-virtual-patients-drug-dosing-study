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
- **2 · Find the dose**: For every population, the dose that puts the most patients in the target range while keeping toxic peaks rare; every dose tried, as curves.
- **3 · Follow one patient (3D)**: One virtual patient on the regimen, live: each organ coloured by the drug's concentration in it, the heart beating, the doses and the moments the plasma level leaves the target range.
- **Maintenance · Make the 3D models**: The patient (a mockup until the anatomy models are in), the platform, the heartbeat and the scene, into the project's assets (models/).

## What you need

- A computer with the EsaSim worker (Project > Compute); it runs on any worker, faster with: NVIDIA GPU.
- Python packages, installed by the worker the first time: `numpy`, `torch`, `scikit-image`.

## How to run it

1. Under **Compute**, check that a worker is available for this project.
2. Press **Run**, choose the worker, and follow the run on the flow's **View** page.
3. The files it made are under **Results**; they stay on the worker until you save or download them.

## What you get

- On the **View** page: *Summary* (1 · Virtual trial), *Outcome (%, mg/L, h)* (1 · Virtual trial), *Concentration over time* (1 · Virtual trial), *By group* (1 · Virtual trial), *Spread* (1 · Virtual trial), *Patients* (1 · Virtual trial), *Groups* (1 · Virtual trial), *The dose for each population* (2 · Find the dose), *Best dose* (2 · Find the dose), *In range by dose* (2 · Find the dose), *Toxic peaks by dose* (2 · Find the dose), *The patient, live* (3 · Follow one patient (3D)).

<!-- esasim: files of the project -->

---

An [EsaSim](https://esacode.com/) project of the workspace Gianluca De Novi's workspace (Research).

### Blueprints in this repository

- [1 · Virtual trial](blueprints/1-virtual-trial.json) · 8 blocks
- [2 · Find the dose](blueprints/2-find-the-dose.json) · 5 blocks
- [3 · Follow one patient (3D)](blueprints/3-follow-one-patient-3d.json) · 4 blocks
- [Maintenance · Make the 3D models](blueprints/maintenance-make-the-3d-models.json) · 1 block · version 1

### Workspace blocks it uses

- [Virtual patients (PBPK dosing)](blocks/ws-virtual-patients-pbpk-dosing) · version 4

Keys pasted in blocks (secret settings) are not saved here. Datasets, results and the team stay in EsaSim.
