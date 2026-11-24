# Virtual patients (PBPK dosing)

A whole-body physiologically based pharmacokinetic model of a drug in thousands of virtual patients: a trial of a dosing regimen, the best dose per population, and one patient followed live in 3D

Run a virtual trial: a population (healthy adults, smokers, elderly, liver cirrhosis, heart failure, kidney impairment, a mixed clinical population) on a regimen; who ends up in the target range, below it, above it, and who reaches toxic peaks, overall and by group. Find the dose: for every population, the dose that puts the most patients in range while keeping toxic peaks rare. Follow one patient (3D): the typical, most exposed, least exposed or a random patient, each organ coloured by the drug's concentration there, live. Make the 3D models: the patient, the platform and the scene for the project's assets. The model is linear, so each run is a few seconds even for 10,000 patients (GPU or CPU).

An [EsaSim](https://esacode.com/) blueprint block: it runs on a worker, in Python (`backend.py`, `def run(ctx)`).

| | |
|---|---|
| Block | `ws-virtual-patients-pbpk-dosing` · version 7 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |
| Runs on | Any worker, NVIDIA GPU |

## Outputs

| Port | Type | Description |
|---|---|---|
| `course` Course | Table (csv) | course.csv: the plasma concentration over time, percentiles over the patients (p05, p25, median, p75, p95), with the target range and the toxic level |
| `patients` Patients | Table (csv) | patients.csv: one row per virtual patient (who they are, clearance, half-life, trough, peak, average, outcome) |
| `groups` By group | Table (csv) | groups.csv: in range, below, above and toxic peaks for each group (smokers, age bands, sex, conditions) |
| `doses` Doses tried | Table (csv) | Find the dose: dose_finding.csv, every dose tried for every population |
| `best` Best doses | Table (csv) | Find the dose: best_doses.csv, the dose chosen for each population |
| `models` 3D models | Files (glb, json, wav) | Make the 3D models: patient.glb, platform.glb, ring.glb, heartbeat.wav, patient.scene.json for the project's assets (models/) |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `mode` What to do | select | `"Run a virtual trial"` |  |
| `drug` Drug | select | `"Theophylline"` |  |
| `population` Population | select | `"Healthy adults"` |  |
| `patients` Virtual patients | number | `2000` |  |
| `female_pct` Women (%) | number | `50` |  |
| `route` Route | select | `"Oral (immediate release)"` |  |
| `dose_mg` Dose | number | `300` |  |
| `dose_basis` Dose in | select | `"mg"` |  |
| `interval_h` Every (h) | number | `12` |  |
| `loading_mg` Extra with the first dose | number | `0` | A loading dose, in the same unit as the dose |
| `infusion_h` IV infusion over (h) | number | `0.5` |  |
| `days` Days | number | `5` |  |
| `target` Target range | select | `"10-20 mg/L (classic)"` |  |
| `population_scope` Find the dose for | select | `"All populations"` |  |
| `max_toxic_pct` Find the dose: toxic peaks at most (%) | number | `5` |  |
| `follow` Follow | select | `"The typical patient"` |  |
| `speed` Follow: simulated hours per second | number | `1` |  |
| `accelerator` Run on | select | `"auto"` |  |
| `seed` Seed | number | `1` |  |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
