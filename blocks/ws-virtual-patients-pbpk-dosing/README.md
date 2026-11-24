# Virtual patients (PBPK dosing)

A whole-body physiologically based pharmacokinetic model of three drugs in thousands of virtual patients: trials, the best dose per population, and a live bedside with vitals and a 3D patient

Drugs: theophylline (liver, smoking), midazolam (sedative, stored in fat) and digoxin (kidneys, heart). Run a virtual trial: a population on a regimen; who ends up in the target range, below it, above it, and who reaches toxic levels, overall and by group. Find the dose: for every population, the dose that puts the most patients in range while keeping toxic levels rare. Bedside (live): one patient at a time, chosen and dosed from the Bedside monitor view (patient, drug, dose, extra dose, stop, speed); the organs coloured by how much drug they hold, the vitals responding to the drug, a healthy adult on the same drug for comparison. Make the 3D models: the patient (from a labelled CT case, or a mockup), the platform and the scene, also streamed for the Model files view. Concentrations are in each drug's unit (mg/L for theophylline, ng/mL for the others).

An [EsaSim](https://esacode.com/) blueprint block: it runs on a worker, in Python (`backend.py`, `def run(ctx)`).

| | |
|---|---|
| Block | `ws-virtual-patients-pbpk-dosing` · version 19 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |
| Runs on | Any worker, NVIDIA GPU |

## Inputs

| Port | Type | Description |
|---|---|---|
| `ct` CT case | Files (gz, nii) | Make the 3D models: a CT case with its labels (ct.nii.gz and segmentations/<structure>.nii.gz, as the TotalSegmentator dataset gives them). Without it, the CT case setting, else a mockup body. Several connections. Optional. |

## Outputs

| Port | Type | Description |
|---|---|---|
| `course` Course | Table (csv) | course.csv: the plasma concentration over time, percentiles over the patients (p05, p25, median, p75, p95), with the target range and the toxic level |
| `patients` Patients | Table (csv) | patients.csv: one row per virtual patient (who they are, clearance, half-life, trough, peak, average, outcome) |
| `groups` By group | Table (csv) | groups.csv: in range, below, above and toxic peaks for each group (smokers, age bands, sex, conditions) |
| `doses` Doses tried | Table (csv) | Find the dose: dose_finding.csv, every dose tried for every population |
| `best` Best doses | Table (csv) | Find the dose: best_doses.csv, the dose chosen for each population |
| `models` 3D models | Files (glb, json) | Make the 3D models: patient.glb, platform.glb, ring.glb, patient.scene.json for the project's assets (models/) |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `mode` What to do | select | `"Run a virtual trial"` |  |
| `drug` Drug | select | `"Theophylline"` |  |
| `population` Population | select | `"Healthy adults"` |  |
| `patients` Virtual patients | number | `2000` |  |
| `female_pct` Women (%) | number | `50` |  |
| `route` Route | select | `"Oral (immediate release)"` |  |
| `dose_mg` Dose | number | `300` | mg (or mg/kg): theophylline 300 mg, midazolam a few mg, digoxin 0.125-0.25 mg |
| `dose_basis` Dose in | select | `"mg"` |  |
| `interval_h` Every (h) | number | `12` |  |
| `loading_mg` Extra with the first dose | number | `0` | A loading dose, in the same unit as the dose |
| `infusion_h` IV infusion over (h) | number | `0.5` |  |
| `days` Days | number | `5` |  |
| `target` Target range (theophylline) | select | `"10-20 mg/L (classic)"` | Midazolam: sedation 70-260 ng/mL; digoxin: 0.5-2 ng/mL |
| `population_scope` Find the dose for | select | `"All populations"` |  |
| `max_toxic_pct` Find the dose: toxic levels at most (%) | number | `5` |  |
| `bedside_patient` Bedside: first patient | select | `"Elderly woman, 84"` | The Bedside monitor view changes the patient and the drug while it runs |
| `bedside_minutes` Bedside: minutes it stays open | number | `20` |  |
| `ct_case` Make the 3D models: CT case | text | `""` | With nothing wired to CT case: a case of the TotalSegmentator dataset (such as s1400), read in place from Zenodo. Empty: the mockup body. |
| `accelerator` Run on | select | `"auto"` |  |
| `seed` Seed | number | `1` |  |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
