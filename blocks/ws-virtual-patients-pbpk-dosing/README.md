# Virtual patients (PBPK dosing)

A whole-body physiologically based pharmacokinetic model of twelve drugs in virtual patients: trials, the best dose per population, and a live bedside where drugs are combined on one patient

Twelve drugs: theophylline, midazolam, digoxin, morphine, naloxone, flumazenil, propranolol, atropine, ondansetron, paracetamol, amlodipine, salbutamol. Run a virtual trial: a population on a regimen; who ends up in the target range, below it, above it, and who reaches toxic levels, overall and by group. Find the dose: for every population, the dose that puts the most patients in range while keeping toxic levels rare. Bedside (live): one patient (fourteen profiles) dosed from the 3D view's controls; any number of drugs at once, each with its own whole-body model, their effects combining on the vitals (antidotes reverse, beta-blockers blunt, side effects add up), the organs coloured by how much of the chosen drug they hold, a healthy adult on the same doses for comparison. Make the 3D models: the patient (from a labelled CT case, or a mockup) and the scene. Concentrations are in each drug's unit (mg/L for theophylline and paracetamol, ng/mL for the others).

An [EsaSim](https://esacode.com/) blueprint block: it runs on a worker, in Python (`backend.py`, `def run(ctx)`).

| | |
|---|---|
| Block | `ws-virtual-patients-pbpk-dosing` · version 41 |
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
| `models` 3D models | Files (glb, json) | Make the 3D models: patient.glb and patient.scene.json, for the project's assets |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `mode` What to do | select | `"Run a virtual trial"` |  |
| `drug` Drug | select | `"Theophylline"` |  |
| `population` Population | select | `"Healthy adults"` |  |
| `patients` Virtual patients | number | `2000` |  |
| `female_pct` Women (%) | number | `50` |  |
| `route` Route | select | `"Oral (immediate release)"` |  |
| `dose_mg` Dose | number | `300` | mg (or mg/kg) per dose: theophylline 300, midazolam 2-5, digoxin 0.125-0.25, morphine 4-10, propranolol 40, paracetamol 1000, amlodipine 5-10... |
| `dose_basis` Dose in | select | `"mg"` |  |
| `interval_h` Every (h) | number | `12` |  |
| `loading_mg` Extra with the first dose | number | `0` | A loading dose, in the same unit as the dose |
| `infusion_h` IV infusion over (h) | number | `0.5` |  |
| `days` Days | number | `5` |  |
| `target` Target range (theophylline) | select | `"10-20 mg/L (classic)"` | Theophylline only; the other drugs use their own range (midazolam sedation 70-260 ng/mL, digoxin 0.5-2 ng/mL, morphine 20-80 ng/mL...) |
| `population_scope` Find the dose for | select | `"All populations"` |  |
| `max_toxic_pct` Find the dose: toxic levels at most (%) | number | `5` |  |
| `bedside_patient` Bedside: first patient | select | `"Healthy man 35"` | The 3D view's Patient control changes it while it runs |
| `bedside_minutes` Bedside: minutes it stays open | number | `20` |  |
| `ct_case` Make the 3D models: CT case | text | `""` | With nothing wired to CT case: a case of the TotalSegmentator dataset (such as s1400), read in place from Zenodo. Empty: the mockup body. |
| `accelerator` Run on | select | `"auto"` |  |
| `seed` Seed | number | `1` |  |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
