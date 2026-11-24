<!-- esasim: cover -->
![Virtual patients: drug dosing study](assets/cover.jpg)
<!-- /esasim: cover -->

<div align="center">

<!-- Banner: upload a picture named banner.png under Data › Readme pictures (the bedside with the 3D patient and the monitor works well), then: ![Virtual patients: drug dosing study](pictures/banner.png) -->

# Virtual patients: drug dosing study

**Twelve drugs in virtual patients: who a dose works for, the best dose for each population, and a live bedside where you give drugs, watch their side effects and correct them with other drugs.**

A whole-body physiologically based pharmacokinetic (PBPK) model follows each drug through the organs of every virtual patient, from their age, size, sex, smoking, and liver, heart and kidney function. Run thousands of patients at once to see who a regimen puts in range, or sit at one patient's bedside with a 3D body built from a real CT scan, a bedside monitor and the drugs' levels.

![type](https://img.shields.io/badge/type-Medical%20%26%20healthcare-blue) ![language](https://img.shields.io/badge/language-Python-3776AB) ![runs on](https://img.shields.io/badge/runs%20on-CPU%20%7C%20NVIDIA-6f42c1) ![method](https://img.shields.io/badge/method-PBPK%20%2B%20PD-orange) ![data](https://img.shields.io/badge/data-CC%20BY%204.0-green) [![made with](https://img.shields.io/badge/made%20with-EsaSim-0b7285)](https://esacode.com/)

</div>

> **For research and teaching only.** Not a validated clinical tool: do not use it to dose real patients.

## Overview

The same dose behaves very differently in different people. A smoker clears theophylline 1.6 times faster than a non-smoker; an 84-year-old with poor kidneys builds up digoxin to toxic levels on a dose that suits a 35-year-old; midazolam, stored in fat, sedates an older patient for much longer. And drugs change each other's effects: an antidote reverses a sedative, a beta-blocker slows the heart that a bronchodilator speeds up but narrows the airways it opens.

This project models both. Each drug moves through 14 organs and the blood in each patient's own body (the pharmacokinetics, mechanistic); its levels at the sites where it acts set the vital signs, with the drugs adding up or cancelling (the pharmacodynamics, simplified). The model is linear in each drug, so it is solved exactly, for one patient live or for thousands at once.

## Flows

### 1 · Virtual trial

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/104)

A population on a regimen (drug, dose, interval, route, loading dose, days): the plasma level over time (5th to 95th percentile), the spread of steady-state levels, who ends up in range, below, above or toxic, by group, and a table of every patient.

### 2 · Find the dose

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/105)

Every population across the drug's dose range: a curve per population of the share in range and the share at toxic levels, and the best dose for each.

### 3 · Bedside: one patient, many drugs, live (3D)

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/106)

| Block | What it does |
|---|---|
| **The bedside** | The patient and every drug given, simulated live, with a healthy adult on the same doses for comparison. |
| **The patient, live** | Sends the patient, the vitals and the events to the views. |
| **The patient (3D)** | The patient built from a CT scan, the organs coloured by the chosen drug, the heart and lungs moving and heard; the controls. |
| **Bedside monitor** | ECG, pulse oximetry and breathing waveforms, the vital signs and the alarms. |
| **Drug levels** | Each drug given, as a share of the top of its target range, over time. |

### Maintenance · Make the 3D models

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/103)

Builds `patient.glb` from a CT case and shows it in a window with a download button, to put in the project's assets.

## The drugs

| Drug | What it does | Side effects in the model | Corrected by |
|---|---|---|---|
| **Theophylline** | opens the airways | fast heart, ectopic beats, nausea, tremor; seizure risk above 40 mg/L | propranolol (but it narrows asthmatic airways) |
| **Midazolam** | sedation | deep sedation, slow breathing, low SpO₂, low blood pressure | flumazenil |
| **Digoxin** | heart failure, atrial fibrillation | slow heart, heart block, ventricular ectopy, nausea, longer QTc | atropine |
| **Morphine** | pain relief | sedation, slow breathing, low SpO₂, nausea | naloxone; ondansetron for the nausea |
| **Naloxone** | reverses opioids | wears off before morphine; pain comes back | – |
| **Flumazenil** | reverses benzodiazepines | wears off before midazolam | – |
| **Propranolol** | slows the heart, lowers blood pressure | narrows the airways in asthma, blunts salbutamol | – |
| **Atropine** | speeds the heart | delirium, hot dry skin at high levels | – |
| **Ondansetron** | stops nausea | longer QTc | – |
| **Paracetamol** | mild pain relief, lowers fever | liver injury risk at high levels | – |
| **Amlodipine** | lowers blood pressure slowly, for days | a reflex rise in heart rate | propranolol |
| **Salbutamol** | opens the airways | fast heart, tremor | – |

Theophylline and paracetamol are in mg/L, the others in ng/mL.

## The patients

Healthy man 35, healthy woman 30, athlete 24 (resting HR 50), smoker 45, obese woman 40 (BMI 38), asthma attack 28, after surgery 50 (pain 7/10), pneumonia with fever 64, elderly woman 84, frail man 91, liver cirrhosis 58, heart failure 72, kidney failure 66 (GFR 20), and a random patient. The trials draw thousands from seven populations: healthy adults, smokers, elderly (65+), liver cirrhosis, heart failure, kidney impairment, and a mixed clinical population.

## What you need

- A machine with the EsaSim worker (Project › Compute). Any machine: an NVIDIA GPU speeds up the large trials but is not needed.
- Python packages, installed by the worker the first time: `numpy`, `torch`, `scipy`, and for the 3D models `SimpleITK`, `scikit-image`, `trimesh`, `fast-simplification`, `remotezip`, `requests`.
- Internet access on the machine only to make the 3D models (the CT case is read from Zenodo).
- No account or key.

## How to run it

1. Under **Compute**, check that a machine is available for this project.
2. Open a flow under **Blueprint**, set the drug, population and regimen on its block, and press **Run**.
3. Open the **View**.

At the bedside, choose the patient and the drug with the controls at the bottom of the 3D view, then **Give dose** or **Regular doses**; **Dose ×** sets a multiple of the usual dose, **Clock** how fast the patient's time runs (15 min to 12 h a second). Click **Sound on** in the 3D view for the heartbeat and the breathing. A few things to try:

- *Elderly woman 84*, morphine three times: breathing slows and SpO₂ falls, with alarms. Then naloxone: she breathes again, until it wears off.
- *Kidney failure 66*, digoxin with regular doses at 4 h a second: over a few days the heart slows and heart block appears. Atropine lifts it for a while.
- *Asthma attack 28*: salbutamol opens the airways and SpO₂ rises; propranolol narrows them again.
- *After surgery 50*: morphine treats the pain but brings nausea; ondansetron stops it and lengthens the QTc.
- *Pneumonia with fever 64*: paracetamol brings the temperature down.

## What you get

- **Virtual trial**: `course.csv` (the plasma level over time, by percentile), `patients.csv` (every patient: who they are, clearance, half-life, trough, peak, average, outcome), `groups.csv` (in range, below, above and toxic, by group), with their charts and tables in the View.
- **Find the dose**: `dose_finding.csv` (every dose tried for every population) and `best_doses.csv`.
- **Bedside**: in the 3D view, each organ coloured by how much of the chosen drug it holds (log scale, indigo to pale yellow, with the legend), the gut glowing while it absorbs, the liver and kidneys while they clear; the heart beating at the live rate (dropped beats in heart block included), the lungs and ribcage moving with each breath; the drugs on board and the events (doses, levels leaving the range, rhythm, sedation, pain and nausea changes, alarms). The monitor shows HR and rhythm, SpO₂, breathing rate, blood pressure, temperature, QTc, pain, sedation (Ramsay), nausea, airways and the alarms.

## What is computed

| Part | How |
|---|---|
| **Organs** | lungs, brain, heart, kidneys, liver, gut, spleen, pancreas, fat, bone, muscle, skin, the rest of the body, arterial and venous blood, the gut contents; the gut, spleen and pancreas drain into the liver |
| **Each body** | organ volumes and blood flows from body surface area; fat (Deurenberg), blood volume (Nadler), muscle the weight left over; heart failure lowers cardiac output |
| **Distribution** | blood-flow limited; tissue:plasma partition from tissue composition (Rodgers & Rowland), scaled to the observed volume of distribution |
| **Elimination** | well-stirred liver on the unbound drug, calibrated to the observed clearance, plus renal excretion scaled by GFR; smoking, cirrhosis and heart failure change the liver's capacity, drug by drug |
| **Solution** | exact 6-minute steps with a matrix exponential (PyTorch, CPU or GPU) |
| **Effects** | receptor occupancy (Emax) at each drug's effect site; competitive antagonists (naloxone, flumazenil); beta-blockade; breathing depression by opioids and benzodiazepines together; airway narrowing; heart rhythm (ectopy, heart block, QTc); pain, fever and nausea |

Checked against the literature for a 73 kg adult: theophylline half-life about 9 h (observed about 8 h), midazolam about 2.7 h (1.8–6.4 h), digoxin about 40–50 h (about 40 h); digoxin 0.25 mg a day gives about 1 ng/mL in a healthy adult and over 2 ng/mL with a GFR of 20–40 mL/min.

## Data, method and licences

- **3D patient**: case s1400 of the **TotalSegmentator dataset v2.0.1, small subset**, [Zenodo record 10047263](https://zenodo.org/records/10047263), licence **CC BY 4.0**. Wasserthal J., Breit H.-C., Meyer M.T., et al. *TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images.* Radiology: Artificial Intelligence 2023; 5(5): e230024. Each labelled structure goes into the organ of the model it belongs to.
- **Sounds**: "Heartbeat (single)" by Universfield and "Deep breath" by hasin2004, from Pixabay.
- **Physiology**: ICRP reference adult organ volumes and blood flows, as tabulated with PK-Sim's values in *Modeling principles for a physiology-based whole-body model of human metabolism* (arXiv:2404.14224). Tissue composition after Ruark et al. (Metrum Research Group, ACoP 2019). Rodgers T, Rowland M. J Pharm Sci 2006;95:1238–57. Deurenberg P et al. (1991); Nadler SB et al. (1962); Mosteller (1987).
- **Theophylline**: Clin Pharmacokinet 1978;3:267–93; Abduljalil K et al. Front Pediatr 2022; Jusko WJ et al. J Pharm Sci 1979;68:1358–66.
- **Midazolam**: Molecules 2020;25:1080; the US prescribing information; Sunzel M et al. Br J Clin Pharmacol 1988.
- **Digoxin**: Iisalo E. Clin Pharmacokinet 1977;2:1–16.
- **The other nine drugs**: their prescribing information and standard pharmacology references.
- The code of the blocks is EsaSim's.

## Notes

- The pharmacokinetics are mechanistic; the effects, the interactions and the baseline vital signs are simplified and illustrative, following published thresholds and half-effect levels where they exist.
- The nine added drugs use typical literature values; strong bases' binding to tissue phospholipids is not modelled, so their organ-by-organ pattern is approximate (their total volume is right). Clearances above the liver's blood flow are capped.
- Drugs interact through their effects only, not their kinetics (no enzyme inhibition or induction). Absorption is first order, elimination linear. Digoxin's binding to the sodium pump is not modelled.
- The 3D body is one scanned person, neck to pelvis: it does not change with the virtual patient, and the brain has no surface.
- The monitor's waveforms run in real time at the patient's current rates while the patient's clock runs faster.

<!-- esasim: files of the project -->

---

An [EsaSim](https://esacode.com/) project of the workspace Gianluca De Novi's workspace (Medical & healthcare).

### Blueprints in this repository

- [1 · Virtual trial](blueprints/1-virtual-trial.json) · 12 blocks · version 3
- [2 · Find the dose](blueprints/2-find-the-dose.json) · 7 blocks · version 1
- [3 · Bedside: one patient, many drugs, live (3D)](blueprints/3-bedside-one-patient-many-drugs-live-3d.json) · 5 blocks · version 6
- [Maintenance · Make the 3D models](blueprints/maintenance-make-the-3d-models.json) · 2 blocks · version 3

### Workspace blocks it uses

- [Virtual patients (PBPK dosing)](blocks/ws-virtual-patients-pbpk-dosing) · version 34
- [Bedside monitor](blocks/ws-bedside-monitor) · version 2
- [Model files (preview and download)](blocks/ws-model-files-preview-and-download) · version 1

### Assets

The files the project's views draw with (3D models, textures, sounds), in [assets/](assets/):

- [covers/fda-project-optimus-selecting-optimal-dose-oncology.webp](assets/covers/fda-project-optimus-selecting-optimal-dose-oncology.webp) · 355 KB
- [hasin2004-deep-breath-247459.mp3](assets/hasin2004-deep-breath-247459.mp3) · 148 KB
- [patient.glb](assets/patient.glb) · 12 MB
- [patient.scene.json](assets/patient.scene.json) · 6.0 KB
- [universfield-heartbeat-single-383748.mp3](assets/universfield-heartbeat-single-383748.mp3) · 35 KB

Keys pasted in blocks (secret settings) are not saved here, nor are the team and the results of runs.
