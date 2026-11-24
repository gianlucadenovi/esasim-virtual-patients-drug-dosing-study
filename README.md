<!-- esasim: cover -->
![Virtual patients: drug dosing study](assets/cover.jpg)
<!-- /esasim: cover -->

# Virtual patients: drug dosing study

A whole-body **physiologically based pharmacokinetic (PBPK) model** run in virtual patients. The same dose behaves very differently in different people, and this project shows why:

- **Who does a regimen work for?** Thousands of patients at once: who ends up in the target range, below it, above it, or toxic, overall and by group.
- **What dose suits each population?** For every population, the dose that puts the most patients in range while keeping toxic levels rare.
- **What happens inside one patient?** A live bedside. You choose the patient and the drug, dose them, and watch where the drug goes and how the vital signs respond.

> This is a teaching and exploration model, not a validated clinical tool. Do not use it to dose real patients.

## Three drugs, three stories

| Drug | Cleared by | What it shows | Unit |
|---|---|---|---|
| **Theophylline** (bronchodilator) | the liver (CYP1A2) | Smokers clear it 1.6× faster; liver cirrhosis and heart failure slow it down, and the same dose becomes toxic. Toxic levels speed up the heart and cause ectopic beats. | mg/L |
| **Midazolam** (sedative) | the liver (CYP3A4) | Very lipophilic, so it is stored in fat. Older people and those with cirrhosis or heart failure clear it more slowly, so the same infusion sedates them more deeply: breathing slows, oxygen saturation and blood pressure fall. | ng/mL |
| **Digoxin** (heart drug) | the kidneys | Takes days to build up and reaches the heart hours after the plasma. With poor kidney function (older age, kidney failure), the usual dose becomes toxic: the heart slows and heart block appears. | ng/mL |

## The flows

| Flow | What it does |
|---|---|
| **1 · Virtual trial** | A population on a regimen (drug, dose, interval, route, loading dose, days). Shows the plasma level over time (5th to 95th percentile), the spread of steady-state levels, outcomes by group, and a table of every patient. |
| **2 · Find the dose** | Every population across the drug's dose range, with a curve per population of the share in range and the share at toxic levels, and the best dose for each. |
| **3 · Bedside: one patient, live (3D)** | Pick a patient (healthy man or woman, smoker, obese, elderly, cirrhosis, heart failure, kidney failure, or a random patient) and a drug. Change the dose, give an extra dose, stop the drug, speed up the clock. You see where the drug is in the body, a bedside monitor with ECG, pulse oximetry and breathing waveforms and alarms, and the plasma level against a healthy adult on the same drug. |
| **Maintenance · Make the 3D models** | Builds the patient from a CT case (below), plus the platform and the scene. The window shows the model and lets you download the files for the project's Assets. |

### Reading the bedside

- **Organ colour** shows how much drug each organ holds, on a log scale from indigo (hardly any) through violet, magenta and orange to pale yellow (the most). The legend on the right gives the levels for the current drug.
- **Glow** means a process is happening right now: the gut while it absorbs a dose, the liver and the kidneys while they clear the drug.
- **The heart** beats, and the heartbeat sounds, at the patient's live heart rate; the lungs breathe at the live breathing rate.
- **The monitor** shows HR, SpO₂, breathing rate, blood pressure, temperature, the drug level, and sedation (Ramsay score) for midazolam. Alarms appear when a vital sign leaves its normal range. The waveforms run in real time while the simulated clock runs faster.
- **The events list** tells the story: doses, when the level enters or leaves the range, rhythm changes, sedation changes and alarms.

## What is modelled, and how

### Pharmacokinetics (where the drug goes): mechanistic

- **Compartments.** Lungs, brain, heart, kidneys, liver, gut, spleen, pancreas, fat, bone, muscle, skin and the rest of the body, plus arterial and venous blood and the gut contents. The gut, spleen and pancreas drain into the liver, so an oral dose passes the liver first.
- **Each patient's body.** Organ volumes and blood flows scale with body surface area. Fat comes from the Deurenberg formula, blood volume from Nadler's, and muscle is the weight left over. Heart failure lowers cardiac output, and with it liver blood flow.
- **Distribution.** Blood-flow limited. Tissue:plasma partition coefficients come from each tissue's composition (water, neutral lipids, phospholipids, albumin) by the Rodgers & Rowland method. They are scaled once so that a reference adult's volume of distribution matches the observed one.
- **Elimination.** Liver metabolism of the unbound drug (well-stirred liver), calibrated to each drug's observed clearance in a reference adult, plus renal excretion scaled by kidney function (GFR). Smoking, cirrhosis and age change the liver's capacity, drug by drug.
- **Solution.** The model is linear, so it is solved exactly over each 6-minute step with a matrix exponential. The bedside steps one patient (and a healthy reference) live; the trials step thousands at once.

**Checked against the literature** for a reference 73 kg adult. Theophylline: clearance 0.040 L/h/kg, Vss 0.52 L/kg, half-life about 9 h (observed about 8 h). Midazolam: half-life about 2.7 h (observed 1.8–6.4 h). Digoxin: half-life about 40–50 h (observed about 40 h). Digoxin at 0.25 mg a day reaches about 1 ng/mL in a healthy adult, and over 2 ng/mL (toxic) with a GFR of 20–40 mL/min.

### Effects on the body (pharmacodynamics): simplified

The effects follow published concentration thresholds and half-effect concentrations, but their exact shapes are simplified:

- **Theophylline:** heart rate rises slightly in the therapeutic range and steeply above 20 mg/L; ectopic beats above 30 mg/L; seizure risk and hypotension above 40 mg/L.
- **Midazolam:** sedation from the published concentrations for Ramsay scores 2–6 (about 6, 71, 171, 260 and 659 ng/mL). Blood pressure and ventilation are half-maximally affected at 50–60 ng/mL. The fall in breathing rate and SpO₂ at deep sedation is illustrative.
- **Digoxin:** driven by the heart's own drug level, so effects lag the plasma level. The heart slows across the therapeutic range; above 2 ng/mL, heart block; above 3 ng/mL, ventricular ectopy.

**Baseline vital signs** are typical values adjusted for age, smoking and heart failure (illustrative).

## The 3D patient

The patient is built from a real CT scan and its experts' labels: case s1400 of the **TotalSegmentator dataset** (v2.0.1, small subset). Every labelled structure goes into the model compartment it belongs to, and each compartment is one part of the 3D model:

- The liver and gallbladder make the liver.
- The oesophagus, stomach, duodenum, small bowel and colon make the gut.
- The aorta and arteries make the arterial blood; the venae cavae, portal vein and other veins make the venous blood.
- Every vertebra, rib and other bone makes the skeleton.
- The skin and the fat under it come from the scan's own densities.

The scan covers the torso, so the brain is listed but has no surface.

**Credit (CC BY 4.0):** TotalSegmentator dataset v2.0.1 (small subset), Zenodo, doi:10.5281/zenodo.10047263. Wasserthal J., Breit H.-C., Meyer M.T., et al. *TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images.* Radiology: Artificial Intelligence 2023; 5(5): e230024.

**Heartbeat sound:** "Heartbeat (single)" by Universfield, trimmed to one beat.

## Sources

- **Physiology:** ICRP reference adult organ volumes and blood flows, as tabulated with PK-Sim's values in *Modeling principles for a physiology-based whole-body model of human metabolism* (arXiv:2404.14224, Table 1).
- **Tissue composition:** the unified tissue composition of Ruark et al. (*Predicting passive and active tissue:plasma partition coefficients: interindividual and interspecies variability*), as presented by Metrum Research Group (ACoP 2019, Table 1).
- **Partitioning:** Rodgers T, Rowland M. *Physiologically based pharmacokinetic modelling 2: predicting the tissue distribution of acids, very weak bases, neutrals and zwitterions.* J Pharm Sci 2006;95:1238–57.
- **Theophylline:** *Clinical pharmacokinetics of theophylline.* Clin Pharmacokinet 1978;3:267–93. Abduljalil K, Gardner I, Jamei M. *Application of a physiologically based pharmacokinetic approach to predict theophylline pharmacokinetics…* Front Pediatr 2022. Jusko WJ et al. *Factors affecting theophylline clearances…* J Pharm Sci 1979;68:1358–66. Physicochemical data from the PKKB database.
- **Midazolam:** logP and pKa from *Molecules* 2020;25:1080 (Table 1); clearance, volume and half-life ranges from the US prescribing information. Sedation concentrations from a population pharmacodynamic study of midazolam given by target-controlled infusion to surgical ICU patients. Sunzel M et al. *Respiratory and cardiovascular effects in relation to plasma levels of midazolam and diazepam.* Br J Clin Pharmacol 1988.
- **Digoxin:** logP from PKKB and QsarDB. Iisalo E. *Clinical pharmacokinetics of digoxin.* Clin Pharmacokinet 1977;2:1–16 (bioavailability, binding, half-life, renal share, tissue stores).
- **Body composition:** Deurenberg P et al. (1991) body fat; Nadler SB et al. (1962) blood volume; Mosteller (1987) body surface area.

## Limits

- Midazolam's blood:plasma ratio (0.66) and the fraction escaping its gut-wall metabolism (0.75) are typical model values, not measured here.
- Digoxin's binding to the sodium pump in muscle and heart is not modelled. Its tissue levels come from lipids and water, scaled to the observed volume, so fat holds more digoxin than it really does.
- Absorption is first order, elimination linear, and there are no drug–drug interactions.
- The effect curves are simplified, and baseline vital signs are typical values, not measurements.
- The 3D anatomy is one scanned person; it does not change with the virtual patient's sex, size or conditions.

<!-- esasim: files of the project -->

---

An [EsaSim](https://esacode.com/) project of the workspace Gianluca De Novi's workspace (Research).

### Blueprints in this repository

- [1 · Virtual trial](blueprints/1-virtual-trial.json) · 12 blocks · version 3
- [2 · Find the dose](blueprints/2-find-the-dose.json) · 7 blocks · version 1
- [3 · Bedside: one patient, many drugs, live (3D)](blueprints/3-bedside-one-patient-many-drugs-live-3d.json) · 4 blocks · version 4
- [Maintenance · Make the 3D models](blueprints/maintenance-make-the-3d-models.json) · 2 blocks · version 3

### Workspace blocks it uses

- [Virtual patients (PBPK dosing)](blocks/ws-virtual-patients-pbpk-dosing) · version 25
- [Model files (preview and download)](blocks/ws-model-files-preview-and-download) · version 1

### Assets

The files the project's views draw with (3D models, textures, sounds), in [assets/](assets/):

- [covers/fda-project-optimus-selecting-optimal-dose-oncology.webp](assets/covers/fda-project-optimus-selecting-optimal-dose-oncology.webp) · 355 KB
- [hasin2004-deep-breath-247459.mp3](assets/hasin2004-deep-breath-247459.mp3) · 148 KB
- [heartbeat.mp3](assets/heartbeat.mp3) · 2.4 KB
- [patient.glb](assets/patient.glb) · 12 MB
- [patient.scene.json](assets/patient.scene.json) · 6.3 KB
- [universfield-heartbeat-single-383748.mp3](assets/universfield-heartbeat-single-383748.mp3) · 35 KB

Keys pasted in blocks (secret settings) are not saved here, nor are the team and the results of runs.
