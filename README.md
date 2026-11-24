<!-- esasim: cover -->
![Virtual patients: drug dosing study](assets/cover.jpg)
<!-- /esasim: cover -->

# Virtual patients: drug dosing study

A whole-body **physiologically based pharmacokinetic (PBPK) model** of twelve drugs, run in virtual patients. The same dose behaves very differently in different people, and drugs change each other's effects. This project shows both:

- **Who does a regimen work for?** Thousands of patients at once: who ends up in the target range, below it, above it, or toxic, overall and by group.
- **What dose suits each population?** For every population, the dose that puts the most patients in range while keeping toxic levels rare.
- **What happens inside one patient, and how do you fix it?** A live bedside. Give a drug, watch its side effects appear on the monitor, then correct them with other drugs: an antidote, a drug that blocks the side effect, or one that treats the symptom.

> This is a teaching and exploration model, not a validated clinical tool. Do not use it to dose real patients.

## Twelve drugs

| Drug | What it does | Side effects in the model | Corrected by | Unit |
|---|---|---|---|---|
| **Theophylline** | opens the airways | fast heart, ectopic beats, nausea, tremor; seizure risk above 40 mg/L | propranolol (slows the heart; but narrows asthmatic airways) | mg/L |
| **Midazolam** | sedation | deep sedation, slow breathing, low SpO₂, low blood pressure | flumazenil (antidote) | ng/mL |
| **Digoxin** | heart failure, atrial fibrillation | slow heart, heart block, ventricular ectopy, nausea | atropine (lifts the block) | ng/mL |
| **Morphine** | pain relief | sedation, slow breathing, low SpO₂, nausea | naloxone (antidote: pain returns); ondansetron (nausea) | ng/mL |
| **Naloxone** | reverses opioids | short-acting: morphine's effects can come back | – | ng/mL |
| **Flumazenil** | reverses benzodiazepines | short-acting: sedation can come back | – | ng/mL |
| **Propranolol** | slows the heart, lowers blood pressure | narrows the airways in asthma, blunts salbutamol | salbutamol (partly) | ng/mL |
| **Atropine** | speeds the heart | delirium, hot dry skin at high levels | – | ng/mL |
| **Ondansetron** | stops nausea | lengthens the QT interval (adds to digoxin's) | – | ng/mL |
| **Paracetamol** | mild pain relief, lowers fever | liver injury risk at high levels | – | mg/L |
| **Amlodipine** | lowers blood pressure slowly, for days | reflex rise in heart rate | propranolol | ng/mL |
| **Salbutamol** | opens the airways | fast heart, tremor | propranolol (but it cancels the airway effect) | ng/mL |

## Fourteen patients

Healthy man 35, healthy woman 30, athlete 24 (resting HR 50), smoker 45, obese woman 40 (BMI 38), asthma attack 28 (narrowed airways), after surgery 50 (pain 7/10), pneumonia with fever 64, elderly woman 84, frail man 91, liver cirrhosis 58, heart failure 72, kidney failure 66 (GFR 20), and a random patient from the mixed clinical population. Each has their own organ sizes, blood flows, liver and kidney function, so each clears every drug at a different speed; some also start with pain, fever or narrowed airways for the drugs to treat.

## The flows

| Flow | What it does |
|---|---|
| **1 · Virtual trial** | A population on a regimen (drug, dose, interval, route, loading dose, days). Shows the plasma level over time (5th to 95th percentile), the spread of steady-state levels, outcomes by group, and a table of every patient. |
| **2 · Find the dose** | Every population across the drug's dose range, with a curve per population of the share in range and the share at toxic levels, and the best dose for each. |
| **3 · Bedside: one patient, many drugs, live (3D)** | One patient, any number of drugs at once. Three windows: the 3D patient with the controls, the bedside monitor, and each drug's level over time. |
| **Maintenance · Make the 3D models** | Builds the patient from a CT case (below) and the scene. The window shows the model and lets you download the files for the project's Assets. |

### Using the bedside

The controls are at the bottom of the 3D view: **Patient**, **Drug**, **Dose ×** (a multiple of the usual dose), **Give dose**, **Regular doses** (the usual regimen for the chosen drug, or a continuous infusion for midazolam), **Stop all drugs**, **Clock** (15 min to 12 h of patient time per second), **Start again** and **Finish now**. Click **Sound on** in the 3D view to hear the heartbeat and the breathing.

A few things to try:

- *Elderly woman 84*, morphine, give three doses: breathing slows and SpO₂ falls, with alarms. Give naloxone: breathing recovers and pain (if any) comes back; naloxone wears off before morphine does.
- *Kidney failure 66*, digoxin with regular doses at 4 h per second: over a few days the level climbs, the heart slows and heart block appears. Atropine speeds the heart back up for a while.
- *Asthma attack 28*: salbutamol opens the airways and SpO₂ rises; propranolol then narrows them again.
- *After surgery 50*: morphine treats the pain but brings nausea; ondansetron stops the nausea and lengthens the QTc.
- *Pneumonia with fever 64*: paracetamol brings the temperature down.

### Reading the windows

- **Organ colour** (3D view) shows how much of the chosen drug each organ holds, on a log scale from indigo (hardly any) to pale yellow (the most); the legend on the left gives the levels. **Glow**: the gut absorbing a dose, the liver and kidneys clearing it.
- **The heart** contracts on each beat the model counts at the live heart rate (dropped beats in heart block included), with the heartbeat sound; **the lungs** fill on each breath at the live breathing rate, with the breath sound.
- **The bedside monitor** shows ECG, pulse oximetry and breathing waveforms, HR and rhythm, SpO₂, breathing rate, blood pressure, temperature, QTc, pain, sedation (Ramsay), nausea, airways, other effects, the drugs on board, and alarms.
- **Drug levels** shows each drug given as a percentage of the top of its target range (100 = the top).
- **The events list** tells the story: doses, levels entering or leaving range, rhythm, sedation, pain and nausea changes, alarms.

## What is modelled, and how

### Pharmacokinetics (where each drug goes): mechanistic

- **Compartments.** Lungs, brain, heart, kidneys, liver, gut, spleen, pancreas, fat, bone, muscle, skin and the rest of the body, plus arterial and venous blood and the gut contents. The gut, spleen and pancreas drain into the liver, so an oral dose passes the liver first.
- **Each patient's body.** Organ volumes and blood flows scale with body surface area. Fat comes from the Deurenberg formula, blood volume from Nadler's, and muscle is the weight left over. Heart failure lowers cardiac output, and with it liver blood flow.
- **Distribution.** Blood-flow limited. Tissue:plasma partition coefficients come from each tissue's composition by the Rodgers & Rowland method, scaled once so that a reference adult's volume of distribution matches the observed one.
- **Elimination.** Liver metabolism of the unbound drug (well-stirred liver), calibrated to each drug's observed clearance, plus renal excretion scaled by GFR. Smoking, cirrhosis and heart failure change the liver's capacity, drug by drug.
- **Solution.** The model is linear, so it is solved exactly over each 6-minute step with a matrix exponential. At the bedside each drug has its own model in the patient and in a healthy reference adult.

**Checked against the literature** for a reference 73 kg adult. Theophylline: half-life about 9 h (observed about 8 h). Midazolam: about 2.7 h (observed 1.8–6.4 h). Digoxin: about 40–50 h (observed about 40 h); 0.25 mg a day reaches about 1 ng/mL in a healthy adult and over 2 ng/mL with a GFR of 20–40 mL/min. The nine other drugs use typical adult values from their prescribing information and standard references; their clearances and volumes are approximate.

### Effects and interactions (pharmacodynamics): simplified

Each drug acts through receptor occupancy (an Emax model) on a few mechanisms, so drugs add up or cancel:

- **Antagonists compete:** naloxone with morphine, flumazenil with midazolam (competitive binding: more agonist needs more antidote).
- **Beta-blockade** (propranolol) lowers heart rate and blunts the heart-speeding effects of theophylline, salbutamol and the reflex to amlodipine; in asthma it narrows the airways and blunts salbutamol.
- **Breathing:** opioids and benzodiazepines depress it, together more than apart. Narrowed airways lower SpO₂ and raise the breathing rate.
- **Heart rhythm:** theophylline above 30 mg/L gives ectopic beats; digoxin above 2 ng/mL gives heart block, lifted by atropine; ondansetron and digoxin lengthen the QTc.
- **Symptoms:** pain (relieved by morphine, less by paracetamol), fever (paracetamol), nausea (theophylline, digoxin and morphine; ondansetron).

Half-effect levels and toxic thresholds follow published values where they exist; the shapes of the curves and the interaction terms are simplified and illustrative. The drugs do not change each other's kinetics (no enzyme inhibition or induction).

## The 3D patient

The patient is built from a real CT scan and its experts' labels: case s1400 of the **TotalSegmentator dataset** (v2.0.1, small subset). Every labelled structure goes into the model compartment it belongs to (the liver and gallbladder make the liver; oesophagus, stomach and bowel the gut; arteries and veins the blood; every bone the skeleton; skin and fat from the scan's densities). The scan covers the torso, so the brain is listed but has no surface.

**Credit (CC BY 4.0):** TotalSegmentator dataset v2.0.1 (small subset), Zenodo, doi:10.5281/zenodo.10047263. Wasserthal J., Breit H.-C., Meyer M.T., et al. *TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images.* Radiology: Artificial Intelligence 2023; 5(5): e230024.

**Sounds:** "Heartbeat (single)" by Universfield and "Deep breath" by hasin2004 (Pixabay).

## Sources

- **Physiology:** ICRP reference adult organ volumes and blood flows, as tabulated with PK-Sim's values in *Modeling principles for a physiology-based whole-body model of human metabolism* (arXiv:2404.14224, Table 1).
- **Tissue composition:** Ruark et al., as presented by Metrum Research Group (ACoP 2019, Table 1).
- **Partitioning:** Rodgers T, Rowland M. J Pharm Sci 2006;95:1238–57.
- **Theophylline:** *Clinical pharmacokinetics of theophylline.* Clin Pharmacokinet 1978;3:267–93. Abduljalil K et al. Front Pediatr 2022. Jusko WJ et al. J Pharm Sci 1979;68:1358–66.
- **Midazolam:** *Molecules* 2020;25:1080; the US prescribing information; Sunzel M et al. Br J Clin Pharmacol 1988.
- **Digoxin:** Iisalo E. *Clinical pharmacokinetics of digoxin.* Clin Pharmacokinet 1977;2:1–16.
- **The other nine drugs:** their prescribing information and standard pharmacology references (clearance, volume, binding, bioavailability, therapeutic and toxic levels).
- **Body composition:** Deurenberg P et al. (1991); Nadler SB et al. (1962); Mosteller (1987).

## Limits

- The nine added drugs use typical literature values; strong bases' binding to tissue phospholipids is not modelled, so their organ-by-organ pattern is approximate (their total volume is right). Clearances above liver blood flow are capped.
- Digoxin's binding to the sodium pump in muscle and heart is not modelled.
- Absorption is first order, elimination linear; drugs interact only through their effects, not their kinetics.
- Effect curves, interactions and baseline vital signs are simplified and illustrative.
- The 3D anatomy is one scanned person; it does not change with the virtual patient's sex, size or conditions.

<!-- esasim: files of the project -->

---

An [EsaSim](https://esacode.com/) project of the workspace Gianluca De Novi's workspace (Research).

### Blueprints in this repository

- [1 · Virtual trial](blueprints/1-virtual-trial.json) · 12 blocks · version 3
- [2 · Find the dose](blueprints/2-find-the-dose.json) · 7 blocks · version 1
- [3 · Bedside: one patient, many drugs, live (3D)](blueprints/3-bedside-one-patient-many-drugs-live-3d.json) · 5 blocks · version 5
- [Maintenance · Make the 3D models](blueprints/maintenance-make-the-3d-models.json) · 2 blocks · version 3

### Workspace blocks it uses

- [Virtual patients (PBPK dosing)](blocks/ws-virtual-patients-pbpk-dosing) · version 28
- [Bedside monitor](blocks/ws-bedside-monitor) · version 2
- [Model files (preview and download)](blocks/ws-model-files-preview-and-download) · version 1

### Assets

The files the project's views draw with (3D models, textures, sounds), in [assets/](assets/):

- [covers/fda-project-optimus-selecting-optimal-dose-oncology.webp](assets/covers/fda-project-optimus-selecting-optimal-dose-oncology.webp) · 355 KB
- [hasin2004-deep-breath-247459.mp3](assets/hasin2004-deep-breath-247459.mp3) · 148 KB
- [patient.glb](assets/patient.glb) · 12 MB
- [patient.scene.json](assets/patient.scene.json) · 5.4 KB
- [universfield-heartbeat-single-383748.mp3](assets/universfield-heartbeat-single-383748.mp3) · 35 KB

Keys pasted in blocks (secret settings) are not saved here, nor are the team and the results of runs.
