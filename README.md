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

This project models both. Each drug moves through 14 organs and the blood in each patient's own body (the pharmacokinetics, mechanistic); its levels at the sites where it acts set the vital signs and the side effects, with the drugs adding up or cancelling (the pharmacodynamics, simplified). The model is linear in each drug, so it is solved exactly, for one patient live or for thousands at once.

## Flows

### 1 · Virtual trial

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/104)

One population on one regimen. On the block, choose the **drug**, the **population**, the number of **virtual patients**, the share of **women**, the **route** (by mouth, sustained release, or IV infusion), the **dose** (in mg or mg/kg), how often (**every h**), an extra **loading dose** with the first one, and how many **days**. The View shows the plasma level over time (the median and the 5th to 95th percentiles against the target range and the toxic level), the spread of steady-state levels, who ends up in range, below, above or toxic in each group (smokers, age bands, sex, conditions), and a table of every patient.

### 2 · Find the dose

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/105)

Every population across the drug's whole dose range: for each, a curve of the share of patients in range and the share reaching toxic levels, and the best dose (the most patients in range while toxic levels stay under the limit you set, 5% by default).

### 3 · Bedside: one patient, many drugs, live (3D)

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/106)

| Block | What it does |
|---|---|
| **The bedside** | The patient and every drug given, simulated live, with a healthy adult on the same doses for comparison. Stays open 20 minutes (Bedside: minutes it stays open). |
| **The patient, live** | Sends the patient, the vitals and the events to the views. |
| **The patient (3D)** | The patient built from a CT scan, the drugs in the organs, the heart beating and the chest breathing, the side effects as they come and go; the controls; what to do next. |
| **Bedside monitor** | ECG, pulse oximetry and breathing waveforms, the vital signs and the alarms. |
| **Drug levels** | Each drug given, as a share of the top of its target range, over the last 90 seconds. |

### Maintenance · Make the 3D models

[**Open the blueprint →**](https://esacode.com/app/projects/29/flows/103)

Rebuilds the 3D patient from a CT case (s1400 of the TotalSegmentator dataset by default) and shows it in a window with a **Download** button. Only needed if the model changes: see *The 3D patient* below.

## Using the bedside

Run flow 3 and open the **View**. Everything is done from the bar under the 3D patient:

| Control | What it does |
|---|---|
| **Patient** | Changes the patient: the clock starts again at Day 1, with no drugs. |
| **Drug** | The drug that **Give dose** and **Regular** act on, and the one whose levels colour the organs. |
| **Dose ×** | A multiple of the usual dose (¼ to 4 times) for the chosen drug, for the next doses. |
| **Give dose** | One usual dose of the chosen drug (by mouth or IV, as it is usually given). |
| **Regular** | The usual regimen of the chosen drug, dose after dose (midazolam: a continuous infusion); press again to stop it. |
| **Stop all** | Stops every regular dose and infusion; what is already in the body stays and is cleared. |
| **Clock** | How fast the patient's time runs: 15 minutes, 1 hour, 4 hours or 12 hours of the patient's time a second. Faster clocks show days of regular dosing in a minute. |
| **Start again** | The same patient from Day 1, with no drugs. |

After each control a message appears across the middle of the 3D view for a few seconds (*Morphine 4 mg IV given*, *Started again: Healthy man 35, Day 1, no drugs*…). It takes a second or two: the control goes to the simulation on the worker and comes back.

**What to do next.** The yellow box at the top of the 3D view reads the alarms and the drugs on board and suggests the next step: at the start, to choose a drug and give it; when morphine slows the breathing, to give naloxone; when digoxin blocks the heart, to give atropine and stop the regular doses; when propranolol narrows an asthmatic's airways, to stop it and give salbutamol; pain, fever, wheeze and nausea each point to the drug that treats them; a racing heart, low blood pressure or a long QTc to what to slow, stop or avoid. When nothing is wrong it says so and suggests raising the dose, adding another drug or a faster clock.

**Sound and the mouse.** Click **Sound on** (top right of the 3D view) to hear the heartbeat and the breathing. Drag with the left button to turn the patient in your hand, double-click to put it back; the wheel zooms and the right button moves the view.

### A few things to try

- *Elderly woman 84*, Morphine, **Give dose** three times: breathing slows and SpO₂ falls, with alarms, and the side effects list them. Then **Naloxone**: she breathes again, until it wears off before the morphine does.
- *Kidney failure 66*, Digoxin, **Regular**, clock 4 h/s: over a few days the level climbs, the heart slows, heart block and yellow-tinged vision appear. **Atropine** lifts the block for a while; **Stop all** lets the digoxin clear, slowly.
- *Asthma attack 28*: **Salbutamol** opens the airways and SpO₂ rises (with a tremor); **Propranolol** narrows them again.
- *After surgery 50*: **Morphine** treats the pain but brings nausea, and constipation if kept on **Regular** for a day; **Ondansetron** stops the nausea and lengthens the QTc.
- *Pneumonia with fever 64*: **Paracetamol** brings the temperature down.
- *Smoker 45* against *Liver cirrhosis 58*, Theophylline **Regular**: the smoker stays low, the cirrhotic climbs into the toxic range (restlessness, fast heart, extra beats, nausea).
- *Heart failure 72*, Amlodipine **Regular**, clock 12 h/s: blood pressure falls over days, and after three days the ankles swell.

### Reading the windows

**The 3D patient**
- **Organ colours**: how much of the chosen drug each organ holds, on a log scale from indigo (hardly any) through violet, magenta and orange to pale yellow (the most). The gut glows while it absorbs a dose, the liver and the kidneys while they clear the drug.
- **Top left, Virtual patient**: the patient, the day and time, and the drugs on board with their levels (*in range*, *low*, *high* or the toxic effect).
- **Left, Where the chosen drug is**: its plasma level against its target and a healthy adult's; this patient's clearance and half-life against a healthy adult's; one bar per organ, on the same scale as the colours, from 0% (hardly any) to 100% (the most), with what 0%, 25%, 50%, 75% and 100% mean in the drug's unit (mg/L or ng/mL) underneath.
- **Top right, Side effects**: each side effect as it starts, top to bottom in the order they appear, with where in the body it acts, the drug causing it and since when (●); when it passes it stays a few hours with its start and end (○).
- **Bottom right**: the events, as they happen (doses, levels entering or leaving the range, side effects starting and passing, alarms).
- **The heart** contracts on every beat the model counts at the live heart rate (a dropped beat in heart block is skipped), with the heartbeat sound. **The lungs and the ribcage** rise on every breath at the live breathing rate, with the breath sound.

**The bedside monitor**
- The ECG, pulse oximetry (pleth) and breathing waveforms, drawn in real time at the patient's current rates. Ectopic beats, dropped beats and a long QT show on the ECG.
- HR with the rhythm, SpO₂, breathing rate (RESP), blood pressure (NIBP, with the mean), temperature, QTc, pain (0–10) and sedation (Ramsay 1–6).
- Below: nausea, the airways, other effects, and the drugs on board. The alarms show at the top in red (dangerous) or amber.
- **Beep** turns on a pulse-oximetry beep, its pitch falling with SpO₂.

**Drug levels**: each drug given, as a percentage of the top of its target range (100 = the top of the range), over the last 90 seconds of the run. Drag on it to zoom in; double-click to follow it live again.

## The drugs

| Drug | What it does | Side effects in the model | Corrected by |
|---|---|---|---|
| **Theophylline** | opens the airways | restlessness, fast heart, extra beats, nausea, tremor; seizure risk above 40 mg/L | propranolol (but it narrows asthmatic airways) |
| **Midazolam** | sedation | deep sedation, slow breathing, low SpO₂, low blood pressure | flumazenil |
| **Digoxin** | heart failure, atrial fibrillation | slow heart, heart block, extra beats, nausea, yellow-tinged vision, longer QTc | atropine |
| **Morphine** | pain relief | sedation, slow breathing, low SpO₂, nausea; constipation after a day | naloxone; ondansetron for the nausea |
| **Naloxone** | reverses opioids | sudden withdrawal (pain, agitation) at high doses; wears off before morphine | – |
| **Flumazenil** | reverses benzodiazepines | wears off before midazolam | – |
| **Propranolol** | slows the heart, lowers blood pressure | narrows the airways in asthma, blunts salbutamol; tiredness and cold hands after a day | – |
| **Atropine** | speeds the heart | delirium, hot dry skin at high levels | – |
| **Ondansetron** | stops nausea | longer QTc | – |
| **Paracetamol** | mild pain relief, lowers fever | liver injury risk above 150 mg/L | – |
| **Amlodipine** | lowers blood pressure slowly, for days | a reflex rise in heart rate; ankle swelling after three days | propranolol |
| **Salbutamol** | opens the airways | fast heart, tremor | – |

Theophylline and paracetamol are in mg/L, the others in ng/mL.

### When the side effects appear

Most side effects follow the drug's level at the site where it acts, so they appear when the level crosses a threshold and pass when it falls back: within minutes for IV drugs (morphine, midazolam, naloxone, salbutamol), within an hour or two for drugs by mouth, and over days for drugs that build up slowly with regular doses (digoxin, amlodipine, theophylline in a slow clearer). Three come only after continued use: **constipation** after a day on morphine, **tiredness and cold hands** after a day on propranolol, **ankle swelling** after three days on amlodipine; they fade as slowly when the drug is stopped.

## The patients

| Patient | What makes them different |
|---|---|
| Healthy man 35, healthy woman 30 | the reference adults |
| Athlete 24 | resting heart rate 50 |
| Smoker 45 | smoking speeds the liver's clearance of theophylline and propranolol |
| Obese woman 40 (BMI 38) | more fat: fat-loving drugs (midazolam, amlodipine) spread wider and last longer |
| Asthma attack 28 | narrowed airways, low SpO₂: for salbutamol and theophylline, and against propranolol |
| After surgery 50 | pain 7/10: for morphine and paracetamol |
| Pneumonia with fever 64 | fever 39.3 °C, some pain and wheeze |
| Elderly woman 84, frail man 91 | small, with weaker kidneys and liver: higher levels on usual doses |
| Liver cirrhosis 58 | slower liver clearance |
| Heart failure 72 | less blood flow to the liver, faster resting heart and breathing |
| Kidney failure 66 (GFR 20) | digoxin and the kidney-cleared drugs build up |
| Random patient | drawn from the mixed clinical population |

The trials draw thousands from seven populations: healthy adults, smokers, elderly (65+), liver cirrhosis, heart failure, kidney impairment, and a mixed clinical population.

## What you need

- A machine with the EsaSim worker (Project › Compute). Any machine: an NVIDIA GPU speeds up the large trials but is not needed.
- Python packages, installed by the worker the first time: `numpy`, `torch`, `scipy`, and for the 3D models `SimpleITK`, `scikit-image`, `trimesh`, `fast-simplification`, `remotezip`, `requests`.
- Internet access on the machine only to make the 3D models (the CT case is read from Zenodo).
- No account or key.

## How to run it

1. Under **Compute**, check that a machine is available for this project.
2. Open a flow under **Blueprint**, set its block (drug, population, regimen; for the bedside, the first patient) and press **Run**.
3. Open the **View**. For the bedside, see *Using the bedside* above.

## What you get

- **Virtual trial**: `course.csv` (the plasma level over time, by percentile), `patients.csv` (every patient: who they are, clearance, half-life, trough, peak, average, outcome), `groups.csv` (in range, below, above and toxic, by group), with their charts and tables in the View.
- **Find the dose**: `dose_finding.csv` (every dose tried for every population) and `best_doses.csv`.
- **Bedside**: the live views above; nothing is saved.

## What is computed

| Part | How |
|---|---|
| **Organs** | lungs, brain, heart, kidneys, liver, gut, spleen, pancreas, fat, bone, muscle, skin, the rest of the body, arterial and venous blood, the gut contents; the gut, spleen and pancreas drain into the liver |
| **Each body** | organ volumes and blood flows from body surface area; fat (Deurenberg), blood volume (Nadler), muscle the weight left over; heart failure lowers cardiac output |
| **Distribution** | blood-flow limited; tissue:plasma partition from tissue composition (Rodgers & Rowland), scaled to the observed volume of distribution |
| **Elimination** | well-stirred liver on the unbound drug, calibrated to the observed clearance, plus renal excretion scaled by GFR; smoking, cirrhosis and heart failure change the liver's capacity, drug by drug |
| **Solution** | exact 6-minute steps with a matrix exponential (PyTorch, CPU or GPU) |
| **Effects** | receptor occupancy (Emax) at each drug's effect site (an effect compartment where the effect lags the plasma; digoxin's level in the heart itself); competitive antagonists (naloxone, flumazenil); beta-blockade; breathing depression by opioids and benzodiazepines together; airway narrowing; heart rhythm (ectopy, heart block, QTc); pain, fever and nausea |
| **Side effects** | 21 side effects, each tied to where it acts in the body and to the drugs that cause it: present while the vitals or the drug's level say so, the three late ones after the hours of exposure they need |
| **Heartbeats and breaths** | counted in real time at the live heart and breathing rates, so the 3D heart, the chest, the sounds and the monitor all follow the same patient |

Checked against the literature for a 73 kg adult: theophylline half-life about 9 h (observed about 8 h), midazolam about 2.7 h (1.8–6.4 h), digoxin about 40–50 h (about 40 h); digoxin 0.25 mg a day gives about 1 ng/mL in a healthy adult and over 2 ng/mL with a GFR of 20–40 mL/min.

## The 3D patient

The patient is a real CT scan and its experts' labels (case s1400 of the TotalSegmentator dataset): each labelled structure goes into the organ of the model it belongs to (the liver and gallbladder make the liver; oesophagus, stomach and bowel the gut; arteries and veins the blood; the ribs, sternum and costal cartilages the ribcage, which moves with the breath; every other bone the skeleton), and the skin and the fat under it come from the scan's own densities.

Its files are in the project's assets:

| Asset | What it is |
|---|---|
| `patient split.glb` | the 3D patient, one part per organ |
| `patient.scene.json` | the 3D view's scene: the colours, the heart and the breathing, the sounds, the controls, the panels |
| `patient-turn.js` | lets the mouse turn the patient instead of the camera |
| `universfield-heartbeat-single-383748.mp3`, `hasin2004-deep-breath-247459.mp3` | the heartbeat and the breath |

To rebuild the model, run **Maintenance · Make the 3D models**, press **Download** next to `patient.glb` in its window, and upload it to the assets in place of `patient split.glb` (under that name, or change the `model` line of `patient.scene.json`).

## Data, method and licences

- **3D patient**: case s1400 of the **TotalSegmentator dataset v2.0.1, small subset**, [Zenodo record 10047263](https://zenodo.org/records/10047263), licence **CC BY 4.0**. Wasserthal J., Breit H.-C., Meyer M.T., et al. *TotalSegmentator: Robust Segmentation of 104 Anatomic Structures in CT Images.* Radiology: Artificial Intelligence 2023; 5(5): e230024.
- **Sounds**: "Heartbeat (single)" by Universfield and "Deep breath" by hasin2004, from Pixabay.
- **Physiology**: ICRP reference adult organ volumes and blood flows, as tabulated with PK-Sim's values in *Modeling principles for a physiology-based whole-body model of human metabolism* (arXiv:2404.14224). Tissue composition after Ruark et al. (Metrum Research Group, ACoP 2019). Rodgers T, Rowland M. J Pharm Sci 2006;95:1238–57. Deurenberg P et al. (1991); Nadler SB et al. (1962); Mosteller (1987).
- **Theophylline**: Clin Pharmacokinet 1978;3:267–93; Abduljalil K et al. Front Pediatr 2022; Jusko WJ et al. J Pharm Sci 1979;68:1358–66.
- **Midazolam**: Molecules 2020;25:1080; the US prescribing information; Sunzel M et al. Br J Clin Pharmacol 1988.
- **Digoxin**: Iisalo E. Clin Pharmacokinet 1977;2:1–16.
- **The other nine drugs and the side effects**: their prescribing information and standard pharmacology references.
- The code of the blocks is EsaSim's.

## Notes

- The pharmacokinetics are mechanistic; the effects, the interactions, the side effects' thresholds and timings, the suggestions and the baseline vital signs are simplified and illustrative, following published thresholds and half-effect levels where they exist.
- Side effects are switched on and off by thresholds; their severity is shown by the vitals, not graded in the list. Inflammation and allergic reactions are not modelled.
- The nine added drugs use typical literature values; strong bases' binding to tissue phospholipids is not modelled, so their organ-by-organ pattern is approximate (their total volume is right). Clearances above the liver's blood flow are capped.
- Drugs interact through their effects only, not their kinetics (no enzyme inhibition or induction). Absorption is first order, elimination linear. Digoxin's binding to the sodium pump is not modelled.
- The 3D body is one scanned person, neck to pelvis: it does not change with the virtual patient, and the brain has no surface. The breathing moves the ribcage and the lungs only; the skin stays still.
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

- [Virtual patients (PBPK dosing)](blocks/ws-virtual-patients-pbpk-dosing) · version 55
- [Bedside monitor](blocks/ws-bedside-monitor) · version 3
- [Model files (preview and download)](blocks/ws-model-files-preview-and-download) · version 1

### Assets

The files the project's views draw with (3D models, textures, sounds), in [assets/](assets/):

- [covers/fda-project-optimus-selecting-optimal-dose-oncology.webp](assets/covers/fda-project-optimus-selecting-optimal-dose-oncology.webp) · 355 KB
- [hasin2004-deep-breath-247459.mp3](assets/hasin2004-deep-breath-247459.mp3) · 148 KB
- [patient split.glb](assets/patient%20split.glb) · 12 MB
- [patient-turn.js](assets/patient-turn.js) · 1.6 KB
- [patient.glb](assets/patient.glb) · 12 MB
- [patient.scene.json](assets/patient.scene.json) · 8.2 KB
- [universfield-heartbeat-single-383748.mp3](assets/universfield-heartbeat-single-383748.mp3) · 35 KB

Keys pasted in blocks (secret settings) are not saved here, nor are the team and the results of runs.
