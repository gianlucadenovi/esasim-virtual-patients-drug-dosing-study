# ── EsaSim block ──────────────────────────────────────────────────────
# Virtual patients (PBPK dosing) (ws-virtual-patients-pbpk-dosing) · version 64 · worker block, Python: def run(ctx)
# A whole-body pharmacokinetic model of twenty drugs in virtual patients: trials, the best dose per population, a live bedside, and an anaesthesia room where a case is run and scored
# Workspace: Gianluca De Novi's workspace · Author: Gianluca De Novi
# Saved: 2026-11-25
#
# Inputs   ctx.inputs["id"] is a list of files: .path .name .size .read_text() .read_json()
#   ct          Files (gz, nii)           CT case · several · optional · Make the 3D models: a CT case with its labels (ct.nii.gz and segmentations/<structure>.nii.gz, as the TotalSegmentator dataset gives them). Without it, the CT case setting, else a mockup body.
# Outputs  ctx.output("id", "name.ext") is the path to write; ctx.write_json / ctx.write_text
#   course      Table (csv)               Course · course.csv: the plasma concentration over time, percentiles over the patients (p05, p25, median, p75, p95), with the target range and the toxic level
#   patients    Table (csv)               Patients · patients.csv: one row per virtual patient (who they are, clearance, half-life, trough, peak, average, outcome)
#   groups      Table (csv)               By group · groups.csv: in range, below, above and toxic peaks for each group (smokers, age bands, sex, conditions)
#   doses       Table (csv)               Doses tried · Find the dose: dose_finding.csv, every dose tried for every population
#   best        Table (csv)               Best doses · Find the dose: best_doses.csv, the dose chosen for each population
#   models      Files (glb, json)         3D models · Make the 3D models: patient.glb and patient.scene.json, for the project's assets
#   telemetry   Live telemetry            ctx.telemetry("channel", data) while it runs
# Settings ctx.params["key"]
#   mode        select    What to do = "Run a virtual trial" of Run a virtual trial | Find the dose | Bedside (live) | Anaest…
#   drug        select    Drug = "Theophylline" of Theophylline | Midazolam | Digoxin | Morphine | Naloxone | Flumazenil | Pr…
#   population  select    Population = "Healthy adults" of Healthy adults | Smokers | Elderly (65+) | Liver cirrhosis | Heart…
#   patients    number    Virtual patients = 2000
#   female_pct  number    Women (%) = 50
#   route       select    Route = "Oral (immediate release)" of Oral (immediate release) | Oral (sustained release) | IV infu…
#   dose_mg     number    Dose = 300
#   dose_basis  select    Dose in = "mg" of mg | mg/kg
#   interval_h  number    Every (h) = 12
#   loading_mg  number    Extra with the first dose = 0
#   infusion_h  number    IV infusion over (h) = 0.5
#   days        number    Days = 5
#   target      select    Target range (theophylline) = "10-20 mg\/L (classic)" of 10-20 mg/L (classic) | 5-15 mg/L (current …
#   population_scope select    Find the dose for = "All populations" of All populations | The population chosen
#   max_toxic_pct number    Find the dose: toxic levels at most (%) = 5
#   bedside_patient select    Bedside: first patient = "Healthy man 35" of Healthy man 35 | Healthy woman 30 | Athlete 24 (restin…
#   bedside_minutes number    Bedside: minutes it stays open = 20
#   an_scenario select    Anaesthesia: the job = "Operation: a four-hour laparotomy" of Operation: a four-hour laparotomy | I…
#   an_patient  select    Anaesthesia: the patient = "Healthy man 35" of Healthy man 35 | Healthy woman 30 | Athlete 24 (rest…
#   an_minutes  number    Anaesthesia: minutes it stays open = 30
#   ct_case     text      Make the 3D models: CT case = ""
#   accelerator select    Run on = "auto" of auto | NVIDIA GPU | CPU
#   seed        number    Seed = 1
# Runs on  Any worker, NVIDIA GPU. Return a dict to fill the "result" output.
# Packages numpy, torch, scikit-image, scipy, SimpleITK, trimesh, fast-simplification, remotezip, requests (installed by the worker with pip when missing)
# EsaSim rewrites these lines from the block's definition each time it is saved: change the definition, not them.
# ── end of EsaSim block header ────────────────────────────────────────

"""Virtual patients: a whole-body physiologically based pharmacokinetic (PBPK) model of a drug in a virtual
population, the trial a dosing regimen gives (who is in the target range, who is under- or over-exposed, who reaches
toxic peaks), the dose that suits each population best, and one patient followed live in the 3D view (each organ
coloured by the drug's concentration in it)."""
import json
import math
import struct

import numpy as np


# ── the body: reference organ volumes (L) and blood flows (L/min) ──────────────────────────────────────────────────
# ICRP reference adults (male 73 kg, 176 cm; female 60 kg, 163 cm), as tabulated with PK-Sim's values in "Modeling
# principles for a physiology-based whole-body model of human metabolism" (arXiv:2404.14224, Table 1).
# The stomach, small and large intestine are one gut compartment; the gonads go with the rest of the body.
REF = {   # organ: (V male, V female, Q male, Q female)
    "adipose": (14.65, 19.42, 0.320, 0.503),
    "brain": (1.509, 1.355, 0.780, 0.706),
    "bone": (11.82, 9.122, 0.325, 0.295),
    "heart": (0.417, 0.329, 0.260, 0.295),
    "kidney": (0.438, 0.404, 1.327, 1.122),
    "liver": (2.377, 1.918, 0.426, 0.386),          # (its own arterial flow; the portal flow comes from the gut, spleen, pancreas)
    "gut": (0.413 + 0.725 + 0.169, 0.417 + 0.695 + 0.163, 0.260 + 0.650 + 0.065, 0.295 + 0.649 + 0.059),
    "lung": (1.215, 0.946, None, None),             # (all of the cardiac output)
    "muscle": (32.65, 20.29, 1.116, 0.666),
    "pancreas": (0.190, 0.170, 0.065, 0.059),
    "skin": (3.763, 2.724, 0.325, 0.296),
    "spleen": (0.207, 0.186, 0.166, 0.150),
    "rest": (0.040, 0.013, 0.003, 0.001),           # (the gonads)
}
REF_BLOOD = {"male": (0.419, 0.964 + 1.037), "female": (0.320, 0.737 + 0.793)}   # arterial, venous (with the portal vein)
REF_CO = {"male": 6.088, "female": 5.482}           # cardiac output (L/min)
REF_BODY = {"male": (73.0, 176.0), "female": (60.0, 163.0)}   # kg, cm

# ── tissue composition, for the partition coefficients ───────────────────────────────────────────────────────────
# Unified tissue composition (Ruark et al., as tabulated by Metrum Research Group, ACoP 2019, Table 1):
# neutral lipid, neutral phospholipid, extracellular and intracellular water fractions, intracellular pH, and the
# albumin tissue:plasma ratio (AR). The gut row stands for the pancreas too, the muscle row for the rest.
COMP = {   # tissue: (f_nl, f_npl, f_ew, f_iw, pH, AR)
    "bone": (0.074, 0.0016, 0.100, 0.346, 7.00, 0.100),
    "brain": (0.045, 0.0553, 0.162, 0.620, 7.10, 0.048),
    "adipose": (0.798, 0.0478, 0.135, 0.017, 7.10, 0.049),
    "heart": (0.089, 0.0079, 0.320, 0.456, 7.10, 0.157),
    "kidney": (0.036, 0.0166, 0.273, 0.483, 7.22, 0.130),
    "gut": (0.0487, 0.0124, 0.282, 0.475, 7.40, 0.158),
    "liver": (0.037, 0.0115, 0.161, 0.573, 7.23, 0.086),
    "lung": (0.003, 0.0056, 0.336, 0.446, 6.60, 0.212),
    "muscle": (0.013, 0.0092, 0.118, 0.630, 6.81, 0.064),
    "skin": (0.036, 0.0502, 0.382, 0.291, 7.00, 0.277),
    "spleen": (0.014, 0.0103, 0.207, 0.579, 7.00, 0.097),
}
COMP["pancreas"] = COMP["gut"]
COMP["rest"] = COMP["muscle"]
PLASMA_NL, PLASMA_NPL, PLASMA_PH = 0.003, 0.005, 7.4

TISSUES = ["lung", "brain", "heart", "kidney", "liver", "gut", "spleen", "pancreas", "adipose", "bone", "muscle", "skin", "rest"]
PORTAL = ["gut", "spleen", "pancreas"]               # (they drain into the liver)
# the state: the amount (mg) in each tissue, in arterial and venous blood, and in the gut's lumen (an oral dose)
STATE = TISSUES + ["arterial", "venous", "lumen"]
IX = {k: i for i, k in enumerate(STATE)}
NS = len(STATE)


# ── drugs ────────────────────────────────────────────────────────────────────────────────────────────────────────
# Three drugs, each with a different story: theophylline (cleared by the liver's CYP1A2, sped up by smoking, slowed by
# liver disease and heart failure), midazolam (a lipophilic sedative cleared by the liver's CYP3A4, stored in fat;
# it slows breathing and lowers blood pressure) and digoxin (cleared mostly unchanged by the kidneys, slow to reach
# the heart; it slows the heart and becomes toxic when the kidneys fail).
# Concentrations are in each drug's own unit (mg/L for theophylline, ng/mL for the others): "unit_k" converts the
# model's mg/L to it. "factors": how smoking, cirrhosis and heart failure change the hepatic intrinsic clearance (heart
# failure's lower cardiac output and liver blood flow also act through the physiology, for every drug).
#
# Theophylline: logP -0.02 and pKa 8.6 (a very weak acid; PKKB); 40 % bound in plasma, mostly to albumin (fu 0.6);
# oral bioavailability about 96-100 %, peak 0.5-2 h after an immediate-release tablet (Clin Pharmacokinet 1978;
# 3:267-93); about 85 % of its clearance hepatic (CYP1A2 about 88 % of that), the renal clearance 0.31 L/h in adults
# (Abduljalil et al. 2022); Vss about 0.5 L/kg and a typical adult non-smoker's clearance about 0.04 L/h/kg.
# Concentrations: adverse effects become common above 20 mg/L, arrhythmias between 20 and 40 mg/L, seizures above
# 40 mg/L.
# Midazolam: logP 4.33, pKa 5.5 (a weak base; Molecules 2020;25:1080, Table 1); about 97 % bound to albumin (fu 0.03);
# total clearance 0.25-0.54 L/h/kg, Vd 1.0-3.1 L/kg, half-life 1.8-6.4 h (mean about 3 h) in healthy adults; clearance
# lower with old age, heart failure and cirrhosis; Vd larger in women, the elderly and the obese (US prescribing
# information). Blood:plasma ratio 0.66 and the fraction reaching the blood after the gut wall's CYP3A4 (0.75): typical
# model values (assumptions). Sedation: Ramsay scores 2-6 at about 6, 71, 171, 260 and 659 ng/mL (a population
# pharmacodynamic study of midazolam by target-controlled infusion in surgical ICU patients); blood pressure and
# ventilation (PaCO2) half-maximally affected at 50-60 ng/mL (Sunzel et al., Br J Clin Pharmacol 1988).
# Digoxin: logP 1.26 (neutral; PKKB, QsarDB); 25 % bound (fu 0.75); oral bioavailability about 70 %; clearance about
# 6 L/h and Vd about 400 L (about 6 L/kg) in adults, half-life about 40 h, mostly renal (25-28 % non-renal; Iisalo,
# Clin Pharmacokinet 1977); high concentrations in the heart and kidneys, skeletal muscle the largest store, fat no
# reservoir; therapeutic 0.5-2 ng/mL, toxic effects (bradycardia, heart block, ventricular arrhythmias) from about
# twice that. Its effects follow its uptake in the heart, after a 6-8 h distribution phase: here the heart's own
# concentration (as the plasma concentration it is in balance with) drives them. (Its binding to the sodium pump in
# muscle and heart is not modelled: the tissues' levels come from their lipids and water, scaled to the observed Vss.)
DRUGS = {
    "Theophylline": {
        "logP": -0.02, "pKa": 8.6, "kind": "acid", "fu": 0.60, "bp": 1.0, "F": 0.96,
        "ka": {"Oral (immediate release)": 1.5, "Oral (sustained release)": 0.25},
        "cl_per_kg": 0.040,        # (L/h/kg) a typical healthy adult non-smoker's total clearance
        "cl_renal": 0.31,          # (L/h) at a GFR of 120 mL/min
        "vss_per_kg": 0.50,        # (L/kg) the observed volume the tissue partitioning is scaled to
        "unit": "mg/L", "unit_k": 1.0,
        "factors": {"smoking": 1.6, "cirrhosis": 0.40, "chf": 0.55},
        "range": {"10-20 mg/L (classic)": (10.0, 20.0), "5-15 mg/L (current guidance)": (5.0, 15.0)},
        "toxic": [(20.0, "adverse effects (nausea, tachycardia)"), (30.0, "arrhythmia risk"), (40.0, "seizure risk")],
        "dose_steps": (25.0, 2000.0, 25.0),
        "effect_at": "plasma", "mw": 180.17,
    },
    "Midazolam": {
        "logP": 4.33, "pKa": 5.5, "kind": "base", "fu": 0.03, "bp": 0.66, "F": 0.75,
        "ka": {"Oral (immediate release)": 3.0, "Oral (sustained release)": 3.0},
        "cl_per_kg": 0.40, "cl_renal": 0.0, "vss_per_kg": 1.5,
        "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"smoking": 1.0, "cirrhosis": 0.50, "chf": 1.0},
        "range": {"Sedation 70-260 ng/mL": (70.0, 260.0)},
        "toxic": [(260.0, "deep sedation"), (659.0, "unrousable: airway at risk")],
        "dose_steps": (0.5, 30.0, 0.5),
        "ramsay_c50": [5.7, 71.0, 171.0, 260.0, 659.0],
        "effect_at": "plasma", "mw": 325.77,
    },
    "Digoxin": {
        "logP": 1.26, "pKa": 7.0, "kind": "neutral", "fu": 0.75, "bp": 1.0, "F": 0.70,
        "ka": {"Oral (immediate release)": 1.0, "Oral (sustained release)": 1.0},
        "cl_per_kg": 6.0 / 73.0, "cl_renal": 4.4, "vss_per_kg": 6.0,
        "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"smoking": 1.0, "cirrhosis": 1.0, "chf": 1.0},
        "range": {"0.5-2 ng/mL": (0.5, 2.0)},
        "toxic": [(2.0, "toxic: bradycardia, heart block"), (3.0, "ventricular arrhythmias")],
        "dose_steps": (0.0625, 1.0, 0.03125),
        "effect_at": "heart", "mw": 780.94,
    },
}


# The other nine drugs: typical adult values from their prescribing information and standard pharmacology references
# (approximate). Clearances above about 1.1 L/h/kg are capped there (the model's liver cannot clear faster than its
# blood flow; naloxone's real clearance is partly outside the liver). For strong bases (pKa above 7) the binding to
# acidic phospholipids is not modelled, so their organ-by-organ pattern is approximate (their total volume is right).
# "give": one dose from the bedside (mg, route, over how many hours); "regular": the usual regimen (mg, every h, route,
# over h). "ke0" (per hour): how fast the effect follows the plasma (an effect compartment); none: at once.
DRUGS.update({
    "Morphine": {   # opioid analgesic: pain relief; sedation, slow breathing, nausea; reversed by naloxone
        "logP": 0.89, "pKa": 8.2, "kind": "base", "fu": 0.65, "bp": 1.0, "F": 1.0, "ka": {"Oral (immediate release)": 1.5},
        "cl_per_kg": 0.85, "cl_renal": 0.0, "vss_per_kg": 3.5, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.6, "chf": 1.0, "smoking": 1.0}, "range": {"Analgesia 20-80 ng/mL": (20.0, 80.0)},
        "toxic": [(150.0, "breathing depressed")], "dose_steps": (1.0, 40.0, 1.0), "ke0": 1.0, "ec50": 30.0,
        "give": (4.0, "iv", 2 / 60), "regular": (4.0, 4.0, "iv", 2 / 60), "mw": 285.3},
    "Naloxone": {   # opioid antagonist: reverses morphine, for about an hour
        "logP": 2.09, "pKa": 7.94, "kind": "base", "fu": 0.55, "bp": 1.0, "F": 0.05, "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 1.1, "cl_renal": 0.0, "vss_per_kg": 2.5, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.7, "chf": 1.0, "smoking": 1.0}, "range": {"Reversal 1-10 ng/mL": (1.0, 10.0)},
        "toxic": [(50.0, "abrupt reversal: pain, agitation")], "dose_steps": (0.1, 2.0, 0.1), "ke0": 6.0, "ki": 1.5,
        "give": (0.4, "iv", 1 / 60), "regular": (0.4, 1.0, "iv", 1 / 60), "mw": 327.4},
    "Flumazenil": {   # benzodiazepine antagonist: reverses midazolam, for about an hour
        "logP": 1.0, "pKa": 1.7, "kind": "neutral", "fu": 0.5, "bp": 1.0, "F": 0.16, "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.9, "cl_renal": 0.0, "vss_per_kg": 1.0, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.5, "chf": 1.0, "smoking": 1.0}, "range": {"Reversal 10-40 ng/mL": (10.0, 40.0)},
        "toxic": [(200.0, "seizure risk")], "dose_steps": (0.1, 2.0, 0.1), "ke0": 6.0, "ki": 1.0,
        "give": (0.2, "iv", 0.25 / 60), "regular": (0.2, 1.0, "iv", 0.25 / 60), "mw": 303.3},
    "Propranolol": {   # beta-blocker: slows the heart, lowers blood pressure; tightens the airways in asthma
        "logP": 3.09, "pKa": 9.5, "kind": "base", "fu": 0.10, "bp": 0.8, "F": 1.0, "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.75, "cl_renal": 0.0, "vss_per_kg": 4.0, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.5, "chf": 1.0, "smoking": 1.2}, "range": {"20-100 ng/mL": (20.0, 100.0)},
        "toxic": [(300.0, "bradycardia, hypotension")], "dose_steps": (10.0, 320.0, 10.0), "ec50": 18.2,
        "give": (1.0, "iv", 1 / 60), "regular": (40.0, 8.0, "oral", 0.0), "mw": 259.3},
    "Atropine": {   # antimuscarinic: speeds the heart, lifts vagal heart block
        "logP": 1.83, "pKa": 9.8, "kind": "base", "fu": 0.56, "bp": 1.0, "F": 0.9, "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.5, "cl_renal": 18.0, "vss_per_kg": 2.5, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.7, "chf": 1.0, "smoking": 1.0}, "range": {"2-8 ng/mL": (2.0, 8.0)},
        "toxic": [(20.0, "delirium, hot dry skin")], "dose_steps": (0.1, 3.0, 0.1), "ec50": 3.0,
        "give": (0.5, "iv", 0.5 / 60), "regular": (0.5, 4.0, "iv", 0.5 / 60), "mw": 289.4},
    "Ondansetron": {   # antiemetic: stops nausea; prolongs the QT interval
        "logP": 2.4, "pKa": 7.4, "kind": "base", "fu": 0.27, "bp": 1.0, "F": 0.6, "ka": {"Oral (immediate release)": 1.2},
        "cl_per_kg": 0.35, "cl_renal": 0.0, "vss_per_kg": 1.9, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.4, "chf": 1.0, "smoking": 1.0}, "range": {"20-100 ng/mL": (20.0, 100.0)},
        "toxic": [(300.0, "QT prolongation")], "dose_steps": (1.0, 32.0, 1.0), "ec50": 20.0,
        "give": (4.0, "iv", 2 / 60), "regular": (4.0, 8.0, "iv", 2 / 60), "mw": 293.4},
    "Paracetamol": {   # analgesic and antipyretic: mild pain relief, lowers fever; liver toxicity at high levels
        "logP": 0.46, "pKa": 9.4, "kind": "acid", "fu": 0.8, "bp": 1.0, "F": 0.85, "ka": {"Oral (immediate release)": 2.0},
        "cl_per_kg": 0.27, "cl_renal": 0.0, "vss_per_kg": 0.9, "unit": "mg/L", "unit_k": 1.0,
        "factors": {"cirrhosis": 0.6, "chf": 1.0, "smoking": 1.0}, "range": {"10-20 mg/L": (10.0, 20.0)},
        "toxic": [(150.0, "liver injury risk")], "dose_steps": (250.0, 4000.0, 250.0), "ec50": 10.0,
        "give": (1000.0, "oral", 0.0), "regular": (1000.0, 6.0, "oral", 0.0), "mw": 151.2},
    "Amlodipine": {   # calcium-channel blocker: lowers blood pressure slowly, for days; a reflex rise in heart rate
        "logP": 3.0, "pKa": 8.6, "kind": "base", "fu": 0.03, "bp": 1.0, "F": 0.64, "ka": {"Oral (immediate release)": 0.12},
        "cl_per_kg": 0.42, "cl_renal": 0.0, "vss_per_kg": 21.0, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.6, "chf": 1.0, "smoking": 1.0}, "range": {"3-11 ng/mL": (3.0, 11.0)},
        "toxic": [(50.0, "severe hypotension")], "dose_steps": (2.5, 20.0, 2.5), "ec50": 5.0,
        "give": (10.0, "oral", 0.0), "regular": (10.0, 24.0, "oral", 0.0), "mw": 408.9},
    "Salbutamol": {   # beta2-agonist: opens the airways; speeds the heart
        "logP": 0.64, "pKa": 9.3, "kind": "base", "fu": 0.9, "bp": 1.0, "F": 0.5, "ka": {"Oral (immediate release)": 1.5},
        "cl_per_kg": 0.4, "cl_renal": 12.0, "vss_per_kg": 2.0, "unit": "ng/mL", "unit_k": 1000.0,
        "factors": {"cirrhosis": 0.7, "chf": 1.0, "smoking": 1.0}, "range": {"3-20 ng/mL": (3.0, 20.0)},
        "toxic": [(30.0, "tachycardia, tremor, low potassium")], "dose_steps": (0.25, 8.0, 0.25), "ec50": 5.0,
        "give": (0.25, "iv", 10 / 60), "regular": (0.25, 4.0, "iv", 10 / 60), "mw": 239.3},
})
# the first three: what the bedside gives
DRUGS["Theophylline"].update({"give": (300.0, "oral", 0.0), "regular": (300.0, 12.0, "oral", 0.0)})
DRUGS["Midazolam"].update({"give": (2.0, "iv", 2 / 60), "regular": (3.0, 1.0, "iv", 1.0), "ke0": 8.0})
DRUGS["Digoxin"].update({"give": (0.25, "oral", 0.0), "regular": (0.25, 24.0, "oral", 0.0)})
for _d in DRUGS.values():
    _d["ka"].setdefault("Oral (sustained release)", _d["ka"]["Oral (immediate release)"])

BEDSIDE_DRUGS = list(DRUGS)          # (the drug-dosing study's twelve: the anaesthetist's are added below)


def kp_rodgers_rowland(drug):
    """tissue:plasma partition coefficients (Rodgers & Rowland 2006: acids, very weak bases, neutrals)"""
    P = 10.0 ** drug["logP"]
    fu = drug["fu"]
    pka = drug["pKa"]
    acid = drug["kind"] == "acid"

    def ion(ph):
        return 1.0 + (10.0 ** (ph - pka) if acid else 10.0 ** (pka - ph)) if drug["kind"] in ("acid", "base") else 1.0
    Y = ion(PLASMA_PH)
    plasma_lip = (P * PLASMA_NL + (0.3 * P + 0.7) * PLASMA_NPL) / Y
    out = {}
    for t in TISSUES:
        nl, npl, ew, iw, ph, ar = COMP[t]
        X = ion(ph)
        # (in fat, the drug's partitioning into vegetable oil rather than octanol: log Pvo = 1.115 log P - 1.35)
        p = 10.0 ** (1.115 * drug["logP"] - 1.35) if t == "adipose" else P
        kpu = X * iw / Y + ew + (p * nl + (0.3 * p + 0.7) * npl) / Y + max(0.0, 1.0 / fu - 1.0 - plasma_lip) * ar
        out[t] = kpu * fu
    return out


# ── the virtual population ───────────────────────────────────────────────────────────────────────────────────────
POPULATIONS = {
    # name: age range, smokers' share, and the share with each condition (heart failure, cirrhosis), kidney function
    "Healthy adults": {"age": (20, 60), "smokers": 0.0},
    "Smokers": {"age": (20, 60), "smokers": 1.0},
    "Elderly (65+)": {"age": (65, 90), "smokers": 0.1},
    "Liver cirrhosis": {"age": (40, 75), "smokers": 0.2, "cirrhosis": 1.0},
    "Heart failure": {"age": (50, 85), "smokers": 0.1, "chf": 1.0},
    "Kidney impairment": {"age": (40, 80), "smokers": 0.1, "gfr": (15, 59)},
    "Mixed clinical population": {"age": (40, 85), "smokers": 0.35, "chf": 0.10, "cirrhosis": 0.05},
}
# what changes a patient's metabolism of the drug (factors on the hepatic intrinsic clearance; approximate, from the
# direction and size the literature reports: smoking induces CYP1A2, cirrhosis and heart failure lower clearance,
# clearance falls with age; Jusko et al. 1979, and the 1978 review)
SMOKING_FACTOR = 1.6
CIRRHOSIS_FACTOR = 0.40
CHF_FACTOR = 0.55
AGE_DECLINE = 0.008        # (per year over 30)
ENZYME_SD = 0.35           # (log-normal spread of CYP1A2 activity between people)


def population(n, name, rng, female_share=0.5):
    """n virtual patients: dict of arrays (sex, age, height, weight, BMI, conditions, kidney function, enzyme activity)"""
    pp = POPULATIONS[name]
    female = rng.random(n) < female_share
    age = rng.uniform(*pp["age"], n)
    height = np.where(female, rng.normal(163.0, 6.5, n), rng.normal(176.0, 7.0, n)).clip(145, 205)
    bmi = np.exp(rng.normal(math.log(25.5), 0.16, n)).clip(17, 45)
    weight = bmi * (height / 100.0) ** 2
    smoker = rng.random(n) < pp.get("smokers", 0.0)
    chf = rng.random(n) < pp.get("chf", 0.0)
    cirr = rng.random(n) < pp.get("cirrhosis", 0.0)
    if "gfr" in pp:
        gfr = rng.uniform(*pp["gfr"], n)
    else:                  # (kidney function falls with age: about 0.8 mL/min a year past 40)
        gfr = (120.0 - 0.8 * np.clip(age - 40.0, 0, None)) * np.exp(rng.normal(0, 0.15, n))
    enzyme = np.exp(rng.normal(0, ENZYME_SD, n))
    return {"female": female, "age": age, "height": height, "weight": weight, "bmi": bmi, "smoker": smoker, "chf": chf,
            "cirrhosis": cirr, "gfr": gfr, "enzyme": enzyme}


def bsa(w, h):
    """body surface area (m², Mosteller)"""
    return np.sqrt(w * h / 3600.0)


def physiology(pop):
    """each patient's organ volumes (L) and blood flows (L/h), from their sex, age, size and body fat"""
    n = len(pop["age"])
    fem = pop["female"]
    w, h, age, bmi = pop["weight"], pop["height"], pop["age"], pop["bmi"]
    ref_w = np.where(fem, REF_BODY["female"][0], REF_BODY["male"][0])
    ref_h = np.where(fem, REF_BODY["female"][1], REF_BODY["male"][1])
    s_bsa = bsa(w, h) / bsa(ref_w, ref_h)
    V, Q = {}, {}
    pick = lambda t, k: np.where(fem, REF[t][k + 1], REF[t][k])     # noqa: E731
    for t in TISSUES:
        V[t] = pick(t, 0) * s_bsa
    V["brain"] = pick("brain", 0) * np.ones(n)
    V["bone"] = pick("bone", 0) * (h / ref_h) ** 2
    # body fat (Deurenberg et al. 1991: BF% = 1.2 BMI + 0.23 age - 10.8 male - 5.4)
    bf = np.clip(1.2 * bmi + 0.23 * age - 10.8 * (~fem) - 5.4, 5, 60) / 100.0
    V["adipose"] = bf * w
    # blood (Nadler et al. 1962), a third of it arterial
    hm = h / 100.0
    vb = np.where(fem, 0.3561 * hm ** 3 + 0.03308 * w + 0.1833, 0.3669 * hm ** 3 + 0.03219 * w + 0.6041)
    V["arterial"], V["venous"] = vb * 0.30, vb * 0.70
    # muscle: what weight is left (kg ~ L), at least half its reference
    other = sum(V[t] for t in TISSUES if t != "muscle") + vb
    V["muscle"] = np.maximum(w - other, 0.5 * pick("muscle", 0))
    co = np.where(fem, REF_CO["female"], REF_CO["male"]) * s_bsa * 60.0      # L/h
    co = co * np.where(pop["chf"], 0.75, 1.0)                                # (heart failure: a lower output)
    ref_co = np.where(fem, REF_CO["female"], REF_CO["male"])
    for t in TISSUES:
        if t != "lung":
            Q[t] = pick(t, 2) / ref_co * co
    Q["lung"] = co
    return V, Q


def patients(pop, drug, route):
    """the per-patient model: volumes, flows, partition coefficients, clearances"""
    V, Q = physiology(pop)
    kp = kp_rodgers_rowland(drug)
    n = len(pop["age"])
    fu = np.full(n, drug["fu"]) * np.where(pop["cirrhosis"], 1.15, 1.0)     # (less albumin: less bound)
    bp = drug["bp"]
    # the tissue partitioning scaled once so that a reference adult's Vss is the observed one
    vss_ref = _vss(kp, {t: REF[t][0] for t in TISSUES}, sum(REF_BLOOD["male"]), bp)
    scale = drug["vss_per_kg"] * REF_BODY["male"][0] / vss_ref
    Kp = {t: kp[t] * scale * np.ones(n) for t in TISSUES}
    # the hepatic intrinsic clearance (unbound), from a reference adult's observed total clearance (well-stirred liver)
    ref_cl_h = drug["cl_per_kg"] * REF_BODY["male"][0] - drug["cl_renal"]
    q_liver_ref = (REF["liver"][2] + sum(REF[t][2] for t in PORTAL)) * 60.0
    fub = drug["fu"] / bp
    clint_ref = ref_cl_h * q_liver_ref / (fub * (q_liver_ref - ref_cl_h))       # per reference liver
    age_f = np.clip(1.0 - AGE_DECLINE * np.clip(pop["age"] - 30.0, 0, None), 0.45, 1.0)
    fx = drug.get("factors", {})
    act = pop["enzyme"] * np.where(pop["smoker"], fx.get("smoking", SMOKING_FACTOR), 1.0) \
        * np.where(pop["cirrhosis"], fx.get("cirrhosis", CIRRHOSIS_FACTOR), 1.0) * np.where(pop["chf"], fx.get("chf", CHF_FACTOR), 1.0) * age_f
    ref_liver = np.where(pop["female"], REF["liver"][1], REF["liver"][0])
    clint = clint_ref * act * V["liver"] / REF["liver"][0] * np.where(pop["female"], REF["liver"][0] / ref_liver, 1.0)
    cl_renal = drug["cl_renal"] * pop["gfr"] / 120.0
    cl_blood = drug.get("cl_blood_per_kg", 0.0) * pop["weight"]      # (esterases in the blood itself: remifentanil)
    return {"V": V, "Q": Q, "Kp": Kp, "fu": fu, "bp": bp, "clint": clint, "clr": cl_renal, "clb": cl_blood,
            "F": drug["F"], "ka": drug["ka"].get(route, 1.5), "k": drug.get("unit_k", 1.0)}


def _vss(kp, V, vblood, bp):
    return sum(kp[t] * V[t] for t in TISSUES) + vblood * bp


def system_matrix(m):
    """dx/dt = A x (x: amounts, mg; per hour): flow-limited tissues, the portal organs draining into the liver,
    metabolism in the liver (well-stirred, on the unbound drug), renal excretion, absorption from the gut's lumen"""
    V, Q, Kp, bp = m["V"], m["Q"], m["Kp"], m["bp"]
    n = len(m["clint"])
    A = np.zeros((n, NS, NS))
    a, v, ln = IX["arterial"], IX["venous"], IX["lumen"]
    q_portal = sum(Q[t] for t in PORTAL)
    q_hep = Q["liver"] + q_portal

    def flow(src, dst, rate):          # rate: per hour, on the amount in src
        A[:, dst, src] += rate
        A[:, src, src] -= rate
    # lungs: all venous blood through them to the arteries
    flow(v, IX["lung"], Q["lung"] / V["venous"])
    flow(IX["lung"], a, Q["lung"] * bp / (Kp["lung"] * V["lung"]))
    for t in TISSUES:
        if t == "lung":
            continue
        i = IX[t]
        flow(a, i, Q[t] / V["arterial"])                                # (arterial blood in)
        out = (q_hep if t == "liver" else Q[t]) * bp / (Kp[t] * V[t])   # (venous blood out, at equilibrium with the tissue)
        flow(i, IX["liver"] if t in PORTAL else v, out)
    # elimination: the liver (unbound concentration there: C * fu / Kp), the kidneys (filtered from the plasma there)
    A[:, IX["liver"], IX["liver"]] -= m["clint"] * m["fu"] / (Kp["liver"] * V["liver"])
    A[:, IX["kidney"], IX["kidney"]] -= m["clr"] / (Kp["kidney"] * V["kidney"])
    # broken down in the blood itself (esterases), at the same rate wherever the blood is
    kb = np.asarray(m.get("clb", 0.0)) / (V["arterial"] + V["venous"])
    A[:, a, a] -= kb
    A[:, v, v] -= kb
    # absorption: the lumen into the gut wall
    flow(ln, IX["gut"], m["ka"] * np.ones(n))
    return A


def discretise(A, dt, torch, dev):
    """exact steps of dt hours: x(t + dt) = Phi x(t) + Gamma r, r a constant infusion (mg/h) into the venous blood"""
    n = A.shape[0]
    M = np.zeros((n, NS + 1, NS + 1))
    M[:, :NS, :NS] = A
    M[:, IX["venous"], NS] = 1.0
    E = torch.linalg.matrix_exp(torch.as_tensor(M * dt, dtype=torch.float64, device=dev))
    return E[:, :NS, :NS].contiguous(), E[:, :NS, NS].contiguous()


def regimen(route, dose_mg, interval_h, days, loading_mg=0.0, infusion_h=0.5):
    """the doses: [(time h, amount mg, kind)] with kind oral or iv (an infusion over infusion_h)"""
    kind = "iv" if route.startswith("IV") else "oral"
    out = []
    t, k = 0.0, 0
    while t < days * 24.0 - 1e-9:
        amt = dose_mg + (loading_mg if k == 0 else 0.0)
        out.append((t, amt, kind))
        t += interval_h
        k += 1
    return out, infusion_h


def simulate(m, doses, days, torch, dev, dt=0.1, F=None, infusion_h=0.5, per_kg=None, keep_tissues=False):
    """the course of every patient: plasma concentration (n, T), and the tissue concentrations if asked (n, T, tissues)
    doses: [(time, mg, kind)]; per_kg: each patient's weight when the doses are mg/kg"""
    A = system_matrix(m)
    Phi, Gam = discretise(A, dt, torch, dev)
    n = A.shape[0]
    T = int(round(days * 24.0 / dt)) + 1
    x = torch.zeros(n, NS, dtype=torch.float64, device=dev)
    scale = torch.ones(n, dtype=torch.float64, device=dev) if per_kg is None else torch.as_tensor(per_kg, dtype=torch.float64, device=dev)
    Fv = float(m["F"] if F is None else F)
    vv = torch.as_tensor(m["V"]["venous"], dtype=torch.float64, device=dev)
    Vt = torch.as_tensor(np.stack([m["V"][t] for t in TISSUES], 1), dtype=torch.float64, device=dev)
    plasma = torch.zeros(n, T, dtype=torch.float32, device=dev)
    tissue = torch.zeros(n, T, len(TISSUES), dtype=torch.float32, device=dev) if keep_tissues else None
    by_step = {}
    for t0, amt, kind in doses:
        by_step.setdefault(int(round(t0 / dt)), []).append((amt, kind))
    inf_steps = max(1, int(round(infusion_h / dt)))
    rate = torch.zeros(n, dtype=torch.float64, device=dev)
    stops = {}
    for k in range(T):
        for amt, kind in by_step.get(k, ()):
            if kind == "oral":
                x[:, IX["lumen"]] += amt * Fv * scale
            else:
                r = amt * scale / (inf_steps * dt)
                rate = rate + r
                stops.setdefault(k + inf_steps, []).append(r)
        for r in stops.pop(k, ()):
            rate = rate - r
        plasma[:, k] = (x[:, IX["venous"]] / vv / m["bp"] * m["k"]).float()
        if keep_tissues:
            tissue[:, k] = (x[:, :len(TISSUES)] / Vt * m["k"]).float()
        x = torch.bmm(Phi, x[:, :, None])[:, :, 0] + Gam * rate[:, None]
    return plasma, tissue


def outcomes(plasma, dt, interval_h, lo, hi, last_h=24.0):
    """per patient, over the last day (steady state): trough, peak, average, the share of time in range, and the
    highest concentration over the whole course"""
    k = max(1, int(round(last_h / dt)))
    ss = plasma[:, -k:]
    cmin, cmax, cavg = ss.min(1).values, ss.max(1).values, ss.mean(1)
    inr = ((ss >= lo) & (ss <= hi)).float().mean(1)
    peak = plasma.max(1).values
    return {"cmin": cmin, "cmax": cmax, "cavg": cavg, "in_range": inr, "peak": peak}


# ── the patient for the 3D view ──────────────────────────────────────────────────────────────────────────────────
# A mockup until the anatomy models come: a translucent body and, inside it, one part per compartment of the model,
# each a mesh of its own with a material named for it, so the 3D view colours it by the drug's concentration there.
# Shapes are signed distance fields (ellipsoids, tapered capsules, swept tubes, joined smoothly), turned into meshes by
# marching cubes. Metres, z up, the feet at 0, facing +y, a 1.76 m man.
# ── signed distances (negative inside) ──────────────────────────────────────────────────────────────────────────
def sd_ellipsoid(P, c, r):
    q = (P - np.asarray(c)) / np.asarray(r)
    k = np.linalg.norm(q, axis=-1)
    return (k - 1.0) * min(r)


def sd_capsule(P, a, b, ra, rb=None):
    """a capsule from a to b, its radius going from ra to rb"""
    rb = ra if rb is None else rb
    a, b = np.asarray(a, float), np.asarray(b, float)
    ab = b - a
    t = np.clip(((P - a) @ ab) / (ab @ ab), 0, 1)
    r = ra + (rb - ra) * t
    return np.linalg.norm(P - (a + t[..., None] * ab), axis=-1) - r


def sd_tube(P, pts, r):
    """a tube along a polyline"""
    d = np.full(P.shape[:-1], 1e9)
    for a, b in zip(pts[:-1], pts[1:]):
        d = np.minimum(d, sd_capsule(P, a, b, r))
    return d


def smin(a, b, k):
    """a smooth union (k: how far the blend reaches, m)"""
    h = np.clip(0.5 + 0.5 * (b - a) / k, 0, 1)
    return b * (1 - h) + a * h - k * h * (1 - h)


def smax(a, b, k):
    return -smin(-a, -b, k)


# ── the body ────────────────────────────────────────────────────────────────────────────────────────────────────
def body_sd(P):
    """the outside of the skin"""
    d = sd_ellipsoid(P, (0, 0.005, 1.635), (0.078, 0.097, 0.115))                    # head
    d = smin(d, sd_capsule(P, (0, -0.005, 1.47), (0, 0.0, 1.56), 0.052), 0.03)        # neck
    d = smin(d, sd_capsule(P, (-0.115, 0, 1.36), (0.115, 0, 1.36), 0.105, 0.105), 0.06)   # chest
    d = smin(d, sd_ellipsoid(P, (0, 0.0, 1.22), (0.155, 0.105, 0.17)), 0.06)           # ribcage, belly
    d = smin(d, sd_ellipsoid(P, (0, -0.005, 1.0), (0.165, 0.11, 0.13)), 0.07)           # hips
    d = smin(d, sd_capsule(P, (-0.17, 0, 1.43), (0.17, 0, 1.43), 0.06), 0.06)           # shoulders
    for s in (-1, 1):
        d = smin(d, sd_capsule(P, (s * 0.19, 0, 1.42), (s * 0.24, 0.0, 1.13), 0.048, 0.04), 0.03)     # upper arm
        d = smin(d, sd_capsule(P, (s * 0.24, 0.0, 1.13), (s * 0.28, 0.03, 0.86), 0.039, 0.03), 0.02)  # forearm
        d = smin(d, sd_ellipsoid(P, (s * 0.29, 0.035, 0.80), (0.022, 0.045, 0.075)), 0.02)            # hand
        d = smin(d, sd_capsule(P, (s * 0.085, 0, 0.93), (s * 0.10, 0.0, 0.50), 0.085, 0.055), 0.05)   # thigh
        d = smin(d, sd_capsule(P, (s * 0.10, 0.0, 0.50), (s * 0.11, -0.01, 0.09), 0.052, 0.035), 0.03)  # shin
        d = smin(d, sd_capsule(P, (s * 0.11, -0.02, 0.045), (s * 0.12, 0.13, 0.03), 0.035, 0.03), 0.02)  # foot
    return d


def gut_path():
    """the small intestine: a tube coiled in the belly, back and forth in rows"""
    pts = []
    rows = 7
    for k in range(rows):
        z = 1.07 - k * 0.022
        xs = np.linspace(-0.085, 0.085, 14)
        if k % 2:
            xs = xs[::-1]
        for j, x in enumerate(xs):
            pts.append((x, 0.045 + 0.018 * math.sin(j * 1.3 + k), z + 0.008 * math.sin(j * 2.1)))
    return np.array(pts)


def colon_path():
    """the large intestine: up the right side, across, down the left"""
    return np.array([(-0.115, 0.03, 0.93), (-0.12, 0.035, 1.03), (-0.10, 0.05, 1.11), (0.0, 0.07, 1.12), (0.10, 0.05, 1.115),
                     (0.12, 0.03, 1.04), (0.11, 0.025, 0.95), (0.05, 0.0, 0.92)])


def skeleton_sd(P):
    d = sd_ellipsoid(P, (0, 0.0, 1.645), (0.068, 0.085, 0.09))                         # skull
    d = smin(d, sd_tube(P, [(0, -0.045, 0.98), (0, -0.05, 1.15), (0, -0.04, 1.35), (0, -0.02, 1.56)], 0.017), 0.01)   # spine
    for k in range(9):                                                                  # ribs
        z = 1.40 - k * 0.028
        w = 0.10 + 0.035 * math.sin(math.pi * (k + 1) / 10)
        ang = np.linspace(0.25 * math.pi, 1.75 * math.pi, 18)
        ring = [(w * math.sin(a), -0.05 + 0.095 * (1 - math.cos(a)) / 2 * 1.4 - 0.0, z - 0.03 * math.sin(a / 2)) for a in ang]
        d = np.minimum(d, sd_tube(P, ring, 0.006))
    d = smin(d, sd_capsule(P, (0, 0.06, 1.38), (0, 0.075, 1.24), 0.012), 0.005)        # sternum
    d = smin(d, sd_ellipsoid(P, (0, -0.01, 0.97), (0.13, 0.06, 0.06)), 0.02)            # pelvis
    for s in (-1, 1):
        d = np.minimum(d, sd_capsule(P, (s * 0.02, 0.02, 1.47), (s * 0.17, -0.01, 1.45), 0.009))        # collarbone
        d = np.minimum(d, sd_capsule(P, (s * 0.19, 0, 1.40), (s * 0.24, 0.0, 1.14), 0.016, 0.012))     # humerus
        d = np.minimum(d, sd_capsule(P, (s * 0.24, 0.0, 1.12), (s * 0.28, 0.03, 0.87), 0.011, 0.009))  # forearm
        d = np.minimum(d, sd_capsule(P, (s * 0.09, 0.0, 0.92), (s * 0.10, 0.0, 0.52), 0.022, 0.018))   # femur
        d = np.minimum(d, sd_capsule(P, (s * 0.10, 0.0, 0.48), (s * 0.11, -0.01, 0.08), 0.017, 0.013))  # tibia
    return d


ORGANS = {
    # name: (signed distance function, the node's centre (heart and lungs scale about it))
    "brain": (lambda P: smax(sd_ellipsoid(P, (0, 0.005, 1.65), (0.062, 0.078, 0.062)), -sd_ellipsoid(P, (0, 0.005, 1.585), (0.07, 0.09, 0.03)), 0.01), (0, 0.005, 1.65)),
    "heart": (lambda P: smin(sd_ellipsoid(P, (0.022, 0.05, 1.31), (0.05, 0.042, 0.062)), sd_capsule(P, (0.01, 0.04, 1.35), (0.0, 0.02, 1.42), 0.014), 0.01), (0.022, 0.05, 1.31)),
    "lung": (lambda P: smax(np.minimum(sd_ellipsoid(P, (-0.075, 0.0, 1.31), (0.062, 0.07, 0.125)), sd_ellipsoid(P, (0.08, -0.005, 1.31), (0.058, 0.07, 0.12))),
                            -sd_ellipsoid(P, (0.022, 0.05, 1.31), (0.06, 0.055, 0.075)), 0.01), (0, 0, 1.31)),
    "liver": (lambda P: smax(sd_ellipsoid(P, (-0.05, 0.03, 1.165), (0.115, 0.075, 0.06)), sd_capsule(P, (-0.2, 0.03, 1.1), (0.08, 0.03, 1.19), 0.08), 0.02), (-0.05, 0.03, 1.165)),
    "gut": (lambda P: smin(np.minimum(sd_tube(P, gut_path(), 0.013), sd_tube(P, colon_path(), 0.022)),
                           smin(sd_ellipsoid(P, (0.07, 0.045, 1.14), (0.055, 0.035, 0.05)), sd_capsule(P, (0.035, 0.03, 1.19), (0.07, 0.045, 1.15), 0.016), 0.02), 0.005), (0, 0.04, 1.03)),
    "spleen": (lambda P: sd_ellipsoid(P, (0.115, -0.04, 1.15), (0.028, 0.04, 0.055)), (0.115, -0.04, 1.15)),
    "pancreas": (lambda P: sd_capsule(P, (-0.03, 0.0, 1.10), (0.09, -0.01, 1.12), 0.016, 0.01), (0.03, 0.0, 1.11)),
    "kidney": (lambda P: np.minimum(sd_ellipsoid(P, (-0.065, -0.055, 1.06), (0.025, 0.02, 0.05)), sd_ellipsoid(P, (0.065, -0.055, 1.08), (0.025, 0.02, 0.05))), (0, -0.055, 1.07)),
    "bone": (skeleton_sd, (0, 0, 0)),
}
# the blood: the aorta (arterial) from the heart down to the legs and up to the head; the venae cavae (venous) beside it
ARTERIAL = [[(0.02, 0.03, 1.36), (0.03, 0.02, 1.42), (0.0, -0.02, 1.40), (-0.01, -0.03, 1.25), (-0.005, -0.03, 1.0), (0.0, -0.02, 0.95)],
            [(0.0, -0.02, 0.95), (0.08, 0.0, 0.88), (0.10, 0.01, 0.50), (0.11, 0.0, 0.12)],
            [(0.0, -0.02, 0.95), (-0.08, 0.0, 0.88), (-0.10, 0.01, 0.50), (-0.11, 0.0, 0.12)],
            [(0.02, 0.0, 1.42), (0.025, 0.0, 1.55), (0.03, 0.0, 1.60)], [(0.0, 0.0, 1.42), (-0.025, 0.0, 1.55), (-0.03, 0.0, 1.60)],
            [(0.03, 0.01, 1.43), (0.17, 0.0, 1.42), (0.24, 0.01, 1.13), (0.28, 0.035, 0.85)],
            [(-0.01, 0.01, 1.43), (-0.17, 0.0, 1.42), (-0.24, 0.01, 1.13), (-0.28, 0.035, 0.85)]]
VENOUS = [[(0.0, 0.04, 1.33), (-0.03, -0.005, 1.25), (-0.03, -0.015, 1.0), (-0.02, -0.01, 0.95)],
          [(-0.02, -0.01, 0.95), (-0.07, 0.015, 0.88), (-0.09, 0.025, 0.50), (-0.10, 0.015, 0.12)],
          [(-0.02, -0.01, 0.95), (0.07, 0.015, 0.88), (0.09, 0.025, 0.50), (0.10, 0.015, 0.12)],
          [(0.0, 0.04, 1.36), (-0.03, 0.02, 1.45), (-0.04, 0.01, 1.60)]]


def mesh(sd, lo, hi, step):
    """a closed surface where sd is 0 inside the box lo-hi, cells of step (m): vertices, faces, normals"""
    from skimage.measure import marching_cubes
    lo, hi = np.asarray(lo, float) - 2 * step, np.asarray(hi, float) + 2 * step
    axes = [np.arange(lo[i], hi[i] + step, step) for i in range(3)]
    G = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    vol = sd(G.reshape(-1, 3)).reshape(G.shape[:3]).astype(np.float32)
    v, f, n, _ = marching_cubes(vol, 0.0, spacing=(step, step, step), gradient_direction="ascent")
    v = v + lo
    return v, f[:, ::-1].copy(), -n          # (outward, counter-clockwise)


def bbox_of(sd, lo, hi, step=0.02):
    axes = [np.arange(lo[i], hi[i] + step, step) for i in range(3)]
    G = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    ins = G[sd(G) < step]
    return ins.min(0) - step, ins.max(0) + step


def parts(step=0.004):
    """[(name, centre, vertices relative to it, faces, normals)] for every part of the patient"""
    out = []
    lo, hi = (-0.36, -0.16, 0.0), (0.36, 0.2, 1.78)
    v, f, n = mesh(body_sd, lo, hi, 0.006)
    out.append(("skin", np.zeros(3), v, f, n))
    for name, (sd, c) in ORGANS.items():
        a, b = bbox_of(sd, lo, hi)
        v, f, n = mesh(sd, a, b, step if name != "bone" else 0.004)
        out.append((name, np.asarray(c, float), v - np.asarray(c, float), f, n))
    for name, paths, r in (("arterial", ARTERIAL, 0.0055), ("venous", VENOUS, 0.0065)):
        sd = lambda P, paths=paths, r=r: np.min([sd_tube(P, p, r) for p in paths], axis=0)   # noqa: E731
        a, b = bbox_of(sd, lo, hi, 0.01)
        v, f, n = mesh(sd, a, b, 0.0035)
        out.append((name, np.zeros(3), v, f, n))
    # muscle and fat: the limbs' and the trunk's flesh under the skin, as shells (fat just under it, muscle deeper)
    inner = lambda P: body_sd(P) + 0.012                                                 # noqa: E731
    cavity = lambda P: smin(sd_ellipsoid(P, (0, 0.0, 1.2), (0.13, 0.085, 0.22)), sd_ellipsoid(P, (0, 0.0, 1.64), (0.075, 0.09, 0.1)), 0.04)   # noqa: E731
    v, f, n = mesh(lambda P: smax(inner(P), -cavity(P), 0.02), lo, hi, 0.008)
    out.append(("muscle", np.zeros(3), v, f, n))
    v, f, n = mesh(lambda P: smax(body_sd(P) + 0.004, -(body_sd(P) + 0.011), 0.002), lo, hi, 0.008)
    out.append(("adipose", np.zeros(3), v, f, n))
    return out


# ── the patient from a CT case: every labelled structure, grouped by the model's compartments ────────────────────────
# A CT scan and its experts' labels (the TotalSegmentator dataset: a mask per structure, segmentations/<name>.nii.gz)
# become the patient's parts: each compartment's structures joined into one surface (the liver and gallbladder; the
# stomach, bowel and oesophagus as the gut; every bone as the skeleton...), the skin and the fat under it from the CT's
# own densities. Metres, z up, facing +y, the scan's middle at TORSO_Z above the platform (it floats, as on a display).
TORSO_Z = 1.15
COMPARTMENT_OF = [   # (compartment, words in the structures' names), the first that matches
    ("venous", ("vena_cava", "portal_vein", "brachiocephalic_vein", "iliac_vena", "subclavian_vein", "jugular")),
    ("arterial", ("aorta", "artery", "brachiocephalic_trunk", "pulmonary_vein", "atrial_appendage")),
    ("heart", ("heart", "myocardium", "atrium", "ventricle")),
    ("lung", ("lung", "trachea", "bronch")),
    ("liver", ("liver", "gallbladder")),
    ("gut", ("stomach", "small_bowel", "duodenum", "colon", "esophagus", "rectum")),
    ("spleen", ("spleen",)),
    ("pancreas", ("pancreas",)),
    ("kidney", ("kidney", "adrenal")),
    ("brain", ("brain",)),
    ("ribcage", ("rib_", "sternum", "costal")),      # (a part of its own, so that it can move with the breath; drawn as bone)
    ("bone", ("vertebra", "hip_", "sacrum", "femur", "humerus", "scapula", "clavicul", "skull")),
    ("muscle", ("gluteus", "iliopsoas", "autochthon", "muscle")),
    ("rest", ("bladder", "prostate", "thyroid", "spinal_cord")),
]
MATERIAL_OF = {"ribcage": "bone"}   # (parts drawn with another part's material)
PART_TRIANGLES = {"skin": 70000, "adipose": 50000, "bone": 90000, "ribcage": 60000, "gut": 60000, "lung": 40000, "muscle": 40000,
                  "arterial": 25000, "venous": 25000, "liver": 16000, "heart": 12000, "kidney": 8000, "spleen": 6000,
                  "pancreas": 6000, "brain": 12000, "rest": 8000}


def compartment(name):
    n = name.lower()
    for comp, words in COMPARTMENT_OF:
        if any(w in n for w in words):
            return comp
    return None


def ct_parts(files, log):
    """the patient's parts from a CT case: [(name, centre, vertices relative to it, faces, normals)]"""
    import SimpleITK as sitk
    from scipy import ndimage
    ct = [f for f in files if str(f.name).lower() in ("ct.nii.gz", "ct.nii")]
    masks = [f for f in files if str(f.name).lower().endswith((".nii.gz", ".nii")) and f not in ct]
    if not ct or not masks:
        raise RuntimeError("A CT case is needed: its scan (ct.nii.gz) and its labels (segmentations/<structure>.nii.gz).")
    img = sitk.ReadImage(str(ct[0].path))
    hu = sitk.GetArrayFromImage(img).astype(np.float32)                   # (z, y, x)
    groups, used = {}, {}
    for f in masks:
        nm = str(f.name).lower().replace(".nii.gz", "").replace(".nii", "")
        comp = compartment(nm)
        if comp is None:
            continue
        m = sitk.GetArrayFromImage(sitk.ReadImage(str(f.path))) > 0
        if m.shape != hu.shape or not m.any():
            continue
        groups[comp] = m if comp not in groups else (groups[comp] | m)
        used.setdefault(comp, []).append(nm)
    for comp in sorted(used):
        log("  %s: %s" % (comp, ", ".join(sorted(used[comp]))))
    # the body (denser than -300 HU, its largest piece, its holes filled) and the fat under its skin (-190 to -30 HU)
    body = hu > -300
    lab, n = ndimage.label(body)
    if n > 1:
        body = lab == 1 + int(np.argmax(ndimage.sum(body, lab, range(1, n + 1))))
    body = ndimage.binary_fill_holes(ndimage.binary_closing(body, iterations=2))
    inner = ndimage.binary_erosion(body, iterations=12)
    fat = body & ~inner & (hu > -190) & (hu < -30)
    fat = ndimage.binary_opening(fat, iterations=1)
    groups["skin"], groups["adipose"] = body, fat
    sp = np.array(img.GetSpacing(), float)
    d = np.array(img.GetDirection(), float).reshape(3, 3)
    org = np.array(img.GetOrigin(), float)
    zz = np.nonzero(body.any((1, 2)))[0]
    mid = (org + (np.array([0, 0, (zz[0] + zz[-1]) / 2.0]) * sp) @ d.T)[2]
    out = []
    for comp, mask in groups.items():
        mesh = _surface(mask, sp, d, org, PART_TRIANGLES.get(comp, 10000), sigma=1.6 if comp in ("skin", "adipose") else 0.9)
        if mesh is None:
            continue
        v = np.asarray(mesh.vertices)
        # LPS (mm) -> the cell (m): the patient's left to -x, the front to +y, up to +z, the scan's middle at TORSO_Z
        v = np.column_stack([-v[:, 0], -v[:, 1], v[:, 2] - mid]) / 1000.0 + np.array([0, 0, TORSO_Z])
        c = v.mean(0)
        if comp == "ribcage":                # (its pivot at the back, by the spine: it widens and rises forwards)
            c[1] = v[:, 1].min()
        nrm = np.asarray(mesh.vertex_normals) * np.array([-1, -1, 1])
        out.append((comp, c, v - c, np.asarray(mesh.faces), nrm))
        log("  %s: %d triangles." % (comp, len(mesh.faces)))
    return out


def _surface(mask, sp, d, org, tri, sigma=0.9):
    import trimesh
    from scipy import ndimage
    from skimage.measure import marching_cubes
    if mask.sum() < 50:
        return None
    sl = ndimage.find_objects(mask.astype(np.uint8))[0]
    sub = ndimage.gaussian_filter(np.pad(mask[sl], 2).astype(np.float32), sigma)
    try:
        v, f, _, _ = marching_cubes(sub, 0.5, allow_degenerate=False)
    except (ValueError, RuntimeError):
        return None
    ijk = (v + np.array([s.start for s in sl], float) - 2.0)[:, ::-1]
    lps = org + (ijk * sp) @ d.T
    mesh = trimesh.Trimesh(lps, f, process=True)
    trimesh.smoothing.filter_taubin(mesh, lamb=0.5, nu=0.53, iterations=12)
    if len(mesh.faces) > tri:
        try:
            mesh = mesh.simplify_quadric_decimation(face_count=tri)
        except Exception:      # noqa: BLE001
            pass
    trimesh.repair.fix_normals(mesh)
    return mesh


# ── a glTF binary with one node per part, each its own material ─────────────────────────────────────────────────
PART_STYLE = {   # name: (base colour, roughness, metalness, opacity)
    "skin": ((0.55, 0.75, 0.95), 0.15, 0.0, 0.16),
    "adipose": ((0.95, 0.85, 0.55), 0.5, 0.0, 0.10),
    "muscle": ((0.75, 0.30, 0.32), 0.5, 0.0, 0.22),
    "bone": ((0.92, 0.90, 0.84), 0.55, 0.0, 0.55),
    "arterial": ((0.80, 0.12, 0.15), 0.35, 0.0, 1.0),
    "venous": ((0.25, 0.30, 0.75), 0.35, 0.0, 1.0),
}


def patient_glb(ps=None, step=0.004):
    ps = parts(step) if ps is None else ps
    blobs, views, accs, meshes, nodes, mats = [], [], [], [], [], []
    off = 0

    def blob(b, target):
        nonlocal off
        views.append({"buffer": 0, "byteOffset": off, "byteLength": len(b), "target": target})
        blobs.append(b + b"\0" * (-len(b) % 4))
        off += len(blobs[-1])
        return len(views) - 1
    for k, (name, c, v, f, n) in enumerate(ps):
        v, n = v.astype(np.float32), (n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-9)).astype(np.float32)
        idx = f.astype(np.uint32).reshape(-1)
        a0 = len(accs)
        accs.append({"bufferView": blob(v.tobytes(), 34962), "componentType": 5126, "count": len(v), "type": "VEC3",
                     "min": v.min(0).tolist(), "max": v.max(0).tolist()})
        accs.append({"bufferView": blob(n.tobytes(), 34962), "componentType": 5126, "count": len(n), "type": "VEC3"})
        accs.append({"bufferView": blob(idx.tobytes(), 34963), "componentType": 5125, "count": len(idx), "type": "SCALAR"})
        mname = MATERIAL_OF.get(name, name)
        col, rough, metal, op = PART_STYLE.get(mname, ((0.8, 0.8, 0.8), 0.45, 0.0, 1.0))
        mt = {"name": mname, "pbrMetallicRoughness": {"baseColorFactor": list(col) + [op], "roughnessFactor": rough, "metallicFactor": metal},
              "doubleSided": op < 1}
        if op < 1:
            mt["alphaMode"] = "BLEND"
        mats.append(mt)
        meshes.append({"name": name, "primitives": [{"attributes": {"POSITION": a0, "NORMAL": a0 + 1}, "indices": a0 + 2, "material": k}]})
        nodes.append({"name": name, "mesh": k, "translation": [float(c[0]), float(c[1]), float(c[2])]})
    # (z up, as the cell's other files: the scene says so with "up": "z")
    js = {"asset": {"version": "2.0", "generator": "EsaSim virtual patient"}, "scene": 0,
          "scenes": [{"nodes": [len(nodes)]}], "nodes": nodes + [{"name": "patient", "children": list(range(len(nodes)))}],
          "meshes": meshes, "materials": mats, "bufferViews": views, "accessors": accs, "buffers": [{"byteLength": off}]}
    j = json.dumps(js, separators=(",", ":")).encode()
    j += b" " * (-len(j) % 4)
    binc = b"".join(blobs)
    return (struct.pack("<III", 0x46546C67, 2, 28 + len(j) + len(binc)) + struct.pack("<I", len(j)) + b"JSON" + j
            + struct.pack("<I", len(binc)) + b"BIN\0" + binc), {p[0]: len(p[3]) for p in ps}


# ── the block: modes ─────────────────────────────────────────────────────────────────────────────────────────────────
TARGETS = list(DRUGS["Theophylline"]["range"])
ROUTES = ["Oral (immediate release)", "Oral (sustained release)", "IV infusion"]
ORGAN_VIEW = ["brain", "heart", "lung", "liver", "gut", "spleen", "pancreas", "kidney", "muscle", "adipose", "bone", "skin"]
DT = 0.1                   # (h) the simulation's step (exact for the linear model: any step gives the same values at it)


def device(ctx):
    import torch
    want = str(ctx.params.get("accelerator") or "auto").lower()
    if want.startswith("cpu") or not torch.cuda.is_available():
        return "cpu", "the CPU"
    return "cuda:0", torch.cuda.get_device_name(0)


def settings(ctx):
    p = ctx.params
    drug = DRUGS[str(p.get("drug") or "Theophylline")]
    target = str(p.get("target") or TARGETS[0])
    lo, hi = drug["range"].get(target, next(iter(drug["range"].values())))
    if target not in drug["range"]:
        target = next(iter(drug["range"]))
    route = str(p.get("route") or ROUTES[0])
    return {
        "drug": drug, "drug_name": str(p.get("drug") or "Theophylline"),
        "population": str(p.get("population") or "Healthy adults"),
        "n": int(p.get("patients") or 2000),
        "dose": float(p.get("dose_mg") or 300.0),
        "per_kg": str(p.get("dose_basis") or "mg").lower().startswith("mg/kg"),
        "interval": float(p.get("interval_h") or 12.0),
        "route": route,
        "loading": float(p.get("loading_mg") or 0.0),
        "infusion_h": float(p.get("infusion_h") or 0.5),
        "days": float(p.get("days") or 5.0),
        "lo": lo, "hi": hi, "target": target,
        "female": float(p.get("female_pct") if p.get("female_pct") not in (None, "") else 50) / 100.0,
        "seed": int(p.get("seed") or 1),
    }


def describe(s):
    unit = "mg/kg" if s["per_kg"] else "mg"
    how = {"Oral (immediate release)": "by mouth", "Oral (sustained release)": "by mouth, sustained release",
           "IV infusion": "IV over %.1f h" % s["infusion_h"]}[s["route"]]
    load = (", %g %s more with the first dose" % (s["loading"], unit)) if s["loading"] else ""
    return "%s %g %s every %g h %s%s, for %g days" % (s["drug_name"], s["dose"], unit, s["interval"], how, load, s["days"])


def unit_course(s, pop, torch, dev, keep_tissues=False):
    """the course of every patient for one unit of dose (1 mg, or 1 mg/kg) of the regimen: the model is linear in the
    dose, so any dose's course is this one scaled (the loading dose a unit dose of its own at time 0)"""
    m = patients(pop, s["drug"], s["route"])
    doses, inf_h = regimen(s["route"], 1.0, s["interval"], s["days"])
    w = pop["weight"] if s["per_kg"] else None
    plasma, tissue = simulate(m, doses, s["days"], torch, dev, dt=DT, infusion_h=s["infusion_h"], per_kg=w, keep_tissues=keep_tissues)
    first = (None, None)
    if s["loading"]:
        first = simulate(m, [(0.0, 1.0, doses[0][2])], s["days"], torch, dev, dt=DT, infusion_h=s["infusion_h"], per_kg=w, keep_tissues=keep_tissues)
    return m, (plasma, tissue), first


def scaled(s, unit, first):
    """the plasma (and tissue) courses at the regimen's dose"""
    pl = unit[0] * s["dose"]
    ts = unit[1] * s["dose"] if unit[1] is not None else None
    if s["loading"]:
        pl = pl + first[0] * s["loading"]
        if ts is not None:
            ts = ts + first[1] * s["loading"]
    return pl, ts


def clearance(m):
    """each patient's total clearance (L/h): the liver's (well-stirred) and the kidneys'"""
    q = (m["Q"]["liver"] + sum(m["Q"][t] for t in PORTAL))
    fub = m["fu"] / m["bp"]
    return q * fub * m["clint"] / (q + fub * m["clint"]) + m["clr"] + np.asarray(m.get("clb", 0.0))


def condition(pop, i):
    c = [nm for k, nm in (("cirrhosis", "liver cirrhosis"), ("chf", "heart failure")) if pop[k][i]]
    if pop["gfr"][i] < 60:
        c.append("kidney impairment")
    return ", ".join(c) if c else "none"


def trial_run(ctx):
    import csv
    import torch
    s = settings(ctx)
    dev, label = device(ctx)
    rng = np.random.default_rng(s["seed"])
    pop = population(s["n"], s["population"], rng, s["female"])
    ctx.log("A virtual trial of %d patients (%s) on %s: %s." % (s["n"], s["population"].lower(), label, describe(s)))
    ctx.progress(10, "simulating %d patients" % s["n"])
    m, unit, first = unit_course(s, pop, torch, dev)
    pl, _ = scaled(s, unit, first)
    o = outcomes(pl, DT, s["interval"], s["lo"], s["hi"])
    o = {k: v.cpu().numpy() for k, v in o.items()}
    cl = clearance(m)
    toxic = s["drug"]["toxic"]
    T = pl.shape[1]
    t = np.arange(T) * DT
    # the course: percentiles of the plasma concentration over time
    q = torch.quantile(pl[:, ::5].float(), torch.tensor([0.05, 0.25, 0.5, 0.75, 0.95], device=pl.device), dim=0).cpu().numpy()
    with open(ctx.output("course", "course.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time_h", "p05", "p25", "median", "p75", "p95", "target_low", "target_high", "toxic"])
        for k in range(q.shape[1]):
            w.writerow([round(t[k * 5], 2)] + [round(float(x), 3) for x in q[:, k]] + [s["lo"], s["hi"], toxic[0][0]])
    with open(ctx.output("patients", "patients.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["patient", "sex", "age", "weight_kg", "height_cm", "bmi", "smoker", "condition", "gfr_ml_min", "clearance_l_h",
                    "half_life_h", "trough_mg_l", "peak_mg_l", "average_mg_l", "time_in_range_pct", "highest_mg_l", "outcome"])
        vss = s["drug"]["vss_per_kg"] * pop["weight"]
        for i in range(s["n"]):
            out = "in range" if s["lo"] <= o["cavg"][i] <= s["hi"] else ("too low" if o["cavg"][i] < s["lo"] else "too high")
            if o["peak"][i] > toxic[0][0]:
                out += "; " + [d for c, d in toxic if o["peak"][i] > c][-1]
            w.writerow([i + 1, "F" if pop["female"][i] else "M", round(float(pop["age"][i])), round(float(pop["weight"][i]), 1),
                        round(float(pop["height"][i])), round(float(pop["bmi"][i]), 1), "yes" if pop["smoker"][i] else "no", condition(pop, i),
                        round(float(pop["gfr"][i])), round(float(cl[i]), 2), round(float(math.log(2) * vss[i] / cl[i]), 1),
                        round(float(o["cmin"][i]), 2), round(float(o["cmax"][i]), 2), round(float(o["cavg"][i]), 2),
                        round(float(o["in_range"][i]) * 100, 1), round(float(o["peak"][i]), 2), out])
    # by group: who is under-, well- and over-exposed
    groups = {"all": np.ones(s["n"], bool), "smokers": pop["smoker"], "non-smokers": ~pop["smoker"], "under 40": pop["age"] < 40,
              "40-64": (pop["age"] >= 40) & (pop["age"] < 65), "65 and over": pop["age"] >= 65, "women": pop["female"], "men": ~pop["female"],
              "heart failure": pop["chf"], "liver cirrhosis": pop["cirrhosis"], "kidney impairment (GFR < 60)": pop["gfr"] < 60}
    rows = []
    for g, sel in groups.items():
        if sel.sum() < 5:
            continue
        rows.append({"group": g, "patients": int(sel.sum()),
                     "in_range_pct": round(100 * float(((o["cavg"] >= s["lo"]) & (o["cavg"] <= s["hi"]))[sel].mean()), 1),
                     "too_low_pct": round(100 * float((o["cavg"] < s["lo"])[sel].mean()), 1),
                     "too_high_pct": round(100 * float((o["cavg"] > s["hi"])[sel].mean()), 1),
                     "toxic_peak_pct": round(100 * float((o["peak"] > toxic[0][0])[sel].mean()), 1),
                     "median_average_mg_l": round(float(np.median(o["cavg"][sel])), 2)})
    with open(ctx.output("groups", "groups.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    a = rows[0]
    summary = {"patients": s["n"], "in_range_pct": a["in_range_pct"], "too_low_pct": a["too_low_pct"], "too_high_pct": a["too_high_pct"],
               "toxic_peak_pct": a["toxic_peak_pct"], "median_average_mg_l": a["median_average_mg_l"],
               "median_clearance_l_h": round(float(np.median(cl)), 2),
               "median_half_life_h": round(float(np.median(math.log(2) * s["drug"]["vss_per_kg"] * pop["weight"] / cl)), 1)}
    # the same results live, for the views (a worker's result files may not reach them)
    for k in range(q.shape[1]):
        ctx.telemetry("course", {"time_h": round(float(t[k * 5]), 2), "p05": round(float(q[0, k]), 3), "p25": round(float(q[1, k]), 3),
                                  "median": round(float(q[2, k]), 3), "p75": round(float(q[3, k]), 3), "p95": round(float(q[4, k]), 3),
                                  "target_low": s["lo"], "target_high": s["hi"], "toxic": toxic[0][0]})
    cnt, edges = np.histogram(o["cavg"], bins=40)
    for c_, a_, b_ in zip(cnt, edges[:-1], edges[1:]):
        ctx.telemetry("spread", {"average_mg_l": round(float(a_ + b_) / 2, 2), "patients": int(c_)})
    for r in rows:
        ctx.telemetry("groups", r)
    for i in range(min(s["n"], 300)):        # (the first 300 patients: the file has them all)
        ctx.telemetry("patients", {"patient": i + 1, "sex": "F" if pop["female"][i] else "M", "age": round(float(pop["age"][i])),
                                   "weight_kg": round(float(pop["weight"][i]), 1), "smoker": "yes" if pop["smoker"][i] else "no",
                                   "condition": condition(pop, i), "clearance_l_h": round(float(cl[i]), 2), "average_mg_l": round(float(o["cavg"][i]), 2),
                                   "peak_mg_l": round(float(o["cmax"][i]), 2), "highest_mg_l": round(float(o["peak"][i]), 2)})
    ctx.telemetry("trial", summary)
    ctx.log("At steady state %.0f%% of the patients are in the target range (%s), %.0f%% below it and %.0f%% above; %.1f%% reach a "
            "peak above %g mg/L (%s)." % (a["in_range_pct"], s["target"], a["too_low_pct"], a["too_high_pct"], a["toxic_peak_pct"], toxic[0][0], toxic[0][1]))
    for r in rows[1:]:
        ctx.log("  %s (%d): in range %.0f%%, below %.0f%%, above %.0f%%, toxic peaks %.1f%%." % (r["group"], r["patients"], r["in_range_pct"],
                                                                                         r["too_low_pct"], r["too_high_pct"], r["toxic_peak_pct"]))
    ctx.progress(100, "%.0f%% in range" % a["in_range_pct"])
    return summary


def dose_run(ctx):
    """the dose that puts the most patients of each population in range, its peaks rarely toxic"""
    import csv
    import torch
    s = settings(ctx)
    dev, label = device(ctx)
    a, b, st = s["drug"].get("dose_steps", (25.0, 2000.0, 25.0))
    steps = np.arange(1, 81) * 0.5 if s["per_kg"] else np.arange(a, b + st / 2, st)     # (the drug's doses, or 0.5-40 mg/kg)
    max_toxic = float(ctx.params.get("max_toxic_pct") or 5.0) / 100.0
    toxic = s["drug"]["toxic"][0][0]
    pops = list(POPULATIONS) if str(ctx.params.get("population_scope") or "All populations").startswith("All") else [s["population"]]
    slug = {"Healthy adults": "healthy", "Smokers": "smokers", "Elderly (65+)": "elderly", "Liver cirrhosis": "cirrhosis",
            "Heart failure": "heart_failure", "Kidney impairment": "kidney", "Mixed clinical population": "mixed"}
    rows, best = [{"dose": float(d)} for d in steps], []
    for k, name in enumerate(pops):
        ctx.progress(100 * k / len(pops), name)
        rng = np.random.default_rng(s["seed"] + k)
        pop = population(s["n"], name, rng, s["female"])
        s2 = dict(s, loading=0.0)
        m, unit, _ = unit_course(s2, pop, torch, dev)
        last = unit[0][:, -int(round(24.0 / DT)):]
        cavg1, peak1 = last.mean(1).cpu().numpy(), unit[0].max(1).values.cpu().numpy()
        pick = None
        for j, d in enumerate(steps):
            cavg, peak = cavg1 * d, peak1 * d
            r = {"population": name, "dose": float(d), "in_range_pct": round(100 * float(((cavg >= s["lo"]) & (cavg <= s["hi"])).mean()), 1),
                 "too_low_pct": round(100 * float((cavg < s["lo"]).mean()), 1), "too_high_pct": round(100 * float((cavg > s["hi"]).mean()), 1),
                 "toxic_peak_pct": round(100 * float((peak > toxic).mean()), 1), "median_average_mg_l": round(float(np.median(cavg)), 2)}
            rows[j][slug[name] + "_in_range"] = r["in_range_pct"]
            rows[j][slug[name] + "_toxic"] = r["toxic_peak_pct"]
            if r["toxic_peak_pct"] <= 100 * max_toxic and (pick is None or r["in_range_pct"] > pick["in_range_pct"]):
                pick = r
        if pick:
            best.append(dict(pick, daily=round(pick["dose"] * 24.0 / s["interval"], 1)))
            ctx.log("%s: %g %s every %g h (%g a day) puts %.0f%% in range, %.0f%% below and %.0f%% above; %.1f%% toxic peaks." % (
                name, pick["dose"], "mg/kg" if s["per_kg"] else "mg", s["interval"], pick["dose"] * 24.0 / s["interval"], pick["in_range_pct"],
                pick["too_low_pct"], pick["too_high_pct"], pick["toxic_peak_pct"]))
    for port, nm, data in (("doses", "dose_finding.csv", rows), ("best", "best_doses.csv", best)):
        with open(ctx.output(port, nm), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0]))
            w.writeheader()
            w.writerows(data)
    for r in rows:
        ctx.telemetry("doses", r)
    for b in best:
        ctx.telemetry("best", b)
    ctx.progress(100, "done")
    return {b["population"]: b["dose"] for b in best}


# ── colours ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
def hexc(rgb, k=1.0):
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(255 * x * k)))) for x in rgb)


# ── the bedside: one patient, any of the drugs, live; chosen and dosed from the 3D view's controls ───────────────────
# A profile: sex, age, height (cm), weight (kg), smoker, heart failure, cirrhosis, GFR (mL/min), and what the person
# brings to the bedside before any drug: pain (0-10), temperature (°C), airway narrowing (0-1, asthma), and a
# resting heart rate if it differs from the typical one. (Illustrative people, not cases.)
PROFILES = {
    "Healthy man 35": dict(sex="m", age=35, h=178, w=75, gfr=115),
    "Healthy woman 30": dict(sex="f", age=30, h=165, w=60, gfr=115),
    "Athlete 24 (resting HR 50)": dict(sex="m", age=24, h=182, w=78, gfr=125, hr=50),
    "Smoker 45 (20 a day)": dict(sex="m", age=45, h=176, w=80, gfr=105, smoker=1),
    "Obese woman 40 (BMI 38)": dict(sex="f", age=40, h=165, w=103, gfr=110),
    "Asthma attack 28": dict(sex="f", age=28, h=166, w=58, gfr=115, airway=0.55),
    "After surgery 50 (pain 7/10)": dict(sex="m", age=50, h=177, w=84, gfr=100, pain=7.0),
    "Pneumonia with fever 64": dict(sex="m", age=64, h=174, w=72, gfr=75, temp=39.3, pain=3.0, airway=0.2),
    "Elderly woman 84": dict(sex="f", age=84, h=158, w=55, gfr=40),
    "Frail man 91": dict(sex="m", age=91, h=168, w=52, gfr=28),
    "Liver cirrhosis 58": dict(sex="m", age=58, h=175, w=70, gfr=90, cirrhosis=1),
    "Heart failure 72": dict(sex="m", age=72, h=172, w=82, gfr=55, chf=1),
    "Kidney failure 66 (GFR 20)": dict(sex="f", age=66, h=162, w=68, gfr=20),
    "Random patient": None,
}
REFERENCE = "Healthy man 35"
TISSUE_LABEL = {"brain": "Brain", "heart": "Heart", "lung": "Lungs", "liver": "Liver", "gut": "Gut", "kidney": "Kidneys",
                "muscle": "Muscle", "adipose": "Fat", "bone": "Bone", "spleen": "Spleen", "pancreas": "Pancreas", "skin": "Skin"}
# the organs' colours: how much of the chosen drug each holds, on a log scale (tissue_scale): indigo, hardly any;
# through violet, magenta and orange to pale yellow, the most
SCALE = [(0.0, (0.27, 0.27, 0.66)), (0.25, (0.38, 0.17, 0.66)), (0.5, (0.78, 0.22, 0.55)), (0.75, (0.98, 0.50, 0.25)),
         (1.0, (0.99, 0.95, 0.65))]


def scale_rgb(v):
    v = min(1.0, max(0.0, v))
    for (a, ca), (b, cb) in zip(SCALE[:-1], SCALE[1:]):
        if v <= b:
            k = (v - a) / (b - a)
            return tuple(x + (y - x) * k for x, y in zip(ca, cb))
    return SCALE[-1][1]


def profile(name, rng):
    """(the person as the model's arrays, the extras: pain, temperature, airway, resting heart rate)"""
    p = PROFILES.get(name)
    if p is None:                         # (someone from the mixed clinical population)
        one = {k: np.asarray(v) for k, v in population(1, "Mixed clinical population", rng, 0.5).items()}
        return one, {"pain": float(rng.choice([0, 0, 0, 3, 6])), "temp": 36.8, "airway": float(rng.choice([0, 0, 0, 0.3]))}
    one = {"female": np.array([p["sex"] == "f"]), "age": np.array([float(p["age"])]), "height": np.array([float(p["h"])]),
           "weight": np.array([float(p["w"])]), "bmi": np.array([p["w"] / (p["h"] / 100.0) ** 2]),
           "smoker": np.array([bool(p.get("smoker"))]), "chf": np.array([bool(p.get("chf"))]),
           "cirrhosis": np.array([bool(p.get("cirrhosis"))]), "gfr": np.array([float(p["gfr"])]), "enzyme": np.array([1.0])}
    return one, {"pain": p.get("pain", 0.0), "temp": p.get("temp", 36.8), "airway": p.get("airway", 0.0), "hr": p.get("hr")}


def baseline_vitals(one, extra):
    """typical resting values for this person (illustrative): age, smoking, heart failure, fever and asthma move them"""
    age, chf, smk = float(one["age"][0]), bool(one["chf"][0]), bool(one["smoker"][0])
    fever = max(0.0, extra["temp"] - 37.0)
    hr = extra.get("hr") or (68.0 + (12.0 if chf else 0.0) + (6.0 if smk else 0.0))
    return {"hr": hr + 10.0 * fever, "sys": 115.0 + 0.5 * max(0.0, age - 30.0) - (10.0 if chf else 0.0),
            "dia": 74.0 + 0.12 * max(0.0, age - 30.0), "rr": 14.0 + (4.0 if chf else 0.0) + 3.0 * fever,
            "spo2": 98.0 - (2.0 if chf else 0.0) - (1.0 if smk else 0.0) - (1.0 if age > 75 else 0.0), "temp": extra["temp"],
            "qtc": 410.0 + (10.0 if one["female"][0] else 0.0) + 0.3 * max(0.0, age - 40.0)}


def occ(c, ec50, n=1.0):
    """the share of the receptors a drug occupies (Emax model with a Hill exponent)"""
    c = max(0.0, c)
    return c ** n / (c ** n + ec50 ** n) if c > 0 else 0.0


def vitals(C, b, extra):
    """the vitals from the drugs' effect-site levels C {drug: level in its unit} and the baseline b. Each drug acts
    through a few mechanisms, so drugs add up or cancel: naloxone and flumazenil compete with morphine and midazolam,
    propranolol blocks the heart-speeding (and airway-opening) effects of theophylline and salbutamol, atropine lifts
    digoxin's heart block, ondansetron stops nausea but lengthens the QT interval. The half-effect levels and
    thresholds follow the literature where it gives them (see DRUGS); the shapes of the curves are simplified."""
    g = lambda k: C.get(k, 0.0)                                           # noqa: E731
    beta_block = occ(g("Propranolol"), DRUGS["Propranolol"]["ec50"])
    opioid = g("Morphine") / (g("Morphine") + DRUGS["Morphine"]["ec50"] * (1 + g("Naloxone") / DRUGS["Naloxone"]["ki"])) if g("Morphine") else 0.0
    benzo = g("Midazolam") / (1 + g("Flumazenil") / DRUGS["Flumazenil"]["ki"])                   # (midazolam-equivalent ng/mL)
    atro = occ(g("Atropine"), DRUGS["Atropine"]["ec50"])
    salb = occ(g("Salbutamol"), DRUGS["Salbutamol"]["ec50"])
    theo = g("Theophylline")
    dig = g("Digoxin")
    amlo = occ(g("Amlodipine"), DRUGS["Amlodipine"]["ec50"])
    para = occ(g("Paracetamol"), DRUGS["Paracetamol"]["ec50"])
    onda = occ(g("Ondansetron"), DRUGS["Ondansetron"]["ec50"])
    notes = []
    # heart rate: the sympathetic drive (theophylline, salbutamol, the reflex to amlodipine) is what beta-blockade cuts
    sym = 0.3 * min(theo, 20.0) + 1.6 * max(0.0, theo - 20.0) + 25.0 * salb + 8.0 * amlo
    hr = (b["hr"] + sym * (1 - 0.8 * beta_block)) * (1 - 0.336 * beta_block)
    hr += 35.0 * atro * (1 - 0.5 * min(1.0, dig / 2.0))
    hr -= (5.0 * min(dig, 2.0) + 14.0 * max(0.0, dig - 2.0)) * (1 - 0.7 * atro)
    hr -= 6.0 * opioid
    hr = max(30.0, min(190.0, hr))
    # blood pressure
    e_mid = benzo / (benzo + 55.0) if benzo else 0.0
    sys = b["sys"] * (1 - 0.20 * e_mid) * (1 - 0.20 * amlo) * (1 - 0.12 * beta_block) - 10.0 * opioid - 1.2 * max(0.0, theo - 40.0)
    dia = b["dia"] * (1 - 0.15 * e_mid) * (1 - 0.15 * amlo) * (1 - 0.08 * beta_block) - 6.0 * opioid - 0.8 * max(0.0, theo - 40.0)
    if hr < 45:
        sys -= (45 - hr) * 1.2
    # airways: asthma narrows them; theophylline and salbutamol open them; propranolol narrows them in asthma
    airway = extra["airway"]
    if airway > 0:
        opened = 1 - (1 - 0.75 * salb * (1 - 0.9 * beta_block)) * (1 - 0.4 * occ(theo, 10.0))
        airway = min(1.0, airway * (1 - opened) + 0.6 * beta_block * extra["airway"] / max(extra["airway"], 0.3))
    # breathing: opioids and benzodiazepines depress it, together more than apart; theophylline and fever drive it
    depress = 1 - (1 - 0.45 * benzo ** 2 / (benzo ** 2 + 260.0 ** 2)) * (1 - 0.85 * opioid ** 1.3)
    rr = (b["rr"] + 0.25 * max(0.0, theo - 20.0) + 10.0 * airway) * (1 - depress)
    spo2 = b["spo2"] - 12.0 * airway - 30.0 * max(0.0, depress - 0.15)
    # sedation (Ramsay): midazolam's levels for scores 2-6, opioid sedation added on top
    sed = 1 + sum(benzo >= x for x in DRUGS["Midazolam"]["ramsay_c50"])
    sed = min(6, sed + (2 if opioid > 0.75 else 1 if opioid > 0.45 else 0))
    # pain, nausea, temperature, QT
    pain = extra["pain"] * (1 - 0.85 * opioid) * (1 - 0.35 * para)
    if g("Naloxone") > 20 and g("Morphine") > 5:
        notes.append("abrupt opioid reversal")
    nausea = min(1.0, max(0.0, (theo - 18.0) / 20.0) + max(0.0, (dig - 1.8) / 1.5) + 0.6 * opioid) * (1 - 0.9 * onda)
    temp = b["temp"] - min(1.5, max(0.0, b["temp"] - 37.0)) * para
    qtc = b["qtc"] + 25.0 * occ(g("Ondansetron"), 120.0) * 2 + 15.0 * occ(dig, 2.0)
    # rhythm
    rhythm, ectopy, dropped = "Sinus rhythm", 0.0, 0.0
    if theo > 30:
        rhythm, ectopy = "Atrial ectopic beats", min(0.4, 0.08 + 0.02 * (theo - 30.0)) * (1 - 0.7 * beta_block)
    if dig > 2.0:
        dropped = min(0.35, 0.10 * (dig - 2.0) + 0.05) * (1 - 0.8 * atro)
        if dropped > 0.03:
            rhythm = "2nd-degree AV block"
    if dig > 3.0:
        rhythm, ectopy = "AV block, ventricular ectopics", max(ectopy, min(0.4, 0.1 * (dig - 3.0) + 0.08))
    if rhythm == "Sinus rhythm":
        rhythm = "Sinus tachycardia" if hr > 100 else ("Sinus bradycardia" if hr < 55 else rhythm)
    if theo > 40:
        notes.append("seizure risk")
    if g("Atropine") > 20:
        notes.append("delirium, hot dry skin")
    if salb > 0.6 or theo > 20:
        notes.append("tremor")
    alarms = []
    for ok, what in ((hr >= 50, "HR LOW"), (hr <= 120, "HR HIGH"), (sys >= 90, "BP LOW"), (spo2 >= 92, "SpO2 LOW"),
                     (rr >= 8, "RESP LOW"), (rr <= 26, "RESP HIGH"), (ectopy + dropped < 0.03, "ARRHYTHMIA"), (qtc < 500, "QTc LONG")):
        if not ok:
            alarms.append(what)
    nz = "none" if nausea < 0.15 else "mild" if nausea < 0.45 else "moderate" if nausea < 0.75 else "severe"
    return {"hr": round(hr), "sys": round(sys), "dia": round(dia), "map": round((sys + 2 * dia) / 3), "rr": round(max(3.0, rr)),
            "spo2": round(max(65.0, min(100.0, spo2))), "temp": round(temp, 1), "qtc": round(qtc), "rhythm": rhythm,
            "ectopy": round(ectopy, 3), "dropped": round(dropped, 3), "sedation": sed,
            "sedation_txt": ["awake, anxious", "awake, calm", "drowsy, responds to voice", "asleep, responds to a tap",
                             "responds only to pain", "unrousable"][sed - 1],
            "pain": round(pain, 1), "nausea": nz, "airway": "clear" if airway < 0.1 else "wheeze" if airway < 0.45 else "severe narrowing",
            "effects": ", ".join(notes) or "none", "alarms": alarms}


def advice(v, C, b):
    """what to look at or try next, from the vitals and the drugs on board (the hint over the 3D view)"""
    def on(n):
        return n in b.drugs and C.get(n, 0.0) > 0.2 * next(iter(DRUGS[n]["range"].values()))[0]
    al = v["alarms"]
    if ("RESP LOW" in al or "SpO2 LOW" in al) and on("Morphine"):
        return "Morphine is slowing the breathing: give Naloxone (its effect is short: watch for it coming back)."
    if ("RESP LOW" in al or "SpO2 LOW" in al or v["sedation"] >= 5) and on("Midazolam"):
        return "Midazolam is sedating too deeply: give Flumazenil, or turn Regular off."
    if v["airway"] != "clear" and on("Propranolol"):
        return "Propranolol is narrowing the airways: stop it (Stop all) and give Salbutamol."
    if "SpO2 LOW" in al and v["airway"] != "clear":
        return "The airways are narrowed: give Salbutamol (or Theophylline)."
    if "ARRHYTHMIA" in al and on("Digoxin") or ("HR LOW" in al and on("Digoxin")):
        return "Digoxin is slowing the heart (heart block): give Atropine and turn Regular off."
    if "HR LOW" in al:
        return "The heart is slow: Atropine speeds it up."
    if "HR HIGH" in al or "ARRHYTHMIA" in al:
        cause = next((n for n in ("Theophylline", "Salbutamol", "Atropine", "Amlodipine") if on(n)), None)
        return ("The heart is racing%s: Propranolol slows it (not in asthma), or stop the drug." % (
            " (%s)" % cause.lower() if cause else ""))
    if "BP LOW" in al:
        return "Blood pressure is low: stop the drugs that lower it (Stop all) and wait."
    if "QTc LONG" in al:
        return "The QTc is long: avoid more Ondansetron; with Digoxin, turn Regular off."
    if v["nausea"] in ("moderate", "severe"):
        return "The patient feels sick: give Ondansetron."
    if v["pain"] >= 4:
        return "Pain %g/10: give Morphine, or Paracetamol for milder pain." % round(v["pain"])
    if v["temp"] >= 38.0:
        return "Fever: give Paracetamol."
    if v["airway"] != "clear":
        return "Wheeze: give Salbutamol or Theophylline."
    if not b.drugs or not any(on(n) for n in b.drugs):
        return "Choose a drug below and press Give dose (or Regular for the usual regimen), then watch the monitor."
    return "Stable. Raise the dose (Dose ×), add another drug, or make the clock faster and watch it build up."


# The side effects the bedside lists, as they appear and go: (key, what, where in the body, the drugs that cause it)
SIDE_EFFECTS = [
    ("breathing", "Slow, shallow breathing", "brainstem, lungs", ("Morphine", "Midazolam")),
    ("oxygen", "Low oxygen in the blood", "lungs, blood", ("Morphine", "Midazolam", "Propranolol")),
    ("sedation", "Drowsiness, sedation", "brain", ("Midazolam", "Morphine")),
    ("bronchospasm", "Narrowed airways (bronchospasm)", "lungs", ("Propranolol",)),
    ("block", "Heart block: dropped beats", "heart", ("Digoxin",)),
    ("ectopy", "Extra heartbeats (ectopics)", "heart", ("Theophylline", "Digoxin")),
    ("fast", "Fast heart rate", "heart", ("Theophylline", "Salbutamol", "Atropine", "Amlodipine")),
    ("slow", "Slow heart rate", "heart", ("Digoxin", "Propranolol", "Morphine")),
    ("qt", "Long QT interval", "heart", ("Ondansetron", "Digoxin")),
    ("bp", "Low blood pressure", "blood vessels", ("Amlodipine", "Midazolam", "Propranolol", "Morphine")),
    ("nausea", "Nausea", "gut, brain", ("Morphine", "Theophylline", "Digoxin")),
    ("vision", "Blurred, yellow-tinged vision", "eyes", ("Digoxin",)),
    ("restless", "Restlessness, insomnia", "brain", ("Theophylline",)),
    ("tremor", "Tremor", "muscles", ("Salbutamol", "Theophylline")),
    ("seizure", "Seizure risk", "brain", ("Theophylline",)),
    ("delirium", "Delirium, hot dry skin", "brain, skin", ("Atropine",)),
    ("withdrawal", "Opioid withdrawal: pain, agitation", "brain", ("Naloxone",)),
    ("liver", "Liver injury risk", "liver", ("Paracetamol",)),
    # (these come after days of use, not at once)
    ("constipation", "Constipation (after a day on it)", "gut", ("Morphine",)),
    ("oedema", "Ankle swelling (after 3 days on it)", "legs", ("Amlodipine",)),
    ("fatigue", "Tiredness, cold hands (after a day on it)", "muscles, skin", ("Propranolol",)),
]
LATE = {"constipation": ("Morphine", 0.3, 24.0), "oedema": ("Amlodipine", 0.3, 72.0), "fatigue": ("Propranolol", 0.4, 24.0)}


def side_effects_now(v, C, b, expo, dt):
    """the side effects present now: {key: the drug causing it}; expo keeps the hours each late one has been building"""
    def lvl(n):
        if n not in b.drugs:
            return 0.0
        return C.get(n, 0.0) / next(iter(DRUGS[n]["range"].values()))[1]          # (1 = the top of the target range)

    def on(n, k=0.1):
        return lvl(n) > k
    al, eff, rhythm = v["alarms"], v["effects"], v["rhythm"]
    p = b.drugs["Paracetamol"].plasma()[0] if "Paracetamol" in b.drugs else 0.0
    here = {"breathing": "RESP LOW" in al or v["rr"] < 10, "oxygen": "SpO2 LOW" in al, "sedation": v["sedation"] >= 3,
            "bronchospasm": v["airway"] != "clear" and on("Propranolol"), "block": "block" in rhythm,
            "ectopy": "ectopic" in rhythm.lower(), "fast": "HR HIGH" in al, "slow": "HR LOW" in al, "qt": v["qtc"] >= 470,
            "bp": "BP LOW" in al, "nausea": v["nausea"] != "none", "vision": C.get("Digoxin", 0.0) > 2.0,
            "restless": C.get("Theophylline", 0.0) > 20.0, "tremor": "tremor" in eff, "seizure": "seizure" in eff,
            "delirium": "delirium" in eff, "withdrawal": "opioid reversal" in eff, "liver": p > 150.0}
    for key, (drug, k, hours) in LATE.items():
        expo[key] = max(0.0, expo.get(key, 0.0) + (dt if on(drug, k) else -dt))
        here[key] = expo[key] >= hours and on(drug, k)
    out = {}
    for key, what, where, causes in SIDE_EFFECTS:
        if here.get(key):
            out[key] = next((n for n in causes if on(n)), None)
    return out


class DrugCourse:
    """one drug in the patient: its own exact linear model (the drugs do not change one another's kinetics here)"""

    def __init__(self, name, pop, torch, t, dt=DT):
        self.name, self.d, self.dt = name, DRUGS[name], dt
        self.m = patients(pop, self.d, "Oral (immediate release)")
        Phi, Gam = discretise(system_matrix(self.m), dt, torch, "cpu")
        self.Phi, self.Gam = Phi.numpy(), Gam.numpy()
        n = len(self.m["clint"])
        self.x = np.zeros((n, NS))
        self.rate = np.zeros(n)
        self.base_rate = 0.0                  # (mg/h: an infusion that runs until it is turned off)
        self.infusions = []
        self.ce = np.zeros(n)
        self.regular, self.next = False, t
        self.cum = np.zeros(3)                # (absorbed, liver, kidneys: mg)

    def give(self, t, mg, kind, over_h):
        if kind == "oral":
            self.x[:, IX["lumen"]] += mg * self.m["F"]
        else:
            h = max(over_h, DT)
            self.rate += mg / h
            self.infusions.append((t + h, mg / h))

    def step(self, t):
        for end, r in list(self.infusions):
            if t >= end - 1e-9:
                self.rate -= r
                self.infusions.remove((end, r))
        m, x = self.m, self.x
        self.cum += self.dt * np.array([m["ka"] * x[0, IX["lumen"]],
                                   m["clint"][0] * m["fu"][0] / (m["Kp"]["liver"][0] * m["V"]["liver"][0]) * x[0, IX["liver"]],
                                   m["clr"][0] / (m["Kp"]["kidney"][0] * m["V"]["kidney"][0]) * x[0, IX["kidney"]]])
        self.x = np.einsum("bij,bj->bi", self.Phi, x) + self.Gam * (self.rate + self.base_rate)[:, None]
        ke0 = self.d.get("ke0")
        cp = self.plasma()
        self.ce = cp if not ke0 else self.ce + (cp - self.ce) * (1 - math.exp(-ke0 * self.dt))

    def plasma(self):
        return self.x[:, IX["venous"]] / self.m["V"]["venous"] / self.m["bp"] * self.m["k"]

    def tissue(self, t):
        return self.x[:, IX[t]] / self.m["V"][t] * self.m["k"]

    def effect(self):
        """the level the drug acts at (patient): its effect compartment, or for digoxin the heart's own level"""
        if self.d.get("effect_at") == "heart":
            return float(self.tissue("heart")[0] / self.m["Kp"]["heart"][0])
        return float(self.ce[0])

    def flows(self):
        m, x = self.m, self.x[0]
        return (m["ka"] * x[IX["lumen"]], m["clint"][0] * m["fu"][0] / (m["Kp"]["liver"][0] * m["V"]["liver"][0]) * x[IX["liver"]],
                m["clr"][0] / (m["Kp"]["kidney"][0] * m["V"]["kidney"][0]) * x[IX["kidney"]])


class Bedside:
    """one patient (and a healthy reference adult, for comparison) with any number of drugs on board"""

    def __init__(self, patient_name, rng, torch, start=0.0):
        self.patient_name, self.rng, self.torch = patient_name, rng, torch
        one, self.extra = profile(patient_name, rng)
        ref, _ = profile(REFERENCE, rng)
        self.one = one
        self.pop = {k: np.concatenate([one[k], ref[k]]) for k in one}
        self.base = baseline_vitals(one, self.extra)
        self.t, self.carry = start, 0.0
        self.drugs = {}                       # name: DrugCourse
        self.dose_x = {}                      # name: multiple of the usual dose

    def course(self, name):
        if name not in self.drugs:
            self.drugs[name] = DrugCourse(name, self.pop, self.torch, self.t)
        return self.drugs[name]

    def give(self, name):
        mg, kind, over = DRUGS[name]["give"]
        mg *= self.dose_x.get(name, 1.0)
        self.course(name).give(self.t, mg, kind, over)
        return mg, kind

    def set_regular(self, name, on):
        c = self.course(name)
        c.regular, c.next = bool(on), self.t

    def stop_all(self):
        for c in self.drugs.values():
            c.regular = False
            c.rate[:] = 0.0
            c.infusions.clear()

    def advance(self, hours, on_dose):
        self.carry += hours
        steps = int(self.carry / DT + 1e-9)
        self.carry -= steps * DT
        for _ in range(steps):
            for name, c in self.drugs.items():
                if c.regular and self.t >= c.next - 1e-9:
                    mg, every, kind, over = DRUGS[name]["regular"]
                    mg *= self.dose_x.get(name, 1.0)
                    c.give(self.t, mg, kind, over)
                    on_dose(name, mg, kind, every, over)
                    c.next += every
                c.step(self.t)
            self.t += DT
        return steps * DT

    def why(self, name):
        """this patient's clearance and half-life of a drug, and a healthy adult's"""
        c = self.course(name)
        cl = clearance(c.m)
        vss = np.array([_vss({t: c.m["Kp"][t][i] for t in TISSUES}, {t: c.m["V"][t][i] for t in TISSUES},
                             c.m["V"]["arterial"][i] + c.m["V"]["venous"][i], c.m["bp"]) for i in range(2)])
        return cl, math.log(2) * vss / cl


def tissue_scale(c, lo, tox):
    """the organs' colour scale: plasma levels from a fifth of the range's low end to one and a half times the first
    toxic level, times the drug's typical tissue:plasma ratio (digoxin's tissues hold many times its plasma level)"""
    k = float(np.median([c.m["Kp"][t][0] for t in ORGAN_VIEW]))
    return lo / 5.0 * k, tox * 1.5 * k


def slug(name):
    return name.lower().replace(" ", "_")


def bedside_run(ctx):
    import time
    import torch
    rng = np.random.default_rng(int(time.time() * 1000) % (2 ** 31))
    focus = str(ctx.params.get("drug") or "Theophylline")
    patient = str(ctx.params.get("bedside_patient") or REFERENCE)
    patient = patient if patient in PROFILES else REFERENCE
    minutes = float(ctx.params.get("bedside_minutes") or 20.0)
    b = Bedside(patient, rng, torch)
    clocks = {"15 min/s": 0.25, "1 h/s": 1.0, "4 h/s": 4.0, "12 h/s": 12.0}       # (patient time per second)
    speed, fps = 1.0, 5.0
    state, rhythm_was, sed_was, alarms_was, pain_was, nausea_was = {}, None, None, (), None, None

    def clock(h):
        mins = int(round(h * 60.0))
        return "Day %d · %02d:%02d" % (mins // 1440 + 1, (mins // 60) % 24, mins % 60)

    def say(text):
        ctx.telemetry("events", {"text": "%s · %s" % (clock(b.t), text)})

    def how(kind, over):
        return "by mouth" if kind == "oral" else ("IV" if over < 0.2 else "IV over %g h" % over)

    def intro():
        one = b.one
        say("%s: a %s of %d, %d kg, %d cm, %s, %s, GFR %d mL/min%s." % (
            b.patient_name, "woman" if one["female"][0] else "man", round(float(one["age"][0])), round(float(one["weight"][0])),
            round(float(one["height"][0])), "smoker" if one["smoker"][0] else "non-smoker", condition(one, 0).replace("none", "no known disease"),
            round(float(one["gfr"][0])), "".join(s for s in (
                ", pain %g/10" % b.extra["pain"] if b.extra["pain"] else "", ", %.1f °C" % b.extra["temp"] if b.extra["temp"] > 37.5 else "",
                ", narrowed airways (asthma)" if b.extra["airway"] >= 0.4 else "") if s)))

    def on_dose(name, mg, kind, every, over):
        if kind == "iv" and over >= every - 1e-9:
            return                                                           # (a continuous infusion)
        say("%s %g mg %s (regular, every %g h)." % (name, mg, how(kind, over), every))

    seen, expo = {}, {}                          # (the side effects: when each started and ended; the hours the late ones have been building)
    shown = {"text": "", "until": 0.0}           # (a message across the 3D view for a few seconds: what a control just did)

    def flash(text):
        shown["text"], shown["until"] = text, time.time() + 3.5

    def handle(name, value):
        nonlocal b, focus, speed, state, rhythm_was, sed_was, alarms_was, pain_was, nausea_was
        name = str(name).strip().lower()
        if name in ("patient", "start again", "restart") and (value in PROFILES or name != "patient"):
            b = Bedside(value if value in PROFILES else b.patient_name, rng, torch)      # (the clock from Day 1 00:00)
            state, rhythm_was, sed_was, alarms_was, pain_was, nausea_was = {}, None, None, (), None, None
            ctx.telemetry("events", {"text": "─── new patient ───"})
            intro()
            flash("Started again: %s, Day 1, no drugs" % b.patient_name)
            seen.clear()
            expo.clear()
        elif name == "drug" and value in DRUGS:
            focus = value
            flash("%s chosen: press Give dose" % focus)
        elif name == "give":
            mg, kind = b.give(focus)
            over = DRUGS[focus]["give"][2]
            say("%s %g mg %s." % (focus, mg, how(kind, over)))
            flash("%s %g mg %s given" % (focus, mg, how(kind, over)))
        elif name == "regular":
            b.set_regular(focus, bool(value))
            mg, every, kind, over = DRUGS[focus]["regular"]
            if value and kind == "iv" and over >= every:
                say("%s infusion started: %g mg/h." % (focus, mg * b.dose_x.get(focus, 1.0) / every))
            elif not value:
                say("%s: regular doses stopped." % focus)
            flash("%s: regular doses %s" % (focus, "on" if value else "off"))
        elif name == "dose":
            b.dose_x[focus] = max(0.0, float(value))
            say("%s dose set to %g%% of the usual." % (focus, round(100 * b.dose_x[focus])))
        elif name in ("stop all", "stop all drugs", "stop_all"):
            b.stop_all()
            say("All drugs stopped.")
            flash("All drugs stopped")
        elif name == "clock":
            speed = clocks.get(str(value), speed)
        elif name in ("finish now", "finish", "stop"):
            return True
        return False

    intro()
    b.course(focus)
    t0, deadline, tick, done = time.time(), time.time() + minutes * 60.0, 0, False
    elapsed, given = 0.0, []                     # (hours simulated since the bedside opened, across restarts: the charts' time; the drugs given)
    beat, breath, last = 0.0, 0.0, time.time()   # (heartbeats and breaths counted in real time at the live rates: the 3D view's heart, lungs and sounds)
    while time.time() < deadline and not done:
        try:
            msgs = ctx.messages("control")
        except Exception:      # noqa: BLE001
            msgs = []
        for msg in msgs or []:
            d = msg.get("data") if isinstance(msg, dict) else None
            if isinstance(d, dict) and "control" in d:
                done = handle(d["control"], d.get("value")) or done
        dt_h = b.advance(speed / fps, on_dose)
        elapsed += dt_h
        C = {n: c.effect() for n, c in b.drugs.items()}
        v = vitals(C, b.base, b.extra)
        fc = b.course(focus)
        fd = fc.d
        lo, hi = next(iter(fd["range"].values()))
        tox = fd["toxic"][0][0]
        s_lo, s_hi = tissue_scale(fc, lo, tox)
        unit = fd["unit"]
        fmt = "%.1f" if unit == "mg/L" else ("%.2f" if hi < 10 else "%.0f")
        cp, cref = (float(x) for x in fc.plasma())
        cl, half = b.why(focus)
        h = b.t
        now_s = time.time()
        beat = (beat + (now_s - last) * v["hr"] / 60.0 * (1.0 - 0.5 * v["dropped"])) % 1200.0
        breath = (breath + (now_s - last) * v["rr"] / 60.0) % 1200.0
        last = now_s
        banner = shown["text"] if now_s < shown["until"] else ""
        item = {"hint": advice(v, C, b), "banner": banner, "banner_on": bool(banner), "t_h": round(elapsed, 2), "day": clock(h), "beat": round(beat, 3), "breath": round(breath, 3), "who": b.patient_name, "patient": b.patient_name, "drug": focus,
                "dose": b.dose_x.get(focus, 1.0), "regular": fc.regular, "clock": next((k for k, x in clocks.items() if x == speed), ""), "unit": unit, "hr": v["hr"], "rr": v["rr"], "glass": "#0a1626",
                "focus_txt": "%s: plasma %s %s (target %g-%g) · a healthy adult: %s %s" % (focus, fmt % cp, unit, lo, hi, fmt % cref, unit),
                "why": "%s: clearance %.1f L/h, half-life %.1f h (healthy adult %.1f h)" % (focus, cl[0], half[0], half[1])}
        # the drugs on board
        lines = []
        for name, c in b.drugs.items():
            dl, dh = next(iter(c.d["range"].values()))
            p = float(c.plasma()[0])
            st = "in range" if dl <= p <= dh else ("low" if p < dl else "high")
            for cut, what in c.d["toxic"]:
                if p > cut:
                    st = what
            f2 = "%.1f" if c.d["unit"] == "mg/L" else ("%.2f" if dh < 10 else "%.0f")
            if p > dl * 0.02 or c.regular or c.rate[0] > 0:
                lines.append("%s %s %s · %s%s" % (name, f2 % p, c.d["unit"], st, " · regular" if c.regular else ""))
            item["lvl_" + slug(name)] = round(100.0 * p / dh, 1)
            if name not in given:
                given.append(name)
            was = state.get(name)
            if was is not None and st != was:
                say("%s %s %s: %s." % (name, f2 % p, c.d["unit"], st))
            state[name] = st
        for name in given:                   # (a drug of a patient before Start again: 0 on the chart)
            item.setdefault("lvl_" + slug(name), 0.0)
        for k in range(6):
            item["d%d" % k] = lines[k] if k < len(lines) else ""
        if not lines:
            item["d0"] = "No drug given yet: choose one and press Give dose."
        legend = [fmt % (s_lo * (s_hi / s_lo) ** (k / 4.0)) for k in range(5)]
        for k in range(5):
            item["leg%d" % k] = legend[k]
        a_r, l_r, k_r = fc.flows()
        mg, every = (fd["regular"][0], fd["regular"][1])
        rate_ref = max(mg * b.dose_x.get(focus, 1.0), 1e-9) / every
        glow = {"gut": a_r / rate_ref, "liver": l_r / rate_ref, "kidney": k_r / rate_ref}
        for t in ORGAN_VIEW:
            ct = float(fc.tissue(t)[0])
            pos = math.log(max(ct, s_lo * 1e-3) / s_lo) / math.log(s_hi / s_lo)
            col = scale_rgb(pos)
            item["p_" + t] = round(min(1.0, max(0.0, pos)), 3)              # (the organ's bar: its place on the colour scale)
            item["c_" + t] = hexc(col)
            item["g_" + t] = hexc(col, 0.12 + 0.88 * min(1.0, glow[t]) if t in glow else 0.12)
            item["t_" + t] = (fmt % ct) + " " + unit
        pcol = scale_rgb(math.log(max(cp, lo * 1e-4) / (lo / 5.0)) / math.log(tox * 7.5 / lo))
        item["c_blood"], item["g_blood"] = hexc(pcol), hexc(pcol, 0.4)
        item["cleared"] = "%s so far: absorbed %s mg · liver %s mg · kidneys %s mg" % (
            (focus,) + tuple(("%.2f" if fd["give"][0] < 1 else "%.0f") % x for x in fc.cum))
        now_se = side_effects_now(v, C, b, expo, dt_h)
        for key, what, where, causes in SIDE_EFFECTS:
            e = seen.get(key)
            if key in now_se and (e is None or e["end"] is not None):
                seen[key] = {"start": b.t, "end": None, "cause": now_se[key]}
                say("Side effect: %s (%s)%s." % (what, where, ", from " + now_se[key].lower() if now_se[key] else ""))
            elif key not in now_se and e is not None and e["end"] is None:
                e["end"] = b.t
                say("%s: gone." % what)
        info = {k: (w, wh) for k, w, wh, _ in SIDE_EFFECTS}
        rows = sorted((e["start"], e["end"] is not None, k) for k, e in seen.items() if e["end"] is None or b.t - e["end"] < 12.0)
        for i in range(7):
            item["se%d" % i] = ""
        for i, (_, ended, k) in enumerate(rows[-7:]):
            e, (what, where) = seen[k], info[k]
            cause = " · " + e["cause"].lower() if e["cause"] else ""
            item["se%d" % i] = ("○ %s (%s)%s · %s to %s" % (what, where, cause, clock(e["start"])[-5:], clock(e["end"])[-5:]) if ended
                                else "● %s (%s)%s · since %s" % (what, where, cause, clock(e["start"])))
        item["se_none"] = "" if rows else "None so far: they appear here as they start, and stay a while after they pass."
        ctx.telemetry("patient", item)
        vit = dict(v, bp="%d/%d" % (v["sys"], v["dia"]), alarm_txt=("⚠ " + " · ".join(v["alarms"])) if v["alarms"] else "",
                   ramsay="Ramsay %d · %s" % (v["sedation"], v["sedation_txt"]), who=b.patient_name, day=item["day"], drugs=lines)
        ctx.telemetry("vitals", vit)
        # the story
        if v["rhythm"] != rhythm_was:
            if rhythm_was is not None:
                say("Rhythm: %s." % v["rhythm"])
            rhythm_was = v["rhythm"]
        if v["sedation"] != sed_was:
            if sed_was is not None:
                say("Sedation: Ramsay %d (%s)." % (v["sedation"], v["sedation_txt"]))
            sed_was = v["sedation"]
        if tuple(v["alarms"]) != alarms_was:
            new = [a for a in v["alarms"] if a not in alarms_was]
            if new:
                say("ALARM: %s." % ", ".join(new))
            alarms_was = tuple(v["alarms"])
        p_now = round(v["pain"])
        if pain_was is not None and abs(p_now - pain_was) >= 2:
            say("Pain %d/10." % p_now)
            pain_was = p_now
        elif pain_was is None:
            pain_was = p_now
        if v["nausea"] != nausea_was:
            if nausea_was not in (None, "none") and v["nausea"] != "none":      # (its start and end: the side effects)
                say("Nausea: %s." % v["nausea"])
            nausea_was = v["nausea"]
        tick += 1
        lag = t0 + tick / fps - time.time()
        if lag > 0:
            time.sleep(lag)
        ctx.progress(min(99.0, 100.0 * (time.time() - t0) / (minutes * 60.0)), "%s · %s" % (item["day"], b.patient_name))
    return {"patient": b.patient_name, "drugs": list(b.drugs), "hours": round(b.t, 1)}


# ── the anaesthetist's simulator ──────────────────────────────────────────────────────────────────────────────
# ── the anaesthetist's drugs ─────────────────────────────────────────────────────────────────────────────────────
# Given intravenously, so their effect follows the plasma with an effect compartment ("ke0", per hour: propofol and
# remifentanil equilibrate in a minute or two, rocuronium in about five). Typical adult values from their prescribing
# information and the standard pharmacokinetic literature (Schnider and Eleveld for propofol, Minto for remifentanil).
# "cl_blood_per_kg" (L/h/kg) is clearance by the blood's own esterases, not the liver: remifentanil is broken down
# wherever blood goes, which is why it stops working within minutes however long it has run. Propofol's clearance is
# about the liver's whole blood flow, so part of it is given as blood clearance too (its real clearance is partly
# outside the liver, in the lungs and kidneys).
# "bolus": (mg, over how many hours) one push from the control bar; "infusion": mg/h the slider's 1x gives;
# "per_kg": both are per kg of the patient's weight.
ANAES_DRUGS = {
    "Propofol": {     # the hypnotic: asleep in a minute, awake again in minutes however long it has run
        "logP": 3.79, "pKa": 11.0, "kind": "neutral", "fu": 0.015, "bp": 1.0, "F": 1.0,
        "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 1.10, "cl_blood_per_kg": 0.80, "cl_renal": 0.0, "vss_per_kg": 4.0,
        "unit": "mg/L", "unit_k": 1.0, "factors": {"cirrhosis": 0.85, "chf": 0.75, "smoking": 1.0},
        "range": {"Anaesthesia 2-6 mg/L": (2.0, 6.0)}, "toxic": [(8.0, "deep: hypotension, burst suppression")],
        "dose_steps": (10.0, 400.0, 10.0), "ke0": 27.4, "ec50": 3.4,
        "bolus": (2.0, 0.5 / 60), "infusion": 6.0, "per_kg": True, "mw": 178.27},
    "Remifentanil": {   # the opioid for the operation: blunts the response to the knife, gone minutes after it stops
        "logP": 1.0, "pKa": 7.07, "kind": "base", "fu": 0.30, "bp": 1.0, "F": 1.0,
        "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.0, "cl_blood_per_kg": 2.3, "cl_renal": 0.0, "vss_per_kg": 0.35,
        "unit": "ng/mL", "unit_k": 1000.0, "factors": {"cirrhosis": 1.0, "chf": 0.85, "smoking": 1.0},
        "range": {"Operating 4-12 ng/mL": (4.0, 12.0)}, "toxic": [(20.0, "chest rigidity, deep bradycardia")],
        "dose_steps": (0.05, 2.0, 0.05), "ke0": 31.0, "ec50": 12.0,
        "bolus": (0.001, 0.5 / 60), "infusion": 0.015, "per_kg": True, "mw": 376.4},
    "Fentanyl": {     # the longer opioid: minutes to work, hours to go; for pain after the operation
        "logP": 4.05, "pKa": 8.4, "kind": "base", "fu": 0.16, "bp": 1.0, "F": 1.0,
        "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.64, "cl_renal": 0.0, "vss_per_kg": 4.0,
        "unit": "ng/mL", "unit_k": 1000.0, "factors": {"cirrhosis": 0.7, "chf": 0.8, "smoking": 1.0},
        "range": {"Analgesia 1-3 ng/mL": (1.0, 3.0)}, "toxic": [(6.0, "breathing depressed after the operation")],
        "dose_steps": (0.025, 0.5, 0.025), "ke0": 8.8, "ec50": 2.0,
        "bolus": (0.0015, 0.5 / 60), "infusion": 0.0, "per_kg": True, "mw": 336.5},
    "Rocuronium": {   # the muscle relaxant: nothing moves, including the breathing
        "logP": -1.5, "pKa": 7.0, "kind": "neutral", "fu": 0.70, "bp": 1.0, "F": 1.0,
        "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.22, "cl_renal": 2.5, "vss_per_kg": 0.25,
        "unit": "mg/L", "unit_k": 1.0, "factors": {"cirrhosis": 0.7, "chf": 0.8, "smoking": 1.0},
        "range": {"Paralysis 0.3-2 mg/L": (0.3, 2.0)}, "toxic": [(4.0, "very deep block: slow to wear off")],
        "dose_steps": (5.0, 100.0, 5.0), "ke0": 14.0, "ec50": 0.35,
        "bolus": (0.6, 0.2 / 60), "infusion": 0.3, "per_kg": True, "mw": 609.7},
    "Sugammadex": {   # wraps rocuronium up and takes it out of play: a deep block gone in two minutes
        "logP": -3.0, "pKa": 3.0, "kind": "neutral", "fu": 1.0, "bp": 1.0, "F": 1.0,
        "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.0, "cl_renal": 7.0, "vss_per_kg": 0.18,
        "unit": "mg/L", "unit_k": 1.0, "factors": {"cirrhosis": 1.0, "chf": 1.0, "smoking": 1.0},
        "range": {"Reversal 20-200 mg/L": (20.0, 200.0)}, "toxic": [(600.0, "far more than is needed")],
        "dose_steps": (50.0, 1200.0, 50.0), "ke0": 60.0,
        "bolus": (2.0, 0.2 / 60), "infusion": 0.0, "per_kg": True, "mw": 2178.0},
    "Neostigmine": {   # the old reversal: works on a block that is already wearing off, not on a deep one
        "logP": -1.0, "pKa": 12.0, "kind": "neutral", "fu": 0.85, "bp": 1.0, "F": 1.0,
        "ka": {"Oral (immediate release)": 1.0},
        "cl_per_kg": 0.72, "cl_renal": 8.0, "vss_per_kg": 0.7,
        "unit": "ng/mL", "unit_k": 1000.0, "factors": {"cirrhosis": 1.0, "chf": 0.9, "smoking": 1.0},
        "range": {"Reversal 20-60 ng/mL": (20.0, 60.0)}, "toxic": [(120.0, "bradycardia, secretions")],
        "dose_steps": (0.5, 5.0, 0.5), "ke0": 7.0, "ec50": 18.0,
        "bolus": (0.04, 0.5 / 60), "infusion": 0.0, "per_kg": True, "mw": 334.4},
    "Phenylephrine": {   # tightens the vessels: the blood pressure up, the heart a little slower
        "logP": 0.3, "pKa": 9.8, "kind": "base", "fu": 0.60, "bp": 1.0, "F": 0.4,
        "ka": {"Oral (immediate release)": 1.5},
        "cl_per_kg": 1.00, "cl_renal": 0.0, "vss_per_kg": 0.34,
        "unit": "ng/mL", "unit_k": 1000.0, "factors": {"cirrhosis": 0.8, "chf": 1.0, "smoking": 1.0},
        "range": {"Pressor 2-10 ng/mL": (2.0, 10.0)}, "toxic": [(30.0, "hypertension, reflex bradycardia")],
        "dose_steps": (0.05, 0.5, 0.05), "ke0": 30.0, "ec50": 3.0,
        "bolus": (0.0013, 0.2 / 60), "infusion": 0.03, "per_kg": True, "mw": 167.2},
    "Ephedrine": {   # the heart faster and the vessels tighter: for the pressure when the heart is slow as well
        "logP": 1.13, "pKa": 9.6, "kind": "base", "fu": 0.75, "bp": 1.0, "F": 0.85,
        "ka": {"Oral (immediate release)": 1.5},
        "cl_per_kg": 0.50, "cl_renal": 10.0, "vss_per_kg": 2.5,
        "unit": "ng/mL", "unit_k": 1000.0, "factors": {"cirrhosis": 0.9, "chf": 1.0, "smoking": 1.0},
        "range": {"Pressor 30-120 ng/mL": (30.0, 120.0)}, "toxic": [(300.0, "tachycardia, hypertension")],
        "dose_steps": (2.5, 30.0, 2.5), "ke0": 20.0, "ec50": 50.0,
        "bolus": (0.08, 0.2 / 60), "infusion": 0.0, "per_kg": True, "mw": 165.2},
}
DRUGS.update(ANAES_DRUGS)
for _d in ANAES_DRUGS.values():
    _d["ka"].setdefault("Oral (sustained release)", _d["ka"]["Oral (immediate release)"])

# Four of the bedside drugs are also used in theatre and in intensive care, where they are given by the kilogram and
# run as infusions rather than prescribed by the dose. These are the doses the anaesthesia room uses for them.
AN_ALSO = {"Midazolam": {"bolus": (0.05, 0.5 / 60), "infusion": 0.06, "per_kg": True},
           "Morphine": {"bolus": (0.1, 1.0 / 60), "infusion": 0.03, "per_kg": True},
           "Ondansetron": {"bolus": (4.0, 2.0 / 60), "infusion": 0.0},
           "Atropine": {"bolus": (0.5, 0.5 / 60), "infusion": 0.0}}
for _n, _extra in AN_ALSO.items():
    DRUGS[_n].update(_extra)

# ── what the anaesthetist sees: depth, paralysis, the airway and the circulation ─────────────────────────────────
# Depth of anaesthesia as a bispectral-index-like number (BIS: 95 wide awake, 40-60 anaesthetised, under 20 burst
# suppression). The hypnotic (propofol, midazolam) sets it; an opioid on its own barely touches it but makes the
# hypnotic go much further, which is the central fact of a balanced anaesthetic: the two together do at a fraction of
# each what neither does alone. The response to the surgeon (the heart and blood pressure rising, the patient moving)
# is the other half: the opioid blunts it, and the hypnotic a little.
# The shapes here follow the published response-surface models in spirit; the numbers are rounded typical ones.
AN_C50 = {"Propofol": 3.4, "Midazolam": 400.0, "Remifentanil": 12.0}     # the level that alone gives half the effect
OPIOID_EQ = {"Remifentanil": 1.0, "Fentanyl": 1.3, "Morphine": 0.02}     # as remifentanil, ng/mL
ROC_MW, SUG_MW = 609.7, 2178.0


def an_depth(C):
    """(BIS without the surgeon's arousal, the hypnotic drive, the opioid drive as remifentanil ng/mL)"""
    uh = C.get("Propofol", 0.0) / AN_C50["Propofol"] + C.get("Midazolam", 0.0) / AN_C50["Midazolam"]
    opi = sum(k * C.get(n, 0.0) for n, k in OPIOID_EQ.items())
    uo = opi / AN_C50["Remifentanil"]
    drive = uh * (1.0 + 0.6 * uo / (1.0 + uo))            # the opioid makes the hypnotic go further
    bis = 95.0 - 93.0 * drive ** 1.6 / (1.0 + drive ** 1.6) if drive > 0 else 95.0
    return bis, uh, opi


def an_vitals(C, b, st):
    """the monitor: depth, train-of-four, oxygen, the circulation. C: the effect-site level of each drug;
    b: this patient's resting values; st: the airway, the oxygen store and how hard the surgeon is working."""
    bis0, uh, opi = an_depth(C)
    uo = (opi / AN_C50["Remifentanil"]) / (1.0 + opi / AN_C50["Remifentanil"])
    uhs = uh / (1.0 + uh)
    # the surgeon's stimulus, and how much of it the anaesthetic blocks
    a = opi / AN_C50["Remifentanil"] + 0.5 * uh
    resp = max(0.0, min(1.0, st["stimulus"] - 1.5 * a / (1.0 + a)))
    bis = min(98.0, bis0 + 14.0 * resp)
    # paralysis: the twitch left, from the rocuronium the sugammadex has not taken away, helped by neostigmine
    roc = max(0.0, st.get("roc_ce", 0.0))
    neo = occ(C.get("Neostigmine", 0.0), 18.0)
    rc = DRUGS["Rocuronium"]["ec50"]
    light = max(0.0, 1.0 - roc / (3.0 * rc))                       # neostigmine only works on a block already lifting
    roc_eff = roc / (1.0 + 3.0 * neo * light)
    twitch = 1.0 / (1.0 + (roc_eff / rc) ** 4.5)
    tof = 0 if twitch < 0.03 else 1 if twitch < 0.12 else 2 if twitch < 0.3 else 3 if twitch < 0.6 else 4
    # breathing: the hypnotic and the opioid stop it, the relaxant stops it altogether
    drive_dep = 1.6 * uh + 1.3 * (opi / AN_C50["Remifentanil"])
    spont = math.exp(-drive_dep ** 2.2) * twitch ** 2
    alv = 1.0 if st["airway"] in ("tube", "mask") else spont
    rr = 0.0 if alv <= 0.02 else (12.0 if st["airway"] == "tube" else b["rr"] * (0.5 + 0.5 * alv))
    o2 = st["o2"]
    spo2 = 100.0 - 65.0 * max(0.0, (0.5 - o2) / 0.5) ** 1.3
    # the circulation
    phe = occ(C.get("Phenylephrine", 0.0), 3.0)
    eph = occ(C.get("Ephedrine", 0.0), 50.0)
    atr = occ(C.get("Atropine", 0.0), 3.0)
    bb = occ(C.get("Propranolol", 0.0), 18.2)
    map0 = (b["sys"] + 2 * b["dia"]) / 3.0
    mean = map0 * (1.0 - 0.32 * uhs - 0.18 * uhs ** 3 - 0.10 * uo) + 38.0 * resp + 32.0 * phe + 18.0 * eph
    # (the cubic term barely shows at surgical depth and bites at burst suppression, where a pressor is the answer)
    mean -= 14.0 * max(0.0, (spo2 - 100.0) / -15.0) if spo2 < 85 else 0.0
    hr = b["hr"] * (1.0 - 0.30 * uo) * (1.0 - 0.33 * bb) + 34.0 * resp + 35.0 * atr + 32.0 * eph \
        - 16.0 * phe - 34.0 * neo * (1.0 - atr)
    hr = max(28.0, min(180.0, hr))
    mean = max(25.0, min(160.0, mean))
    sys = mean + 0.42 * mean * 0.55
    dia = mean - 0.29 * mean * 0.55
    moving = resp > 0.35 and twitch > 0.25 and bis > 55
    notes = []
    if moving:
        notes.append("moving")
    if resp > 0.25:
        notes.append("sweating, tearing")
    if C.get("Remifentanil", 0.0) > 20.0 and twitch > 0.5:
        notes.append("chest rigidity")
    if neo > 0.3 and atr < 0.2:
        notes.append("neostigmine with no atropine: the heart slows")
    return {"bis": round(bis), "tof": tof, "tof_ratio": round(twitch, 3), "hr": round(hr), "map": round(mean),
            "sys": round(sys), "dia": round(dia), "rr": round(rr), "spo2": round(spo2), "temp": b["temp"],
            "response": round(resp, 3), "moving": moving, "spont": round(spont, 3), "alv": alv,
            "uh": uh, "uo": uo, "opioid_eq": round(opi, 1), "notes": ", ".join(notes) or "none",
            "ramsay": 1 if bis > 85 else 2 if bis > 75 else 3 if bis > 65 else 4 if bis > 55 else 5 if bis > 40 else 6}

PHASE_NAMES = {"start": "Before we start", "induced": "Asleep", "surgery": "The operation", "emergence": "Waking up",
               "recovery": "Recovery", "settle": "Settling", "hold": "Holding", "wake": "The wake-up test",
               "resedate": "Settling again"}
RAMSAY = ["awake, anxious", "awake, calm", "drowsy, responds to voice", "asleep, responds to a tap",
          "responds only to pain", "unrousable"]


def mg_text(name, mg):
    return "%.0f micrograms" % (mg * 1000.0) if mg < 1.0 else "%.0f mg" % mg


def rate_text(name, mg_h, w):
    if name in ("Remifentanil", "Phenylephrine"):
        return "%.2f micrograms/kg/min" % (mg_h * 1000.0 / 60.0 / w)
    if name == "Propofol":
        return "%.1f mg/kg/h" % (mg_h / w)
    return "%.1f mg/h" % mg_h


def level_text(name, c):
    d = DRUGS[name]
    return ("%.2f %s" if d["unit"] == "mg/L" else "%.1f %s") % (c, d["unit"])


def dosing_text(name, T):
    d, bits = DRUGS[name], []
    if d.get("bolus"):
        bits.append("Give: %s" % mg_text(name, T.amount(name, "bolus") * T.dose_x.get(name, 1.0)))
    if d.get("infusion"):
        bits.append("Infusion 1x: %s" % rate_text(name, T.amount(name, "infusion"), T.weight))
    lo, hi = next(iter(d["range"].values()))
    bits.append("aiming %g-%g %s" % (lo, hi, d["unit"]))
    return " · ".join(bits)


def rhythm_of(v):
    return "Sinus tachycardia" if v["hr"] > 100 else "Sinus bradycardia" if v["hr"] < 55 else "Sinus rhythm"


def alarms_of(v, T):
    lo, hi = T.sc["band"]
    out = []
    if v["hr"] < 45:
        out.append("HR LOW")
    if v["hr"] > 110:
        out.append("HR HIGH")
    if v["map"] < max(55, 0.75 * T.base_map):
        out.append("BP LOW")
    if v["spo2"] < 92:
        out.append("SpO2 LOW")
    if v["bis"] > 65 and v["tof"] <= 1 and T.scored():
        out.append("AWAKE AND PARALYSED")
    if T.scored() and v["bis"] < lo - 10:
        out.append("TOO DEEP")
    if v["rr"] < 6 and T.airway == "none":
        out.append("NOT BREATHING")
    return out

# ── the scenarios ────────────────────────────────────────────────────────────────────────────────────────────────
# Each is a job with a target depth to hold, a few things that have to be done in the right order, and a score at the
# end. "profile": how hard the surgeon is working, as fractions of the operation's length.
AN_DT = 1.0 / 360.0                     # (h) the step the anaesthesia model takes: ten seconds
AN_SCENARIOS = {
    "Operation: a four-hour laparotomy": {
        "kind": "theatre", "band": (40, 60), "hours": 4.0, "start_airway": "none",
        "brief": "A four-hour open abdominal operation. Send the patient to sleep, take over the breathing, keep them "
                 "there while the surgeon works, then wake them up and take the tube out.",
        "profile": [(0.00, 0.03, 1.00), (0.03, 0.10, 0.75), (0.10, 0.28, 0.95), (0.28, 0.55, 0.70),
                    (0.55, 0.62, 1.00), (0.62, 0.82, 0.70), (0.82, 0.94, 0.50), (0.94, 1.01, 0.30)],
        "steps": ["Preoxygenate", "Asleep (BIS under 60)", "Paralysed (train-of-four 0-1)", "Intubated",
                  "Through the operation", "Paralysis reversed (train-of-four ratio 0.9)", "Awake and breathing",
                  "Tube out"],
    },
    "Intensive care: sedation and the wake-up test": {
        "kind": "icu", "band": (58, 80), "hours": 8.0, "start_airway": "tube",
        "brief": "A patient already on the ventilator. Settle them at light sedation (Ramsay 3-4), hold it for eight "
                 "hours, then do the daily wake-up test: stop the sedation until they are rousable, and settle them "
                 "again. Propofol lets go in minutes; midazolam, after hours of it, does not.",
        "profile": [(0.0, 1.01, 0.25)],
        "steps": ["Sedation started", "At the target depth", "Eight hours held", "Woken for the test (Ramsay 2)",
                  "Sedation restarted", "Settled again"],
    },
    "Therapeutic coma for status epilepticus": {
        "kind": "coma", "band": (12, 25), "hours": 24.0, "start_airway": "tube",
        "brief": "Fits that will not stop. Put the brain deeply to sleep (burst suppression, BIS 12-25), hold it for "
                 "twenty-four hours, then bring the patient back up. Deep enough costs blood pressure: expect to need "
                 "a pressor.",
        "profile": [(0.0, 1.01, 0.15)],
        "steps": ["Sedation started", "Burst suppression reached", "Twenty-four hours held", "Weaned off", "Awake"],
    },
}
# the right protocol for each job, written as the actions the trolley offers: what "show me the right sequence"
# puts in front of someone who has never given an anaesthetic, and what the Demo does for them if they ask
# (each line also names the step of the job it achieves, so the sequence can tick itself off with the time it
#  happened; a line with None is an aside that nothing marks done)
AN_PROTOCOL = {
    "theatre": [
        ("Preoxygenate.", "Preoxygenate"),
        ("Propofol 2 mg/kg, one push.", "Asleep (BIS under 60)"),
        ("Remifentanil 1 microgram/kg, one push.", "Asleep (BIS under 60)"),
        ("Wait until the patient is asleep (BIS under 60).", "Asleep (BIS under 60)"),
        ("Rocuronium 0.6 mg/kg, one push: this stops the breathing too.", "Paralysed (train-of-four 0-1)"),
        ("Airway: Bag-mask, and breathe for them while the relaxant works.", "Paralysed (train-of-four 0-1)"),
        ("Wait for train-of-four 0, then Airway: Intubated.", "Intubated"),
        ("Propofol infusion x1 and Remifentanil infusion x1.", "Through the operation"),
        ("Clock 1 h/s, and hold BIS 40-60 for the four hours. A push of Phenylephrine if the mean pressure "
         "falls under 65.", "Through the operation"),
        ("When the surgeon has finished: Stop every infusion, and Clock back to 1 min/s.",
         "Paralysis reversed (train-of-four ratio 0.9)"),
        ("Sugammadex 2 mg/kg, one push: the paralysis goes in two minutes.",
         "Paralysis reversed (train-of-four ratio 0.9)"),
        ("Wait for BIS over 80 and the breathing back, then Airway: No support.", "Tube out"),
    ],
    "icu": [
        ("Propofol infusion x0.65.", "Sedation started"),
        ("Clock 1 h/s, and move the infusion up or down until BIS sits at 58-80.", "At the target depth"),
        ("Hold it there for eight hours.", "Eight hours held"),
        ("Stop every infusion, and wait until they are rousable (BIS over 80): about half an hour on propofol.",
         "Woken for the test (Ramsay 2)"),
        ("Propofol infusion x0.35 again, and settle them back in the band for two hours.", "Settled again"),
        ("Try the same case on Midazolam to see why the choice of drug matters.", None),
    ],
    "coma": [
        ("Propofol infusion x1.2.", "Sedation started"),
        ("Clock 1 h/s, and raise the infusion until BIS sits at 12-25 (burst suppression).",
         "Burst suppression reached"),
        ("When the mean pressure falls under 70, Phenylephrine infusion x0.5.", None),
        ("Hold it there for twenty-four hours.", "Twenty-four hours held"),
        ("Stop every infusion and wait: after a day of propofol they take hours to come back.", "Awake"),
    ],
}
# the things that cost marks beyond the time at depth: what went wrong, and what each mistake is worth
PENALTY_FOR = {"unsafe": 10.0, "awake": 10.0, "too long to wake": 15.0, "poor intubating": 8.0}
AN_PATIENT_DRUGS = ["Propofol", "Remifentanil", "Fentanyl", "Rocuronium", "Sugammadex", "Neostigmine",
                    "Phenylephrine", "Ephedrine", "Atropine", "Midazolam", "Morphine", "Ondansetron"]
AIRWAYS = {"No support": "none", "Bag-mask": "mask", "Intubated": "tube"}


class Theatre:
    """one patient through one scenario: the drugs given, the airway, the clock and the score"""

    def __init__(self, scenario, patient_name, rng, torch):
        self.scenario = scenario if scenario in AN_SCENARIOS else list(AN_SCENARIOS)[0]
        self.sc = AN_SCENARIOS[self.scenario]
        self.patient_name = patient_name if patient_name in PROFILES else "Healthy man 35"
        one, self.extra = profile(self.patient_name, rng)
        self.one, self.pop, self.torch = one, one, torch
        self.weight = float(one["weight"][0])
        self.base = baseline_vitals(one, self.extra)
        self.base_map = (self.base["sys"] + 2.0 * self.base["dia"]) / 3.0
        self.t, self.carry = 0.0, 0.0
        self.drugs, self.inf, self.dose_x = {}, {}, {}
        self.given = set()                    # (what has actually been given: a drug merely shown has a course too)
        self.roc_ce = 0.0
        self.airway = self.sc["start_airway"]
        self.preox = self.airway == "tube"
        self.o2 = 1.0 if self.airway == "tube" else 0.52
        self.phase, self.phase_t, self.surg_t = "start", 0.0, 0.0
        self.spike_until = -1.0
        self.done = {}                        # the steps ticked off
        self.notes = []                       # what went wrong, for the score
        self.sa = dict(t=0.0, band=0.0, map_ok=0.0, hypox=0.0, severe=0.0, aware=0.0, deep=0.0, moves=0)
        self.moved_at = -1.0
        self.finished = False
        self.v = None

    # ── the drugs ──
    def course(self, name):
        if name not in self.drugs:
            self.drugs[name] = DrugCourse(name, self.pop, self.torch, self.t, dt=AN_DT)
        return self.drugs[name]

    @staticmethod
    def bolus_of(d):
        """(mg or mg/kg, how long the injection takes) — from the anaesthesia dose, or the bedside one"""
        if "bolus" in d:
            return d["bolus"][0], d["bolus"][1]
        if "give" in d:
            return d["give"][0], d["give"][2]
        return 0.0, 0.01

    def amount(self, name, what):
        d = DRUGS[name]
        mg = self.bolus_of(d)[0] if what == "bolus" else d.get("infusion", 0.0)
        return mg * (self.weight if d.get("per_kg") else 1.0)

    def give(self, name):
        mg = self.amount(name, "bolus") * self.dose_x.get(name, 1.0)
        self.course(name).give(self.t, mg, "iv", self.bolus_of(DRUGS[name])[1])
        self.given.add(name)
        return mg

    def set_infusion(self, name, mult):
        self.inf[name] = max(0.0, float(mult))
        self.course(name).base_rate = self.amount(name, "infusion") * self.inf[name]

    def stop_infusions(self):
        for name in list(self.inf):
            self.set_infusion(name, 0.0)

    # ── the surgeon ──
    def stimulus(self):
        if self.t < self.spike_until:
            return 1.2
        if self.phase == "surgery":
            f = self.surg_t / max(self.sc["hours"], 1e-6)
            for a, b, lev in self.sc["profile"]:
                if a <= f < b:
                    return lev
            return 0.3
        if self.sc["kind"] in ("icu", "coma") and self.phase in ("hold", "resedate", "sedate"):
            return self.sc["profile"][0][2]
        return 0.02

    def levels(self):
        C = {n: c.effect() for n, c in self.drugs.items()}
        C["Rocuronium"] = self.roc_ce
        return C

    def step(self, dt):
        for name, c in self.drugs.items():
            c.step(self.t)
        cp_roc = float(self.drugs["Rocuronium"].plasma()[0]) if "Rocuronium" in self.drugs else 0.0
        cp_sug = float(self.drugs["Sugammadex"].plasma()[0]) if "Sugammadex" in self.drugs else 0.0
        free = max(0.0, cp_roc / ROC_MW - cp_sug / SUG_MW) * ROC_MW      # sugammadex wraps it up, one for one
        self.roc_ce += (free - self.roc_ce) * (1.0 - math.exp(-DRUGS["Rocuronium"]["ke0"] * dt))
        self.t += dt
        self.phase_t += dt
        if self.phase == "surgery":
            self.surg_t += dt

    def oxygen(self, alv, dt):
        if alv >= 0.98 or (self.preox and alv > 0.25):
            self.o2 += (1.0 - self.o2) * (1.0 - math.exp(-40.0 * dt))
        else:
            self.o2 = max(0.0, self.o2 - (1.0 - alv) * 6.0 * dt)
            # (a fully preoxygenated adult has five or six minutes before the saturation starts to fall; one who was
            #  breathing air has about a minute and a half, which is the whole reason for preoxygenating)

    def advance(self, hours, on_event):
        """run the model forward, and move the scenario on; returns the hours actually taken"""
        self.carry += hours
        steps = int(self.carry / AN_DT + 1e-9)
        self.carry -= steps * AN_DT
        for _ in range(steps):
            self.step(AN_DT)
            st = {"airway": self.airway, "o2": self.o2, "stimulus": self.stimulus(), "roc_ce": self.roc_ce}
            v = an_vitals(self.levels(), self.base, st)
            self.v = v
            self.oxygen(v["alv"], AN_DT)
            self.watch(v, AN_DT, on_event)
        return steps * AN_DT

    # ── the score ──
    def scored(self):
        return self.phase in ("surgery", "hold", "resedate")

    def watch(self, v, dt, on_event):
        lo, hi = self.sc["band"]
        if self.scored():
            s = self.sa
            s["t"] += dt
            s["band"] += dt if lo <= v["bis"] <= hi else 0.0
            s["map_ok"] += dt if abs(v["map"] - self.base_map) <= 0.25 * self.base_map or 65 <= v["map"] <= 100 else 0.0
            s["deep"] += dt if v["map"] < 55 else 0.0
            if v["bis"] > 65 and v["tof"] <= 1:
                s["aware"] += dt
        if v["spo2"] < 90:
            self.sa["hypox"] += dt
        if v["spo2"] < 85:
            self.sa["severe"] += dt
        if v["moving"] and self.t - self.moved_at > 0.05:
            self.moved_at = self.t
            self.sa["moves"] += 1
            on_event("The patient moves.")
        self.phases(v, on_event)

    def tick(self, step, on_event, text=None):
        if step not in self.done:
            self.done[step] = self.t
            on_event(text or (step + "."))

    def phases(self, v, on_event):
        k, lo, hi = self.sc["kind"], self.sc["band"][0], self.sc["band"][1]
        if k == "theatre":
            if self.phase == "start":
                if v["bis"] < 70:
                    self.tick("Asleep (BIS under 60)", on_event, "The patient is asleep.")
                    self.phase, self.phase_t = "induced", 0.0
            elif self.phase == "induced":
                if v["tof"] <= 1:
                    self.tick("Paralysed (train-of-four 0-1)", on_event, "Fully relaxed: nothing twitches.")
                if self.airway == "tube":
                    self.phase, self.phase_t, self.surg_t = "surgery", 0.0, 0.0
                    on_event("The surgeon starts.")
            elif self.phase == "surgery":
                if self.surg_t >= self.sc["hours"]:
                    self.tick("Through the operation", on_event, "The surgeon has finished.")
                    self.phase, self.phase_t = "emergence", 0.0
            elif self.phase == "emergence":
                if v["tof_ratio"] >= 0.9:
                    self.tick("Paralysis reversed (train-of-four ratio 0.9)", on_event, "The paralysis has worn off.")
                if v["bis"] >= 80 and v["rr"] >= 8:
                    self.tick("Awake and breathing", on_event, "Awake, and breathing for themselves.")
                if self.airway == "none":
                    self.phase, self.phase_t = "recovery", 0.0
            elif self.phase == "recovery" and self.phase_t > 0.08:
                self.finish(on_event)
        else:
            if self.phase == "start":
                if v["bis"] < 90:
                    self.tick("Sedation started", on_event)
                    self.phase, self.phase_t = "settle", 0.0
            elif self.phase == "settle":
                if lo <= v["bis"] <= hi:
                    self.tick("Burst suppression reached" if k == "coma" else "At the target depth", on_event)
                    self.phase, self.phase_t = "hold", 0.0
            elif self.phase == "hold":
                if self.phase_t >= self.sc["hours"]:
                    self.tick("Twenty-four hours held" if k == "coma" else "Eight hours held", on_event)
                    self.phase, self.phase_t = "wake", 0.0
                    on_event("Time for the wake-up test: stop the sedation." if k == "icu"
                             else "Twenty-four hours done: bring the patient back up.")
            elif self.phase == "wake":
                if not any(x > 0 for x in self.inf.values()):
                    self.tick("Weaned off", on_event, "The sedation is off.")
                if v["bis"] >= 80:
                    self.tick("Woken for the test (Ramsay 2)" if k == "icu" else "Awake", on_event,
                              "The patient is rousable.")
                    if k == "coma":
                        self.finish(on_event)
                    else:
                        self.phase, self.phase_t = "resedate", 0.0
                        on_event("Now settle them again at the target.")
                elif self.phase_t > (2.0 if k == "icu" else 8.0):
                    self.notes.append("too long to wake: the test could not be done")
                    on_event("Too slow to wake: the test is over. A drug that lets go faster was the one to use.")
                    if k == "coma":
                        self.finish(on_event)
                    else:
                        self.phase, self.phase_t = "resedate", 0.0
            elif self.phase == "resedate":
                if lo <= v["bis"] <= hi and any(x > 0 for x in self.inf.values()):
                    self.tick("Sedation restarted", on_event)
                if self.phase_t >= 2.0:
                    self.tick("Settled again", on_event)
                    self.finish(on_event)

    def finish(self, on_event):
        if not self.finished:
            self.finished = True
            on_event("Finished. " + self.report()["line"])

    def report(self):
        s, steps = self.sa, self.sc["steps"]
        done = sum(1 for x in steps if x in self.done)
        band = 40.0 * (s["band"] / s["t"] if s["t"] > 0 else 0.0)
        haem = 20.0 * (s["map_ok"] / s["t"] if s["t"] > 0 else 0.0)
        oxy = max(0.0, 15.0 - 2.5 * s["hypox"] * 60.0)
        mile = 25.0 * done / len(steps)
        pen = min(30.0, 8.0 * s["aware"] * 60.0) + min(20.0, 5.0 * s["severe"] * 60.0) + 4.0 * s["moves"] \
            + min(12.0, 3.0 * s["deep"] * 60.0) + sum(PENALTY_FOR[w] for n in self.notes for w in PENALTY_FOR if w in n)
        total = max(0.0, min(100.0, band + haem + oxy + mile - pen))
        grade = "excellent" if total >= 85 else "good" if total >= 70 else "fair" if total >= 55 else "unsafe"
        return {"total": round(total), "grade": grade, "depth": round(band), "haemo": round(haem),
                "oxygen": round(oxy), "steps": round(mile), "penalty": round(pen), "done": done, "of": len(steps),
                "in_band_pct": round(100.0 * s["band"] / s["t"]) if s["t"] > 0 else 0,
                "map_pct": round(100.0 * s["map_ok"] / s["t"]) if s["t"] > 0 else 0,
                "hypox_min": round(s["hypox"] * 60.0, 1), "aware_min": round(s["aware"] * 60.0, 1),
                "moves": s["moves"], "notes": self.notes,
                "line": "Score %d out of 100 (%s): %d%% of the time at the right depth, %d%% with the blood pressure "
                        "held, %d of %d steps done." % (round(total), grade, round(100.0 * s["band"] / s["t"]) if s["t"] > 0 else 0,
                                                        round(100.0 * s["map_ok"] / s["t"]) if s["t"] > 0 else 0, done, len(steps))}

def an_hint(T, v):
    """the line across the top: what this scenario wants next"""
    if T.finished:
        return "Finished. Press Start again for another go, or change the scenario or the patient."
    k, lo, hi = T.sc["kind"], T.sc["band"][0], T.sc["band"][1]
    if v["spo2"] < 92:
        return "The oxygen is falling: ventilate (Airway: Bag-mask, or Intubated if they are relaxed enough)."
    if v["map"] < 0.8 * T.base_map:
        return "The blood pressure is low: a push of Phenylephrine (or Ephedrine if the heart is slow too), and less hypnotic."
    if k == "theatre":
        if T.phase == "start":
            if "Preoxygenate" not in T.done:
                return "Preoxygenate first: it buys you minutes once they stop breathing."
            return "Send them to sleep: Propofol, with Remifentanil for the laryngoscope."
        if T.phase == "induced":
            if v["tof"] > 1:
                return "Asleep. Now relax them for the tube: Rocuronium."
            return "Relaxed and asleep: set Airway to Intubated."
        if T.phase == "surgery":
            if v["bis"] > hi:
                return "Too light for the knife: more Propofol (and Remifentanil for the response)."
            if v["bis"] < lo:
                return "Deeper than it needs to be: turn the Propofol down."
            if v["response"] > 0.25:
                return "The patient is responding to the surgeon: more Remifentanil."
            return "Holding nicely. Keep BIS %d-%d and the blood pressure within a fifth of where it started." % (lo, hi)
        if T.phase == "emergence":
            if sum(T.inf.values()) > 0:
                return "The operation is over: turn the infusions off (Stop infusions)."
            if v["tof_ratio"] < 0.9:
                return "Still paralysed: reverse it with Sugammadex (Neostigmine only works on a block already lifting)."
            if v["bis"] < 80 or v["rr"] < 8:
                return "Reversed. Now wait for them to wake up and breathe before the tube comes out."
            return "Awake, breathing and reversed: set Airway to No support to take the tube out."
        return "In recovery."
    if T.phase == "start":
        return "Start the sedation: Propofol runs best as an infusion (Midazolam if you want to see why it is harder to wake)."
    if T.phase == "settle":
        return ("Deepen to burst suppression, BIS %d-%d." if k == "coma" else "Settle them at BIS %d-%d (Ramsay 3-4).") % (lo, hi)
    if T.phase == "hold":
        rate = T.amount("Propofol", "infusion") * T.inf.get("Propofol", 0.0) / max(T.weight, 1.0)
        if rate > 4.0 and T.phase_t > 6.0:
            return ("Over 4 mg/kg/h of propofol for hours risks propofol infusion syndrome: add Midazolam to spare it, "
                    "and keep BIS %d-%d.") % (lo, hi)
        return "Hold BIS %d-%d. Watch the blood pressure: deep sedation costs it." % (lo, hi)
    if T.phase == "wake":
        return "Stop the sedation and let them come up (Stop infusions)."
    if T.phase == "resedate":
        return "Settle them again at BIS %d-%d." % (lo, hi)
    return ""


def organ_colours(item, c, d, drug):
    """the 3D body: each organ coloured and barred by how much of the chosen drug it holds"""
    lo, hi = next(iter(d["range"].values()))
    tox = d["toxic"][0][0]
    s_lo, s_hi = tissue_scale(c, lo, tox)
    unit = d["unit"]
    fmt = "%.2f" if hi < 10 else "%.1f" if unit == "mg/L" else "%.0f"
    for k in range(5):
        item["leg%d" % k] = fmt % (s_lo * (s_hi / s_lo) ** (k / 4.0))
    item["unit"] = unit
    a_r, l_r, k_r = c.flows()
    rate_ref = max(a_r + l_r + k_r, 1e-9)
    glow = {"gut": a_r / rate_ref, "liver": l_r / rate_ref, "kidney": k_r / rate_ref}
    for t in ORGAN_VIEW:
        ct = float(c.tissue(t)[0])
        pos = math.log(max(ct, s_lo * 1e-3) / s_lo) / math.log(s_hi / s_lo)
        col = scale_rgb(pos)
        item["c_" + t] = hexc(col)
        item["g_" + t] = hexc(col, 0.12 + 0.88 * min(1.0, glow[t]) if t in glow else 0.12)
        item["p_" + t] = round(min(1.0, max(0.0, pos)), 3)
        item["t_" + t] = (fmt % ct) + " " + unit
    cp = float(c.plasma()[0])
    pcol = scale_rgb(math.log(max(cp, lo * 1e-4) / (lo / 5.0)) / math.log(tox * 7.5 / lo))
    item["c_blood"], item["g_blood"] = hexc(pcol), hexc(pcol, 0.4)


def anaesthesia_run(ctx):
    """the anaesthetist's simulator: one patient, one job, scored"""
    import time
    import torch
    rng = np.random.default_rng(int(time.time() * 1000) % (2 ** 31))
    scenario = str(ctx.params.get("an_scenario") or list(AN_SCENARIOS)[0])
    patient = str(ctx.params.get("an_patient") or "Healthy man 35")
    minutes = float(ctx.params.get("an_minutes") or 20.0)
    focus = "Propofol"
    T = Theatre(scenario, patient, rng, torch)
    clocks = {"Real time": 1.0 / 3600.0, "1 min/s": 1.0 / 60.0, "5 min/s": 5.0 / 60.0, "15 min/s": 0.25, "1 h/s": 1.0}
    speed, fps = 1.0 / 60.0, 5.0
    shown = {"text": "", "until": 0.0}
    demo = {"on": False, "next": 0.0}            # (the room playing the case itself, a move at a time)

    def clock(h):
        m = int(round(h * 60.0))                     # (round the minutes first, so 59.7 becomes 1:00, not 0:00)
        return "T + %d:%02d" % (m // 60, m % 60)

    def say(text):
        ctx.telemetry("events", {"text": "%s · %s" % (clock(T.t), text)})

    def flash(text):
        shown["text"], shown["until"] = text, time.time() + 3.5

    def restart(scen=None, who=None):
        nonlocal T
        T = Theatre(scen or T.scenario, who or T.patient_name, rng, torch)
        demo["on"] = False
        ctx.telemetry("events", {"text": "─── %s ───" % T.scenario})
        say("%s. %s" % (T.patient_name, T.sc["brief"]))
        flash("Ready: %s, %s" % (T.patient_name, T.scenario))

    def airway_to(want):
        v = T.v or an_vitals(T.levels(), T.base, {"airway": T.airway, "o2": T.o2, "stimulus": 0.0, "roc_ce": 0.0})
        if want == T.airway:
            return
        if want == "tube":
            if v["bis"] > 80:
                say("The patient is wide awake: the laryngoscope will not go in.")
                flash("Too awake to intubate")
                T.notes.append("tried to intubate an awake patient")
                T.spike_until = T.t + 1.0 / 60.0
                return
            bad = []
            if v["bis"] > 70:
                bad.append("too light (BIS %d)" % v["bis"])
            if v["tof"] > 1:
                bad.append("not relaxed (train-of-four %d)" % v["tof"])
            T.airway = "tube"
            T.preox = True
            T.spike_until = T.t + 1.0 / 60.0                 # (the laryngoscope itself is a stimulus)
            T.tick("Intubated", say, "The tube is in." if not bad else
                   "The tube is in, but the patient coughs and bucks: %s." % ", ".join(bad))
            if bad:
                T.notes.append("poor intubating conditions: " + ", ".join(bad))
                T.spike_until = T.t + 2.0 / 60.0
            flash("Intubated" if not bad else "Intubated: " + ", ".join(bad))
            return
        if T.airway == "tube" and want == "none":
            safe = v["tof_ratio"] >= 0.9 and v["bis"] >= 80 and v["rr"] >= 8
            T.airway = "none"
            T.preox = False
            if safe:
                T.tick("Tube out", say, "The tube is out, and the patient is comfortable.")
                flash("Extubated")
            else:
                why = []
                if v["tof_ratio"] < 0.9:
                    why.append("still weak")
                if v["bis"] < 80:
                    why.append("not awake")
                if v["rr"] < 8:
                    why.append("barely breathing")
                T.notes.append("unsafe extubation: " + ", ".join(why))
                T.tick("Tube out", say, "The tube is out too early: %s." % ", ".join(why))
                flash("Unsafe extubation: " + ", ".join(why))
            return
        T.airway = want
        if want == "mask":
            T.preox = True
            flash("Bag-mask ventilation")
        else:
            T.preox = False
            flash("No airway support")

    def handle(name, value):
        nonlocal focus, speed
        name = str(name).strip().lower()
        if name.startswith("give "):                 # a drug with a button of its own: "give propofol"
            want = next((d for d in AN_PATIENT_DRUGS if d.lower() == name[5:].strip()), None)
            if want:
                focus = want
                mg = T.give(want)
                say("%s %s given." % (want, mg_text(want, mg)))
                flash("%s %s" % (want, mg_text(want, mg)))
            return False
        if name == "administer" and isinstance(value, dict):      # the drug trolley: drug, how much, and how
            want = next((d for d in AN_PATIENT_DRUGS if d.lower() == str(value.get("drug", "")).lower()), None)
            if want:
                focus = want
                x = max(0.0, float(value.get("x", 1.0) or 0.0))
                if str(value.get("route", "bolus")) == "infusion":
                    T.set_infusion(want, x)
                    r = T.amount(want, "infusion") * x
                    say("%s infusion %s." % (want, rate_text(want, r, T.weight)) if r > 0 else "%s infusion off." % want)
                    flash("%s %s" % (want, rate_text(want, r, T.weight)) if r > 0 else "%s infusion off" % want)
                else:
                    T.dose_x[want] = x
                    mg = T.give(want)
                    say("%s %s given." % (want, mg_text(want, mg)))
                    flash("%s %s" % (want, mg_text(want, mg)))
            return False
        if name == "demo":
            demo["on"] = not demo["on"]
            say("Demo: the room plays the case itself." if demo["on"] else "Demo off: the controls are yours.")
            flash("Demo on · watch it once, then press Demo again to take over" if demo["on"] else "Demo off")
            return False
        if name == "scenario" and value in AN_SCENARIOS:
            restart(scen=value)
        elif name == "patient" and value in PROFILES:
            restart(who=value)
        elif name in ("start again", "restart"):
            restart()
        elif name == "drug" and value in AN_PATIENT_DRUGS:
            focus = value
            flash("%s chosen" % focus)
        elif name == "dose":
            T.dose_x[focus] = max(0.0, float(value))
            flash("%s dose set to %g x the usual" % (focus, T.dose_x[focus]))
        elif name == "give":
            mg = T.give(focus)
            say("%s %s given." % (focus, mg_text(focus, mg)))
            flash("%s %s" % (focus, mg_text(focus, mg)))
        elif name == "infusion":
            T.set_infusion(focus, value)
            r = T.amount(focus, "infusion") * T.inf[focus]
            if r > 0:
                say("%s infusion %s." % (focus, rate_text(focus, r, T.weight)))
                flash("%s %s" % (focus, rate_text(focus, r, T.weight)))
            else:
                say("%s infusion off." % focus)
                flash("%s infusion off" % focus)
        elif name in ("stop infusions", "stop all"):
            T.stop_infusions()
            say("All infusions off.")
            flash("All infusions off")
        elif name == "preoxygenate":
            T.preox = True
            T.tick("Preoxygenate", say, "Preoxygenating: the lungs fill with oxygen.")
            flash("Preoxygenating")
        elif name == "airway" and value in AIRWAYS:
            airway_to(AIRWAYS[value])
        elif name == "clock":
            speed = clocks.get(str(value), speed)
        return False

    def demo_step():
        """plays the case the way the README describes it: one move whenever the patient is ready for it,
        read off the monitor rather than off a script, so it answers to the patient in front of it"""
        k = T.sc["kind"]
        v = T.v or an_vitals(T.levels(), T.base, {"airway": T.airway, "o2": T.o2, "stimulus": 0.0, "roc_ce": 0.0})
        if k == "theatre":
            if "Preoxygenate" not in T.done:
                return handle("preoxygenate", None)
            if "Propofol" not in T.given:
                handle("give propofol", None)
                return handle("give remifentanil", None)
            if v["bis"] > 70:
                return None                                    # (waiting for them to go to sleep)
            if "Rocuronium" not in T.given:
                handle("give rocuronium", None)
                return handle("airway", "Bag-mask")
            if T.airway != "tube":
                return handle("airway", "Intubated") if v["tof"] <= 1 else None
            if T.phase == "surgery":
                if not T.inf.get("Propofol"):
                    handle("drug", "Propofol")
                    handle("infusion", 1.0)
                    handle("drug", "Remifentanil")
                    handle("infusion", 1.0)
                    handle("drug", "Propofol")
                    return handle("clock", "1 h/s")
                if v["map"] < 65:
                    return handle("give phenylephrine", None)
                x = T.inf.get("Propofol", 1.0)
                handle("drug", "Propofol")
                if v["bis"] > T.sc["band"][1]:
                    return handle("infusion", round(min(2.5, x + 0.25), 2))
                if v["bis"] < T.sc["band"][0]:
                    return handle("infusion", round(max(0.0, x - 0.25), 2))
                return None
            if T.phase == "emergence":
                if sum(T.inf.values()) > 0:
                    handle("stop infusions", None)
                    return handle("clock", "1 min/s")
                if "Sugammadex" not in T.given:
                    return handle("give sugammadex", None)
                if v["bis"] >= 80 and v["rr"] >= 8:
                    return handle("airway", "No support")
            return None
        # intensive care and the therapeutic coma: settle at the depth asked for, hold it, then let go
        lo, hi = T.sc["band"]
        if T.phase in ("start", "settle", "hold", "resedate"):
            if not T.inf.get("Propofol"):
                handle("drug", "Propofol")
                handle("infusion", 0.65 if k == "icu" else 1.2)
                return handle("clock", "1 h/s")
            if k == "coma" and v["map"] < 70 and not T.inf.get("Phenylephrine"):
                handle("drug", "Phenylephrine")
                handle("infusion", 0.5)
                return handle("drug", "Propofol")
            x = T.inf.get("Propofol", 0.0)
            handle("drug", "Propofol")
            if v["bis"] > hi:
                return handle("infusion", round(min(4.0, x + 0.15), 2))
            if v["bis"] < lo:
                return handle("infusion", round(max(0.0, x - 0.15), 2))
            return None
        if T.phase == "wake" and sum(T.inf.values()) > 0:
            return handle("stop infusions", None)
        return None

    restart()
    t0, deadline, tick = time.time(), time.time() + minutes * 60.0, 0
    beat, breath, last = 0.0, 0.0, time.time()
    while time.time() < deadline:
        try:
            msgs = ctx.messages("control")
        except Exception:      # noqa: BLE001
            msgs = []
        for msg in msgs or []:
            d = msg.get("data") if isinstance(msg, dict) else None
            if isinstance(d, dict) and "control" in d:
                handle(d["control"], d.get("value"))
        if demo["on"] and not T.finished and time.time() >= demo["next"]:
            demo["next"] = time.time() + 0.9          # (a move at a time, slow enough to follow)
            demo_step()
        if not T.finished:
            T.advance(speed / fps, say)
        v = T.v or an_vitals(T.levels(), T.base, {"airway": T.airway, "o2": T.o2, "stimulus": 0.0, "roc_ce": 0.0})
        C = T.levels()
        now_s = time.time()
        # (the view animates on the fraction of these, so wrapping them at a whole number is seamless and keeps a
        #  long case from counting into numbers the view draws badly)
        beat = (beat + (now_s - last) * v["hr"] / 60.0) % 1200.0
        breath = (breath + (now_s - last) * max(v["rr"], 0.0) / 60.0) % 1200.0
        last = now_s
        rep = T.report()
        fc = T.course(focus)
        item = {"t_h": round(T.t, 3), "day": clock(T.t), "who": T.patient_name, "scenario": T.scenario,
                "phase": PHASE_NAMES.get(T.phase, T.phase), "brief": T.sc["brief"], "hint": an_hint(T, v),
                "beat": round(beat, 3), "breath": round(breath, 3), "rr": v["rr"], "glass": "#0a1626",
                "drug": focus, "dose": T.dose_x.get(focus, 1.0), "infusion": T.inf.get(focus, 0.0),
                "airway": {"none": "No support", "mask": "Bag-mask", "tube": "Intubated"}[T.airway],
                "banner": shown["text"] if now_s < shown["until"] else "",
                "banner_on": bool(shown["text"] if now_s < shown["until"] else ""),
                "target": "Target BIS %d-%d" % T.sc["band"],
                "score": rep["total"], "grade": rep["grade"],
                "score_txt": ("Score %d / 100 · %s" % (rep["total"], rep["grade"])) if T.sa["t"] > 0
                             else "Not scored yet",
                "in_band": ("%d%% at the right depth · %d%% with the pressure held" % (rep["in_band_pct"], rep["map_pct"]))
                           if T.sa["t"] > 0 else "The marks start once you are holding the depth.",
                "focus_txt": "%s: effect site %s, plasma %s" % (
                    focus, level_text(focus, C.get(focus, 0.0)), level_text(focus, float(fc.plasma()[0]))),
                "why": "%s: %s" % (focus, dosing_text(focus, T))}
        item["inf"] = {n: round(x, 2) for n, x in T.inf.items() if x > 0}        # (what is running, for the trolley)
        if tick % 25 == 0:                                        # (every five seconds: the trolley keeps the last)
            item["weight"] = round(T.weight, 1)
            item["protocol"] = [{"do": t, "step": s} for t, s in AN_PROTOCOL.get(T.sc["kind"], [])]
            item["steps"] = list(T.sc["steps"])
        item["done_at"] = {s: clock(h) for s, h in T.done.items()}      # (when each step was reached)
            item["doses"] = {n: {"bolus": mg_text(n, T.amount(n, "bolus")),
                                 "rate": rate_text(n, T.amount(n, "infusion"), T.weight) if DRUGS[n].get("infusion") else "",
                                 "now": round(T.inf.get(n, 0.0), 2)} for n in AN_PATIENT_DRUGS}
        organ_colours(item, fc, DRUGS[focus], focus)
        nxt = next((x for x in T.sc["steps"] if x not in T.done), None)
        for i, s in enumerate(T.sc["steps"]):
            item["s%d" % i] = (("✓ " if s in T.done else "▶ " if s == nxt else "· ") + s
                               + ("  " + clock(T.done[s]) if s in T.done else ""))
        for i in range(len(T.sc["steps"]), 8):
            item["s%d" % i] = ""
        lines = []
        for name, c in T.drugs.items():
            ce = C.get(name, 0.0)
            r = T.amount(name, "infusion") * T.inf.get(name, 0.0)
            if ce > 1e-4 or r > 0:
                lines.append("%s %s%s" % (name, level_text(name, ce), " · " + rate_text(name, r, T.weight) if r > 0 else ""))
            dl, dh = next(iter(DRUGS[name]["range"].values()))
            item["lvl_" + slug(name)] = round(100.0 * ce / dh, 1)
        for k in range(6):
            item["d%d" % k] = lines[k] if k < len(lines) else ""
        if not lines:
            item["d0"] = "Nothing given yet."
        ctx.telemetry("patient", item)
        ctx.telemetry("vitals", dict(v, bp="%d/%d" % (v["sys"], v["dia"]), who=T.patient_name, day=clock(T.t),
                                     t_h=round(T.t, 3), drugs=lines, rhythm=rhythm_of(v), effects=v["notes"], nausea="none",
                                     airway={"none": "no support", "mask": "bag-mask", "tube": "ventilated"}[T.airway],
                                     sedation=v["ramsay"], sedation_txt=RAMSAY[v["ramsay"] - 1], pain=0,
                                     qtc=T.base["qtc"], alarms=alarms_of(v, T), ectopy=0.0, dropped=0.0))
        tick += 1
        lag = t0 + tick / fps - time.time()
        if lag > 0:
            time.sleep(lag)
        ctx.progress(min(99.0, 100.0 * (time.time() - t0) / (minutes * 60.0)),
                     "%s · %s · BIS %d" % (T.patient_name, PHASE_NAMES.get(T.phase, T.phase), v["bis"]))
    r = T.report()
    ctx.log("%s · %s · %s" % (T.scenario, T.patient_name, r["line"]))
    return {"scenario": T.scenario, "patient": T.patient_name, "score": r["total"], "grade": r["grade"],
            "steps_done": "%d of %d" % (r["done"], r["of"]), "at_target_pct": r["in_band_pct"],
            "pressure_held_pct": r["map_pct"], "minutes_hypoxic": r["hypox_min"], "times_moved": r["moves"]}

AN_OFFSET = (0.0061, -0.1476, -1.1497)      # (what puts the middle of patient split.glb at the scene's origin)


def an_scene_file(offset=AN_OFFSET, prefix=""):
    """the anaesthesia room's 3D view: the same patient, coloured by the drug chosen, with the anaesthetist's controls,
    the job in hand, the steps to tick off and the score so far"""
    beat, breath = "@patient.beat", "@patient.breath"
    mats = {"skin": {"opacity": 0.12, "emissive": "@patient.glass"}}
    for t in ORGAN_VIEW:
        if t != "skin":
            mats[t] = {"color": "@patient.c_%s" % t, "emissive": "@patient.g_%s" % t}
    mats["bone"]["opacity"] = 1
    mats["arterial"] = {"color": "@patient.c_blood", "emissive": "@patient.g_blood"}
    mats["venous"] = {"color": "@patient.c_blood", "emissive": "@patient.g_blood"}
    pulse = "max(0, 1 - 7 * abs((%s - floor(%s)) * 60 / max(@vitals.hr, 20) - 0.27))" % (beat, beat)
    fill = "sin(3.1416 * min(1, max(0, ((%s - floor(%s)) * 60 / max(@patient.rr, 3) - 0.3) / 2.6)))" % (breath, breath)
    parts = {"heart": {"scale": {"expr": "1 + 0.09 * %s" % pulse}},
             "lung": {"scale": [{"expr": "1 + %g * %s" % (k, fill)} for k in (0.012, 0.012, 0.03)]},
             "ribcage": {"scale": [{"expr": "1 + %g * %s" % (k, fill)} for k in (0.02, 0.03, 0.005)]}}
    organs = [{"type": "bar", "label": TISSUE_LABEL[t], "value": "@patient.p_%s" % t, "min": 0, "max": 1,
               "format": "%", "width": 180, "color": "#fa8040"}
              for t in ("brain", "heart", "lung", "liver", "kidney", "muscle", "adipose")]
    steps = [{"type": "text", "text": "@patient.s%d" % k, "size": "small"} for k in range(8)]
    return _full_width({
        "format": "esasim.scene/2", "title": "Anaesthesia", "units": "m", "up": "z",
        "renderer": {"bloom": {"strength": 0.3, "radius": 0.4, "threshold": 0.8}},
        "data": {"delay": 0.6, "smooth": 0.2, "events": "events"},
        "environment": {"sky": "space", "ground": {"type": "none"}},
        "lights": [{"type": "sun", "intensity": 1.1, "color": "#dfe8ff", "elevation": 60, "azimuth": 120},
                   {"type": "spot", "intensity": 22, "color": "#ffffff", "at": [0.0, 2.0, 1.8], "aim": [0, 0, 0], "angle": 26},
                   {"type": "spot", "intensity": 14, "color": "#5c8dff", "at": [0.0, -1.8, 1.4], "aim": [0, 0, 0], "angle": 30},
                   {"type": "ambient", "intensity": 0.3, "color": "#9fb4ff"}],
        "models": [{"id": "patient", "model": prefix + "patient split.glb", "up": "z",
                    "position": [round(v, 4) for v in offset], "materials": mats, "parts": parts}],
        "scripts": [prefix + TURN_SCRIPT],
        "sounds": [{"id": "heartbeat", "sound": prefix + HEART_SOUND, "volume": 0.8,
                    "on": {"expr": "%s - floor(%s) < 0.5" % (beat, beat)}},
                   {"id": "breath", "sound": prefix + BREATH_SOUND, "volume": 0.6,
                    "on": {"expr": "%s - floor(%s) < 0.5" % (breath, breath)}}],
        "cameras": [{"id": "front", "type": "orbit", "position": [0.55, 1.45, 0.35], "look_at": [0, 0, 0]}],
        "start": "front",
        # no control bar under the patient: the Anaesthetic trolley beside the view is where the room is worked from,
        # with the drugs, the airway and the case in groups instead of seventeen controls in a row
        "controls": [],
        "overlay": [
            {"type": "panel", "at": "top", "css": "width:250px;max-width:none", "items": [
                {"type": "text", "text": "@patient.hint", "color": "#ffe08a"}]},
            {"type": "caption", "at": "center", "text": "@patient.banner", "visible": "@patient.banner_on",
             "style": "broadcast"},
            {"type": "panel", "at": "top-left", "style": "none", "css": "width:232px;max-width:none", "items": [
                {"type": "panel", "title": "The job", "style": "clean", "css": "width:210px;max-width:none", "items": [
                    {"type": "text", "text": "@patient.who"},
                    {"type": "text", "text": "@patient.phase", "size": "large"},
                    {"type": "text", "text": "@patient.day", "color": "#9aa3b2", "size": "small"},
                    {"type": "text", "text": "@patient.target", "color": "#7dffb0", "size": "small"}] + steps + [
                    {"type": "text", "text": "@patient.score_txt", "color": "#ffe08a"},
                    {"type": "text", "text": "@patient.in_band", "color": "#9aa3b2", "size": "small"}]},
                {"type": "panel", "title": "Given", "style": "clean", "css": "width:210px;max-width:none;margin-top:10px",
                 "items": [{"type": "text", "text": "@patient.d%d" % k, "size": "small"} for k in range(6)]}]},
            # (on the right: the job and its steps already make the left side tall, and two tall panels on one side
            #  sit on top of each other when the window is short)
            {"type": "panel", "at": "top-right", "title": "Where the chosen drug is",
             "css": "width:210px;max-width:none;margin-top:36px",      # (clear of Sound, 3D mouse, Picture, Record)
             "items": [{"type": "text", "text": "@patient.focus_txt", "size": "small"},
                       {"type": "text", "text": "@patient.why", "color": "#9aa3b2", "size": "small"}] + organs},
            {"type": "events", "at": "bottom-right", "count": 4}]})


# ── the 3D view's files ────────────────────────────────────────────────────────────────────────────────────────────────────────────
HEART_SOUND = "universfield-heartbeat-single-383748.mp3"     # (the project's assets: one heartbeat, one deep breath)
BREATH_SOUND = "hasin2004-deep-breath-247459.mp3"
TURN_SCRIPT = "patient-turn.js"        # (an asset: the mouse turns the patient instead of moving the camera)


def scene_file(offset=(0.0, 0.0, 0.042), look_z=1.0, prefix="models/"):
    """the 3D view's scene: the patient, each organ coloured by how much of the chosen drug it holds, the gut, liver
    and kidneys glowing while they absorb or clear it, the heart beating (and heard) on each beat the model counts, the
    lungs filling (and heard) on each breath; the drugs on board, the events, and the controls (the vitals: the Bedside
    monitor view)"""
    beat = "@patient.beat"
    mats = {"skin": {"opacity": 0.12, "emissive": "@patient.glass"}}
    for t in ORGAN_VIEW:
        if t != "skin":
            mats[t] = {"color": "@patient.c_%s" % t, "emissive": "@patient.g_%s" % t}
    mats["bone"]["opacity"] = 1      # (one see-through layer fewer: see-through layers are drawn in an order that flips as the body turns)
    mats["arterial"] = {"color": "@patient.c_blood", "emissive": "@patient.g_blood"}
    mats["venous"] = {"color": "@patient.c_blood", "emissive": "@patient.g_blood"}
    # the heart contracts as the recording's first sound plays (0.2-0.35 s into it); the lungs fill while its breath is
    # drawn in (0.3-1.6 s into the breath recording) and empty as it is let out
    pulse = "max(0, 1 - 7 * abs((%s - floor(%s)) * 60 / max(@vitals.hr, 20) - 0.27))" % (beat, beat)
    breath = "@patient.breath"
    fill = "sin(3.1416 * min(1, max(0, ((%s - floor(%s)) * 60 / max(@patient.rr, 3) - 0.3) / 2.6)))" % (breath, breath)
    parts = {"heart": {"scale": {"expr": "1 + 0.09 * %s" % pulse}},
             "lung": {"scale": [{"expr": "1 + %g * %s" % (k, fill)} for k in (0.012, 0.012, 0.03)]}}
    # the ribcage rises with each breath (from its pivot at the spine): forwards, and a little sideways and upwards
    parts["ribcage"] = {"scale": [{"expr": "1 + %g * %s" % (k, fill)} for k in (0.02, 0.03, 0.005)]}
    organs = [{"type": "bar", "label": TISSUE_LABEL[t], "value": "@patient.p_%s" % t, "min": 0, "max": 1, "format": "%",
               "width": 210, "color": "#fa8040"} for t in ("brain", "heart", "lung", "liver", "gut", "kidney", "muscle", "adipose", "bone")]
    legend = [
        {"type": "text", "text": "■ %d%%  {@patient.leg%d} {@patient.unit}%s" % (25 * k, k, ("  hardly any", "", "", "", "  the most")[k]),
         "color": hexc(scale_rgb(k / 4.0)), "size": "small"} for k in range(5)]
    z = look_z
    return _full_width({
        "format": "esasim.scene/2", "title": "Virtual patient", "units": "m", "up": "z",
        "renderer": {"bloom": {"strength": 0.3, "radius": 0.4, "threshold": 0.8}},
        "data": {"delay": 0.6, "smooth": 0.2, "events": "events"},
        "environment": {"sky": "space", "ground": {"type": "none"}},
        "lights": [{"type": "sun", "intensity": 1.1, "color": "#dfe8ff", "elevation": 60, "azimuth": 120},
                   {"type": "spot", "intensity": 22, "color": "#ffffff", "at": [0.0, 2.0, z + 1.8], "aim": [0, 0, z], "angle": 26},
                   {"type": "spot", "intensity": 14, "color": "#5c8dff", "at": [0.0, -1.8, z + 1.4], "aim": [0, 0, z], "angle": 30},
                   {"type": "ambient", "intensity": 0.3, "color": "#9fb4ff"}],
        "models": [{"id": "patient", "model": prefix + "patient.glb", "up": "z", "position": [round(v, 4) for v in offset],
                    "materials": mats, "parts": parts}],
        "scripts": [prefix + TURN_SCRIPT],
        "sounds": [{"id": "heartbeat", "sound": prefix + HEART_SOUND, "volume": 0.8, "on": {"expr": "%s - floor(%s) < 0.5" % (beat, beat)}},
                   {"id": "breath", "sound": prefix + BREATH_SOUND, "volume": 0.6, "on": {"expr": "%s - floor(%s) < 0.5" % (breath, breath)}}],
        "cameras": [{"id": "front", "type": "orbit", "position": [0.55, 1.45, z + 0.35], "look_at": [0, 0, z]}], "start": "front",
        "controls": [
            {"type": "select", "name": "patient", "label": "Patient", "options": list(PROFILES)},
            {"type": "select", "name": "drug", "label": "Drug", "options": BEDSIDE_DRUGS},
            {"type": "slider", "name": "dose", "label": "Dose ×", "min": 0.25, "max": 4, "value": 1, "step": 0.25},
            {"type": "button", "label": "Give dose", "send": "give"},
            {"type": "toggle", "name": "regular", "label": "Regular"},
            {"type": "button", "label": "Stop all", "send": "stop all"},
            {"type": "select", "name": "clock", "label": "Clock",
             "options": ["15 min/s", "1 h/s", "4 h/s", "12 h/s"], "value": "1 h/s"},
            {"type": "button", "label": "Start again", "send": "start again"}],
        "overlay": [
            {"type": "panel", "at": "top", "css": "width:230px;max-width:none", "items": [{"type": "text", "text": "@patient.hint", "color": "#ffe08a"}]},
            {"type": "caption", "at": "center", "text": "@patient.banner", "visible": "@patient.banner_on", "style": "broadcast"},
            # (the two panels on the left stacked in one, so that they never overlap)
            {"type": "panel", "at": "top-left", "style": "none", "css": "width:232px;max-width:none", "items": [
                {"type": "panel", "title": "Virtual patient", "style": "clean", "css": "width:210px;max-width:none", "items": [
                    {"type": "text", "text": "@patient.who"},
                    {"type": "text", "text": "@patient.day", "size": "large"},
                    {"type": "text", "text": "Drugs on board", "color": "#9aa3b2", "size": "small"},
                    {"type": "text", "text": "@patient.d0", "size": "small"}, {"type": "text", "text": "@patient.d1", "size": "small"},
                    {"type": "text", "text": "@patient.d2", "size": "small"}, {"type": "text", "text": "@patient.d3", "size": "small"},
                    {"type": "text", "text": "@patient.d4", "size": "small"}, {"type": "text", "text": "@patient.d5", "size": "small"}]},
                {"type": "panel", "title": "Where the chosen drug is", "style": "clean", "css": "width:210px;max-width:none;margin-top:10px", "items": [
                    {"type": "text", "text": "@patient.focus_txt", "size": "small"},
                    {"type": "text", "text": "@patient.why", "color": "#9aa3b2", "size": "small"}] + organs + legend}]},
            {"type": "panel", "at": "top-right", "title": "Side effects", "css": "top:48px;width:250px;max-width:none;text-align:left", "items": [
                {"type": "text", "text": "@patient.se_none", "color": "#9aa3b2", "size": "small"}] + [
                {"type": "text", "text": "@patient.se%d" % k, "color": "#ffb4a8", "size": "small"} for k in range(7)]},
            {"type": "events", "at": "bottom-right", "count": 4}]})


def _full_width(scene):
    """the items inside the overlay's panels as wide as their panel (the 3D view caps each at under half of it)"""
    def free(items):
        for it in items:
            it.setdefault("css", "max-width:none")
            free(it.get("items") or [])
    for w in scene["overlay"]:
        free(w.get("items") or [])
    return scene


class _CaseFile:
    def __init__(self, path):
        import os
        self.path, self.name = path, os.path.basename(path)


def fetch_ct_case(case, log, progress=None, record="10047263"):
    """a CT case of the TotalSegmentator dataset (Zenodo record 10047263, CC BY 4.0), read in place from its zip: only
    the case's scan and labels are downloaded (for when the CT case cannot be wired in)"""
    import os
    import requests
    from remotezip import RemoteZip
    rec = requests.get("https://zenodo.org/api/records/%s" % record, timeout=60).json()
    url = next(f["links"]["self"] for f in rec["files"] if f["key"].lower().endswith(".zip"))
    import tempfile
    out = tempfile.mkdtemp(prefix="ct-case-")
    with RemoteZip(url, support_suffix_range=False) as z:
        names = [n for n in z.namelist() if ("/" + n).find("/%s/" % case) >= 0 and n.lower().endswith(".nii.gz")]
        if not names:
            raise RuntimeError("No case %s in the TotalSegmentator dataset." % case)
        for k, n in enumerate(names):
            z.extract(n, out)
            if progress:
                progress(40 * (k + 1) / len(names), "reading CT case %s (%d of %d files)" % (case, k + 1, len(names)))
    log("CT case %s read from the TotalSegmentator dataset (Zenodo %s): %d files." % (case, record, len(names)))
    return [_CaseFile(os.path.join(out, n)) for n in names]


def stream_files(ctx, items, channel="model_files", piece=40000, wait_s=180):
    """the files on a telemetry channel too (gzip, base64, in pieces), for the Model files view to put back together,
    show and offer to download: on a worker whose result files do not reach EsaSim, that is how they get out. Then,
    while the view asks for pieces it missed, they are sent again, until it has them all (or wait_s goes by)"""
    import base64
    import gzip
    import time
    packed = {}
    for name, data in items:
        z = gzip.compress(data, 6)
        packed[name] = [base64.b64encode(z[i:i + piece]).decode("ascii") for i in range(0, len(z), piece)] or [""]

    def send(name, k):
        ctx.telemetry(channel, {"file": name, "part": k, "parts": len(packed[name]), "size": len(dict(items)[name]),
                                "gz": True, "b64": packed[name][k]})
    total = sum(len(v) for v in packed.values())
    sent = 0
    for name in packed:
        for k in range(len(packed[name])):
            send(name, k)
            sent += 1
            ctx.progress(90 + 9 * sent / total, "streaming the files to the view")
            time.sleep(0.04)
    ctx.log("Streamed %d files (%d pieces) on the channel %s, for the Model files view." % (len(packed), total, channel))
    end = time.time() + wait_s
    while time.time() < end:
        try:
            msgs = ctx.messages(channel)
        except Exception:      # noqa: BLE001
            msgs = []
        for msg in msgs or []:
            d = msg.get("data") if isinstance(msg, dict) else None
            if not isinstance(d, dict):
                continue
            if d.get("complete"):
                ctx.log("The view has every file.")
                return
            for name, parts in (d.get("missing") or {}).items():
                for k in parts:
                    if name in packed and 0 <= int(k) < len(packed[name]):
                        send(name, int(k))
                        time.sleep(0.05)
        time.sleep(1.0)


def models_run(ctx):
    ct = list(ctx.inputs.get("ct", []) or [])
    if not ct and str(ctx.params.get("ct_case") or "").strip():
        ct = fetch_ct_case(str(ctx.params["ct_case"]).strip(), ctx.log, ctx.progress)
    if ct:
        ctx.log("Making the patient from the CT case wired in: its labelled structures, grouped by the model's compartments.")
        ps = ct_parts(ct, ctx.log)
    else:
        ctx.log("Making the patient (a mockup: a translucent body, its organs inside it, each a part of its own).")
        ps = parts()
    # the model's origin at the middle of the body, so that it turns about itself
    skin = next(p for p in ps if p[0] == "skin")
    mid = skin[1] + 0.5 * (skin[2].min(0) + skin[2].max(0))
    ps = [(n, c - mid, v, f, nrm) for n, c, v, f, nrm in ps]
    glb, tri = patient_glb(ps)
    made = [("patient.glb", glb), ("patient.scene.json", json.dumps(scene_file((0.0, 0.0, 0.0), 0.0, ""), indent=1).encode("utf-8"))]
    for name, data in made:
        with open(ctx.output("models", name), "wb") as fh:
            fh.write(data)
    ctx.log("patient.glb: %d parts, %d triangles, %d KB; patient.scene.json." % (
        len(tri), sum(tri.values()), len(glb) // 1024))
    stream_files(ctx, made)
    return {"parts": len(tri), "triangles": int(sum(tri.values()))}


def run(ctx):
    what = str(ctx.params.get("mode") or "Run a virtual trial").lower()
    if "anaesth" in what or "anesth" in what:
        return anaesthesia_run(ctx)
    if "model" in what:
        return models_run(ctx)
    if "dose" in what:
        return dose_run(ctx)
    if "bedside" in what or "follow" in what or "patient" in what:
        return bedside_run(ctx)
    return trial_run(ctx)
