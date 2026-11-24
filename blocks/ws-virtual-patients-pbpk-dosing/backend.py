# ── EsaSim block ──────────────────────────────────────────────────────
# Virtual patients (PBPK dosing) (ws-virtual-patients-pbpk-dosing) · version 8 · worker block, Python: def run(ctx)
# A whole-body physiologically based pharmacokinetic model of a drug in thousands of virtual patients: a trial of a dosing regimen, the best dose per population, and one patient followed live in 3D
# Workspace: Gianluca De Novi's workspace · Author: Gianluca De Novi
# Saved: 2026-11-24
#
# Inputs   ctx.inputs["id"] is a list of files: .path .name .size .read_text() .read_json()
#   ct          Files (nii)               CT case · several · optional · Make the 3D models: a CT case with its labels (ct.nii.gz and segmentations/<structure>.nii.gz, as the TotalSegmentator dataset gives them). Without it, the patient is a mockup.
# Outputs  ctx.output("id", "name.ext") is the path to write; ctx.write_json / ctx.write_text
#   telemetry   Live telemetry            ctx.telemetry("channel", data) while it runs
# Settings ctx.params["key"]
#   (none)
# Runs on  Any worker. Return a dict to fill the "result" output.
# Packages numpy, torch, scikit-image, scipy, SimpleITK, trimesh, fast-simplification (installed by the worker with pip when missing)
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


# ── drugs ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Theophylline: logP -0.02 and pKa 8.6 (a very weak acid; PKKB); 40 % bound in plasma, mostly to albumin (fu 0.6);
# oral bioavailability about 96-100 %, peak 0.5-2 h after an immediate-release tablet ("Clinical pharmacokinetics of
# theophylline", Clin Pharmacokinet 1978;3:267-93); about 85 % of
# its clearance hepatic (CYP1A2 about 88 % of that), the renal clearance 0.31 L/h in adults (Abduljalil et al. 2022);
# Vss about 0.5 L/kg and a typical adult non-smoker's clearance about 0.04 L/h/kg (half-life about 8 h).
# Concentrations: bronchodilation grows with log concentration over 5-20 mg/L; adverse effects become common above
# 20 mg/L, arrhythmias between 20 and 40 mg/L, seizures above 40 mg/L (the same 1978 review).
DRUGS = {
    "Theophylline": {
        "logP": -0.02, "pKa": 8.6, "kind": "acid", "fu": 0.60, "bp": 1.0, "F": 0.96,
        "ka": {"Oral (immediate release)": 1.5, "Oral (sustained release)": 0.25},
        "cl_per_kg": 0.040,        # (L/h/kg) a typical healthy adult non-smoker's total clearance
        "cl_renal": 0.31,          # (L/h) at a GFR of 120 mL/min
        "vss_per_kg": 0.50,        # (L/kg) the observed volume the tissue partitioning is scaled to
        "range": {"10-20 mg/L (classic)": (10.0, 20.0), "5-15 mg/L (current guidance)": (5.0, 15.0)},
        "toxic": [(20.0, "adverse effects (nausea, tachycardia)"), (30.0, "arrhythmia risk"), (40.0, "seizure risk")],
        "ec50": 10.0,              # (mg/L) the bronchodilation curve's half-way point (illustrative)
        "mw": 180.17,
    },
}


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
    act = pop["enzyme"] * np.where(pop["smoker"], SMOKING_FACTOR, 1.0) * np.where(pop["cirrhosis"], CIRRHOSIS_FACTOR, 1.0) \
        * np.where(pop["chf"], CHF_FACTOR, 1.0) * age_f
    ref_liver = np.where(pop["female"], REF["liver"][1], REF["liver"][0])
    clint = clint_ref * act * V["liver"] / REF["liver"][0] * np.where(pop["female"], REF["liver"][0] / ref_liver, 1.0)
    cl_renal = drug["cl_renal"] * pop["gfr"] / 120.0
    return {"V": V, "Q": Q, "Kp": Kp, "fu": fu, "bp": bp, "clint": clint, "clr": cl_renal, "F": drug["F"],
            "ka": drug["ka"].get(route, 1.5)}


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
        plasma[:, k] = (x[:, IX["venous"]] / vv / m["bp"]).float()
        if keep_tissues:
            tissue[:, k] = (x[:, :len(TISSUES)] / Vt).float()
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
    ("bone", ("vertebra", "rib_", "sternum", "costal", "hip_", "sacrum", "femur", "humerus", "scapula", "clavicul", "skull")),
    ("muscle", ("gluteus", "iliopsoas", "autochthon", "muscle")),
    ("rest", ("bladder", "prostate", "thyroid", "spinal_cord")),
]
PART_TRIANGLES = {"skin": 70000, "adipose": 50000, "bone": 140000, "gut": 60000, "lung": 40000, "muscle": 40000,
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
        col, rough, metal, op = PART_STYLE.get(name, ((0.8, 0.8, 0.8), 0.45, 0.0, 1.0))
        mt = {"name": name, "pbrMetallicRoughness": {"baseColorFactor": list(col) + [op], "roughnessFactor": rough, "metallicFactor": metal},
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
    lo, hi = drug["range"].get(target, drug["range"][TARGETS[0]])
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
    return q * fub * m["clint"] / (q + fub * m["clint"]) + m["clr"]


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
    steps = np.arange(1, 81) * (0.5 if s["per_kg"] else 25.0)       # (25 mg to 2 g, or 0.5 to 40 mg/kg, a dose)
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


# ── one patient, live in the 3D view ─────────────────────────────────────────────────────────────────────────────
def heat(c, lo, hi, toxic):
    """a concentration's colour: dark blue (none), blue to teal (below range), green (in range), amber (above),
    red (toxic) and hot white-pink (far above)"""
    stops = [(0.0, (0.07, 0.10, 0.22)), (0.5 * lo, (0.10, 0.35, 0.85)), (lo, (0.10, 0.75, 0.75)), (0.5 * (lo + hi), (0.25, 0.85, 0.35)),
             (hi, (0.95, 0.80, 0.15)), (toxic * 1.5 if toxic * 1.5 > hi else hi * 1.5, (0.95, 0.20, 0.12)), (toxic * 2.0, (1.0, 0.75, 0.85))]
    c = max(0.0, float(c))
    for (a, ca), (b, cb) in zip(stops[:-1], stops[1:]):
        if c <= b:
            k = (c - a) / max(b - a, 1e-9)
            return tuple(x + (y - x) * k for x, y in zip(ca, cb))
    return stops[-1][1]


def hexc(rgb, k=1.0):
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(round(255 * x * k)))) for x in rgb)


def follow_run(ctx):
    import time
    import torch
    s = settings(ctx)
    dev, label = device(ctx)
    rng = np.random.default_rng(s["seed"])
    n = max(50, min(s["n"], 4000))
    pop = population(n, s["population"], rng, s["female"])
    m, unit, first = unit_course(s, pop, torch, dev)
    pl, _ = scaled(s, unit, first)
    o = outcomes(pl, DT, s["interval"], s["lo"], s["hi"])
    who = str(ctx.params.get("follow") or "The typical patient")
    cavg = o["cavg"].cpu().numpy()
    if who.startswith("The most exposed"):
        i = int(np.argmax(o["peak"].cpu().numpy()))
    elif who.startswith("The least"):
        i = int(np.argmin(cavg))
    elif who.startswith("A random"):
        i = int(rng.integers(n))
    else:
        i = int(np.argsort(cavg)[n // 2])
    one = {k: v[i:i + 1] for k, v in pop.items()}
    m1 = patients(one, s["drug"], s["route"])
    doses, _ = regimen(s["route"], s["dose"], s["interval"], s["days"], s["loading"], s["infusion_h"])
    w = one["weight"] if s["per_kg"] else None
    plasma, tissue = simulate(m1, doses, s["days"], torch, "cpu", dt=DT, infusion_h=s["infusion_h"], per_kg=w, keep_tissues=True)
    plasma, tissue = plasma[0].numpy(), tissue[0].numpy()
    kp = np.array([float(m1["Kp"][t][0]) for t in TISSUES])
    # what each organ is exposed to: its concentration as the plasma concentration it is in balance with (mg/L)
    expo = tissue / kp[None]
    cl = float(clearance(m1)[0])
    toxic = s["drug"]["toxic"]
    sex = "woman" if one["female"][0] else "man"
    info = {"who": "%s: a %s of %d, %d kg, %d cm" % (who, sex, round(one["age"][0]), round(one["weight"][0]), round(one["height"][0])),
            "smoker": "smoker" if one["smoker"][0] else "non-smoker", "condition": condition(one, 0), "gfr": round(float(one["gfr"][0])),
            "clearance": round(cl, 2), "half_life": round(math.log(2) * s["drug"]["vss_per_kg"] * float(one["weight"][0]) / cl, 1),
            "regimen": describe(s), "target": "%g-%g mg/L" % (s["lo"], s["hi"])}
    ctx.log("Following %s (%s, %s; GFR %d mL/min; clearance %.2f L/h, half-life %.1f h) on %s." % (
        info["who"], info["smoker"], info["condition"], info["gfr"], cl, info["half_life"], info["regimen"]))
    speed = float(ctx.params.get("speed") or 1.0)          # (simulated hours per second)
    fps = 5.0
    hr = 72 + (12 if one["chf"][0] else 0) + (6 if one["smoker"][0] else 0)
    rr = 14 + (4 if one["chf"][0] else 0)
    ch = "patient"
    dose_at = {int(round(t0 / DT)): (amt * (float(one["weight"][0]) if s["per_kg"] else 1.0), kind) for t0, amt, kind in doses}
    state = "start"
    t0 = time.time()
    total = plasma.shape[0]
    k = 0
    ctx.telemetry("events", {"text": "%s, %s, %s; clearance %.2f L/h (half-life %.1f h)." % (info["who"], info["smoker"], info["condition"], cl, info["half_life"])})
    while k < total:
        h = k * DT
        for j in range(max(0, k - int(round(speed / fps / DT)) + 1), k + 1):
            if j in dose_at:
                amt, kind = dose_at[j]
                ctx.telemetry("events", {"text": "Day %d %02d:00 · %s %.0f mg %s" % (j * DT // 24 + 1, int(j * DT % 24), s["drug_name"], amt,
                                                                                  "by mouth" if kind == "oral" else "IV")})
        c = float(plasma[k])
        now = "in range" if s["lo"] <= c <= s["hi"] else ("below range" if c < s["lo"] else "above range")
        for cut, what in toxic:
            if c > cut:
                now = what
        if now != state and k > 0:
            ctx.telemetry("events", {"text": "Day %d %02d:00 · plasma %.1f mg/L: %s" % (h // 24 + 1, int(h % 24), c, now)})
        state = now
        item = {"t_h": round(h, 2), "day": "Day %d · %02d:%02d" % (h // 24 + 1, int(h % 24), int(round((h % 1) * 60)) % 60),
                "plasma": round(c, 2), "state": now, "hr": hr, "rr": rr, "glass": hexc((0.35, 0.6, 1.0), 0.25)}
        item.update(info)
        for j, t in enumerate(TISSUES):
            e = float(expo[k, j])
            col = heat(e, s["lo"], s["hi"], toxic[0][0])
            item["c_" + t] = hexc(col)
            item["g_" + t] = hexc(col, 0.35 + 0.65 * min(1.0, max(0.0, (e - s["hi"]) / max(toxic[0][0], 1e-6))))
            item[t] = round(e, 2)
        col = heat(c, s["lo"], s["hi"], toxic[0][0])
        item["c_blood"], item["g_blood"] = hexc(col), hexc(col, 0.5)
        ctx.telemetry(ch, item)
        k += max(1, int(round(speed / fps / DT)))
        lag = t0 + (k / max(1, int(round(speed / fps / DT)))) / fps - time.time()
        if lag > 0:
            time.sleep(lag)
        ctx.progress(100 * k / total, "%s · %.1f mg/L" % (item["day"], c))
        if finish_asked(ctx):
            break
    return {"patient": info["who"], "clearance_l_h": info["clearance"], "half_life_h": info["half_life"],
            "average_mg_l": round(float(o["cavg"][i]), 2), "peak_mg_l": round(float(o["peak"][i]), 2)}


def finish_asked(ctx, ch="control"):
    try:
        for msg in ctx.messages(ch):
            d = msg.get("data") if isinstance(msg, dict) else None
            if isinstance(d, dict) and str(d.get("control", "")).lower() in ("finish now", "finish", "stop"):
                return True
    except Exception:      # noqa: BLE001
        pass
    return False


# ── the 3D view's files ──────────────────────────────────────────────────────────────────────────────────────────
def platform_glb():
    """the round platform the patient stands on: dark glass, a glowing ring at its edge"""
    tri, nrm = [], []
    seg, R, r_in, h = 96, 0.55, 0.52, 0.04
    for k in range(seg):
        a, b = 2 * math.pi * k / seg, 2 * math.pi * (k + 1) / seg
        ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
        tri += [(0, 0, h), (r_in * ca, r_in * sa, h), (r_in * cb, r_in * sb, h)]
        nrm += [(0, 0, 1)] * 3
        tri += [(r_in * ca, r_in * sa, h), (R * ca, R * sa, h), (R * cb, R * sb, h), (r_in * ca, r_in * sa, h), (R * cb, R * sb, h), (r_in * cb, r_in * sb, h)]
        nrm += [(0, 0, 1)] * 6
        tri += [(R * ca, R * sa, 0), (R * cb, R * sb, 0), (R * cb, R * sb, h), (R * ca, R * sa, 0), (R * cb, R * sb, h), (R * ca, R * sa, h)]
        nrm += [(ca, sa, 0)] * 6
    return _glb_parts([("platform", np.array(tri), np.array(nrm), (0.05, 0.06, 0.08), 0.15, 0.6, 1.0, None)])


def ring_glb():
    tri, nrm = [], []
    seg, R, w, h = 128, 0.555, 0.012, 0.042
    for k in range(seg):
        a, b = 2 * math.pi * k / seg, 2 * math.pi * (k + 1) / seg
        p = [(R * math.cos(a), R * math.sin(a)), (R * math.cos(b), R * math.sin(b)), ((R - w) * math.cos(b), (R - w) * math.sin(b)), ((R - w) * math.cos(a), (R - w) * math.sin(a))]
        tri += [(p[0][0], p[0][1], h), (p[1][0], p[1][1], h), (p[2][0], p[2][1], h), (p[0][0], p[0][1], h), (p[2][0], p[2][1], h), (p[3][0], p[3][1], h)]
        nrm += [(0, 0, 1)] * 6
    return _glb_parts([("ring", np.array(tri), np.array(nrm), (0.3, 0.6, 1.0), 0.3, 0.0, 1.0, (0.3, 0.6, 1.0))])


def _glb_parts(items):
    """a glTF binary of separate parts: [(name, triangle corners (n, 3), normals, colour, roughness, metalness,
    opacity, emissive or None)]"""
    blobs, views, accs, meshes, nodes, mats = [], [], [], [], [], []
    off = 0

    def blob(b, target):
        nonlocal off
        views.append({"buffer": 0, "byteOffset": off, "byteLength": len(b), "target": target})
        blobs.append(b + b"\0" * (-len(b) % 4))
        off += len(blobs[-1])
        return len(views) - 1
    for k, (name, p, n, col, rough, metal, op, emi) in enumerate(items):
        p, n = p.astype(np.float32), n.astype(np.float32)
        idx = np.arange(len(p), dtype=np.uint32)
        a0 = len(accs)
        accs.append({"bufferView": blob(p.tobytes(), 34962), "componentType": 5126, "count": len(p), "type": "VEC3", "min": p.min(0).tolist(), "max": p.max(0).tolist()})
        accs.append({"bufferView": blob(n.tobytes(), 34962), "componentType": 5126, "count": len(n), "type": "VEC3"})
        accs.append({"bufferView": blob(idx.tobytes(), 34963), "componentType": 5125, "count": len(idx), "type": "SCALAR"})
        mt = {"name": name, "pbrMetallicRoughness": {"baseColorFactor": list(col) + [op], "roughnessFactor": rough, "metallicFactor": metal}}
        if emi:
            mt["emissiveFactor"] = list(emi)
        mats.append(mt)
        meshes.append({"name": name, "primitives": [{"attributes": {"POSITION": a0, "NORMAL": a0 + 1}, "indices": a0 + 2, "material": k}]})
        nodes.append({"name": name, "mesh": k})
    js = {"asset": {"version": "2.0"}, "scene": 0, "scenes": [{"nodes": list(range(len(nodes)))}], "nodes": nodes, "meshes": meshes,
          "materials": mats, "bufferViews": views, "accessors": accs, "buffers": [{"byteLength": off}]}
    j = json.dumps(js, separators=(",", ":")).encode()
    j += b" " * (-len(j) % 4)
    binc = b"".join(blobs)
    return (struct.pack("<III", 0x46546C67, 2, 28 + len(j) + len(binc)) + struct.pack("<I", len(j)) + b"JSON" + j
            + struct.pack("<I", len(binc)) + b"BIN\0" + binc)


def heartbeat_wav(rate=32000):
    """one heartbeat: the two soft, low thumps of a heart heard through a stethoscope (lub, then dub)"""
    t = np.arange(int(0.45 * rate)) / rate
    rng = np.random.default_rng(3)

    def thump(t0, f, a, d):
        tt, on = np.clip(t - t0, 0, None), (t >= t0).astype(float)
        return a * on * np.sin(2 * np.pi * f * tt) * np.exp(-tt / d) * (1 - np.exp(-tt / 0.004))
    x = thump(0.0, 55, 1.0, 0.045) + thump(0.0, 110, 0.3, 0.025) + thump(0.21, 70, 0.65, 0.035) + thump(0.21, 140, 0.2, 0.02)
    x = x + 0.03 * rng.standard_normal(len(t)) * np.exp(-t / 0.05)
    x = x / (np.abs(x).max() + 1e-9) * 0.9
    pcm = (x * 32767).astype("<i2").tobytes()
    return (b"RIFF" + struct.pack("<I", 36 + len(pcm)) + b"WAVEfmt " + struct.pack("<IHHIIHH", 16, 1, 1, rate, rate * 2, 2, 16)
            + b"data" + struct.pack("<I", len(pcm)) + pcm)


def scene_file(torso=False):
    """the 3D view's scene (models/patient.scene.json): the patient on the platform, each part coloured by the drug's
    concentration there (from the channel patient), the heart beating and the lungs breathing at the patient's rates"""
    beat = "max(0, sin(@t * 6.2832 * (@patient.hr or 72) / 60))"
    mats = {"skin": {"opacity": 0.13, "emissive": "@patient.glass"}}
    for t in ORGAN_VIEW:
        if t in ("skin",):
            continue
        mats[t] = {"color": "@patient.c_%s" % t, "emissive": "@patient.g_%s" % t}
    mats["arterial"] = {"color": "@patient.c_blood", "emissive": "@patient.g_blood"}
    mats["venous"] = {"color": "@patient.c_blood", "emissive": "@patient.g_blood"}
    parts = {"heart": {"scale": {"expr": "1 + 0.07 * %s * %s" % (beat, beat)}},
             "lung": {"scale": {"expr": "1 + 0.035 * sin(@t * 6.2832 * (@patient.rr or 14) / 60)"}}}
    organ_bars = [{"type": "bar", "label": nm.capitalize(), "value": "@patient.%s" % t, "unit": "mg/L", "min": 0, "max": 30, "width": 150,
                   "color": "#5c8dff"} for nm, t in (("brain", "brain"), ("heart", "heart"), ("lungs", "lung"), ("liver", "liver"), ("gut", "gut"),
                                                     ("kidneys", "kidney"), ("muscle", "muscle"), ("fat", "adipose"), ("bone", "bone"))]
    return {
        "format": "esasim.scene/2", "title": "Virtual patient", "units": "m", "up": "z",
        "renderer": {"bloom": {"strength": 0.55, "radius": 0.45, "threshold": 0.55}},
        "data": {"delay": 1.0, "smooth": 0.2, "events": "events"},
        "environment": {"sky": "space", "ground": {"type": "none"}},
        "lights": [{"type": "sun", "intensity": 1.2, "color": "#dfe8ff", "elevation": 60, "azimuth": 120},
                   {"type": "spot", "intensity": 25, "color": "#ffffff", "at": [0.0, 2.2, 3.2], "aim": [0, 0, 1.0], "angle": 28},
                   {"type": "spot", "intensity": 18, "color": "#5c8dff", "at": [0.0, -2.0, 2.6], "aim": [0, 0, 1.1], "angle": 30},
                   {"type": "ambient", "intensity": 0.25, "color": "#9fb4ff"}],
        "models": [{"id": "platform", "model": "models/platform.glb", "up": "z", "position": [0, 0, 0]},
                   {"id": "ring", "model": "models/ring.glb", "up": "z", "position": [0, 0, 0],
                    "materials": {"ring": {"emissive": "@patient.c_blood"}}},
                   {"id": "patient", "model": "models/patient.glb", "up": "z", "position": [0, 0, 0.042], "materials": mats, "parts": parts}],
        "sounds": [{"id": "heartbeat", "sound": "models/heartbeat.wav", "volume": 0.25,
                    "on": {"expr": "sin(@t * 6.2832 * (@patient.hr or 72) / 60) > 0.97"}}],
        "cameras": [{"id": "front", "type": "orbit", "position": [0.45, 1.25, 1.3] if torso else [0.9, 2.6, 1.35],
                     "look_at": [0, 0, TORSO_Z] if torso else [0, 0, 1.0]}], "start": "front",
        "overlay": [
            {"type": "panel", "at": "top-left", "title": "Virtual patient", "items": [
                {"type": "text", "text": "@patient.who"},
                {"type": "text", "text": "{@patient.smoker} · {@patient.condition} · GFR {@patient.gfr} mL/min", "color": "#9aa3b2"},
                {"type": "text", "text": "Clearance {@patient.clearance} L/h · half-life {@patient.half_life} h", "color": "#9aa3b2"},
                {"type": "text", "text": "@patient.regimen", "color": "#9aa3b2"}]},
            {"type": "panel", "at": "bottom-left", "items": [
                {"type": "text", "text": "@patient.day", "size": "large"},
                {"type": "tiles", "items": [{"label": "Plasma", "value": "@patient.plasma", "unit": "mg/L", "format": "0.0"},
                                            {"label": "Now", "value": "@patient.state"},
                                            {"label": "Target", "value": "@patient.target"}]},
                {"type": "chart", "title": "Plasma concentration (mg/L)", "value": "@patient.plasma", "seconds": 90, "width": 300, "height": 90}]},
            {"type": "panel", "at": "right", "title": "Organs (mg/L in balance)", "items": organ_bars},
            {"type": "events", "at": "bottom-right", "count": 5}]}


def models_run(ctx):
    ct = list(ctx.inputs.get("ct", []) or [])
    if ct:
        ctx.log("Making the patient from the CT case wired in: its labelled structures, grouped by the model's compartments.")
        glb, tri = patient_glb(ct_parts(ct, ctx.log))
    else:
        ctx.log("Making the patient (a mockup: a translucent body, its organs inside it, each a part of its own).")
        glb, tri = patient_glb()
    for name, data in (("patient.glb", glb), ("platform.glb", platform_glb()), ("ring.glb", ring_glb()), ("heartbeat.wav", heartbeat_wav())):
        with open(ctx.output("models", name), "wb") as fh:
            fh.write(data)
    with open(ctx.output("models", "patient.scene.json"), "w", encoding="utf-8") as fh:
        json.dump(scene_file(torso=bool(ct)), fh, indent=1)
    ctx.log("patient.glb: %d parts, %d triangles, %d KB; platform.glb, ring.glb, heartbeat.wav, patient.scene.json." % (
        len(tri), sum(tri.values()), len(glb) // 1024))
    return {"parts": len(tri), "triangles": int(sum(tri.values()))}


def run(ctx):
    what = str(ctx.params.get("mode") or "Run a virtual trial").lower()
    if "model" in what:
        return models_run(ctx)
    if "dose" in what:
        return dose_run(ctx)
    if "follow" in what or "patient" in what:
        return follow_run(ctx)
    return trial_run(ctx)
