# ── EsaSim block ──────────────────────────────────────────────────────
# Temporary: patient models into worker storage (ws-temporary-patient-models-into-worker) · version 1 · worker block, Python: def run(ctx)
# One-off: makes the virtual patient's 3D files from a CT case already on this worker, into its storage folder
# Workspace: Gianluca De Novi's workspace · Author: Gianluca De Novi
# Saved: 2026-11-24
#
# Inputs   ctx.inputs["id"] is a list of files: .path .name .size .read_text() .read_json()
#   (none)
# Outputs  ctx.output("id", "name.ext") is the path to write; ctx.write_json / ctx.write_text
#   telemetry   Live telemetry            ctx.telemetry("channel", data) while it runs
# Settings ctx.params["key"]
#   folder      text      Storage folder = "virtual-patient-models"
# Runs on  Any worker. Return a dict to fill the "result" output.
# Packages numpy, torch, scikit-image, scipy, SimpleITK, trimesh, fast-simplification (installed by the worker with pip when missing)
# EsaSim rewrites these lines from the block's definition each time it is saved: change the definition, not them.
# ── end of EsaSim block header ────────────────────────────────────────

"""One-off: on a worker whose run files EsaSim cannot pick up, make the virtual patient's 3D files from the files a
run already left on this machine (the CT case and the block's own code), into a folder of the worker's storage."""
import shutil
from pathlib import Path


class _File:
    def __init__(self, p):
        self.path, self.name = str(p), p.name


class _Ctx:
    def __init__(self, ctx, out, ct):
        self.params, self.inputs, self._ctx, self._out = {"mode": "Make the 3D models", "accelerator": "CPU"}, {"ct": ct}, ctx, out

    def output(self, port, name):
        return str(self._out / name)

    def log(self, *a):
        self._ctx.log(*a)

    def telemetry(self, *a, **k):
        pass

    def progress(self, *a, **k):
        pass


def run(ctx):
    here = Path.cwd().resolve()
    data = next(p for p in here.parents if (p / "jobs").is_dir() and (p / "storage").exists())
    ctx.log("data folder: %s" % data)
    jobs = sorted((p for p in (data / "jobs").iterdir() if p.is_dir() and p.name.isdigit()), key=lambda p: int(p.name), reverse=True)
    ct, code = None, None
    for j in jobs:
        if ct is None:
            hit = list(j.glob("*/outputs/files/s1400/ct.nii.gz"))
            if hit:
                ct = hit[0].parent
        if code is None:
            for f in j.glob("models-*/*.py"):
                t = f.read_text(encoding="utf-8", errors="ignore")
                if "def ct_parts" in t:
                    code = f
                    break
        if ct is not None and code is not None:
            break
    if ct is None or code is None:
        raise RuntimeError("Not found on this machine: CT case %s, block code %s" % (ct, code))
    ctx.log("CT case: %s\nblock code: %s" % (ct, code))
    files = [_File(ct / "ct.nii.gz")] + [_File(p) for p in sorted((ct / "segmentations").glob("*.nii.gz"))]
    out = data / "storage" / str(ctx.params.get("folder") or "virtual-patient-models")
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    ns = {"__name__": "pbpk"}
    exec(compile(code.read_text(encoding="utf-8"), str(code), "exec"), ns)
    res = ns["run"](_Ctx(ctx, out, files))
    ctx.log("written to %s: %s" % (out, ", ".join(sorted(p.name for p in out.iterdir()))))
    return res
