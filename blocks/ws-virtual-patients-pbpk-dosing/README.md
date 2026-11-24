# Virtual patients (PBPK dosing)

A whole-body physiologically based pharmacokinetic model of a drug in thousands of virtual patients: a trial of a dosing regimen, the best dose per population, and one patient followed live in 3D

An [EsaSim](https://esacode.com/) blueprint block: it runs on a worker, in Python (`backend.py`, `def run(ctx)`).

| | |
|---|---|
| Block | `ws-virtual-patients-pbpk-dosing` · version 8 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |
| Runs on | Any worker |

## Inputs

| Port | Type | Description |
|---|---|---|
| `ct` CT case | Files (nii) | Make the 3D models: a CT case with its labels (ct.nii.gz and segmentations/<structure>.nii.gz, as the TotalSegmentator dataset gives them). Without it, the patient is a mockup. Several connections. Optional. |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
