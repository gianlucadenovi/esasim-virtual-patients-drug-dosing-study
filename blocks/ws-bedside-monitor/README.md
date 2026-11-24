# Bedside monitor

A bedside patient monitor (ECG, pulse oximetry and breathing waveforms, vitals, alarms, beep) with controls to choose the patient and the drug and dose them live

Wire a Telemetry block (or a block's Live telemetry) in. Items on the channel (default vitals): hr, sys, dia, map, rr, spo2, temp, rhythm, ectopy (share of ectopic beats), dropped (share of dropped beats), sedation, effect, alarms [..], plasma, unit, state, target, drug, who, regimen, day, speed, patients [..], drugs [..], dose_x, stopped. The controls send {patient}, {drug}, {dose_x}, {give}, {stop}, {restart}, {speed} (each with a seq number) on the control channel (default bedside) to the running blocks. The waveforms are drawn in real time at the current rates.

An [EsaSim](https://esacode.com/) blueprint block: a view that draws in the browser, in JavaScript (`frontend.js`, `function render(ctx)`).

| | |
|---|---|
| Block | `ws-bedside-monitor` · version 1 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |

## Inputs

| Port | Type | Description |
|---|---|---|
| `stream` Live telemetry | Stream | The vitals |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `channel` Vitals channel | text | `"vitals"` |  |
| `control_channel` Controls channel | text | `"bedside"` |  |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
