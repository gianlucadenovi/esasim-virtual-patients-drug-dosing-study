# Bedside monitor

A bedside patient monitor: ECG, pulse oximetry and breathing waveforms, vitals (HR, SpO2, resp, BP, temperature, QTc, pain, sedation), alarms and the drugs on board

Wire a Telemetry block (or a block's Live telemetry) in. Items on the channel (default vitals): hr, sys, dia, map, rr, spo2, temp, qtc, pain, rhythm, ectopy (share of ectopic beats), dropped (share of dropped beats), sedation, sedation_txt, nausea, airway, effects, alarms [..], who, day, drugs [..]. The waveforms are drawn in real time at the current rates; a long QTc moves the T wave.

An [EsaSim](https://esacode.com/) blueprint block: a view that draws in the browser, in JavaScript (`frontend.js`, `function render(ctx)`).

| | |
|---|---|
| Block | `ws-bedside-monitor` · version 3 |
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

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
