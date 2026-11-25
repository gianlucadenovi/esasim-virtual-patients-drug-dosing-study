# Anaesthetic trolley

The anaesthetist's controls as a trolley: pick a drug, set how much, choose a push or an infusion, administer; the airway and the oxygen; the job, the clock and the demo

Everything the anaesthesia room can be told, grouped the way it is done rather than as a row of controls: a drug with the dose this patient's weight comes to, given as one push or as an infusion; preoxygenating and the airway; the job, the clock, the demo and starting again. It reads the run's patient channel for the doses, for what is already running and for the line of advice at the top, and sends one command for each action. Use it beside the 3D patient in the anaesthesia flow.

An [EsaSim](https://esacode.com/) blueprint block: a view that draws in the browser, in JavaScript (`frontend.js`, `function render(ctx)`).

| | |
|---|---|
| Block | `ws-anaesthetic-trolley` · version 1 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |

## Inputs

| Port | Type | Description |
|---|---|---|
| `stream` Live telemetry | Stream | The run's telemetry: the patient channel carries the doses, the infusions running and the advice |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `channel` Commands go on | text | `"control"` |  |
| `data_channel` Patient channel | text | `"patient"` |  |
| `to` Only to the block | text | `""` | A block's id on the grid; empty = every running block |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
