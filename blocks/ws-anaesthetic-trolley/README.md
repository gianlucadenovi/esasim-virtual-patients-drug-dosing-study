# Anaesthetic trolley

The anaesthetist's trolley: the briefing first (the job, what you have to achieve, do it yourself or be shown the sequence), then the drugs, the airway and the case

Opens on the briefing: the patient, the job in a sentence, the depth to hold and the steps that have to be ticked off, with two ways to take it on — I will do it, or Show me the right sequence, which lists the protocol and offers to play it. After that it is the place the case is worked from: a drug with the dose this patient's weight comes to, given as one push or as an infusion; preoxygenating and the airway; the job, the clock, the demo and starting again. The line of advice sits at the top. It reads the run's patient channel and sends one command for each action. Use it beside the 3D patient in the anaesthesia flow.

An [EsaSim](https://esacode.com/) blueprint block: a view that draws in the browser, in JavaScript (`frontend.js`, `function render(ctx)`).

| | |
|---|---|
| Block | `ws-anaesthetic-trolley` · version 3 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |

## Inputs

| Port | Type | Description |
|---|---|---|
| `stream` Live telemetry | Stream | The run's telemetry: the patient channel carries the briefing, the doses, the infusions running and the advice |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `channel` Commands go on | text | `"control"` |  |
| `data_channel` Patient channel | text | `"patient"` |  |
| `to` Only to the block | text | `""` | A block's id on the grid; empty = every running block |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
