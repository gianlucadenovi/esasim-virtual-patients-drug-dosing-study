# Temporary: patient models into worker storage

One-off: makes the virtual patient's 3D files from a CT case already on this worker, into its storage folder

An [EsaSim](https://esacode.com/) blueprint block: it runs on a worker, in Python (`backend.py`, `def run(ctx)`).

| | |
|---|---|
| Block | `ws-temporary-patient-models-into-worker` · version 1 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |
| Runs on | Any worker |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `folder` Storage folder | text | `"virtual-patient-models"` |  |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
