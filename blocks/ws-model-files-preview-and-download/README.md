# Model files (preview and download)

Puts back together the files a block streams on a telemetry channel, shows the 3D model among them and offers every file to download

Wire a block's Live telemetry (or a Telemetry block) in. The block sends each file on the channel as pieces: {"file", "part", "parts", "size", "gz", "b64"} (gzip then base64). The view shows the glb named in Show, and a Download button per file. Pieces missing a while after the last one are asked again with ctx.send(channel, {"missing": {file: [parts]}}), and {"complete": true} when every file is in. Useful when a worker's result files cannot reach EsaSim.

An [EsaSim](https://esacode.com/) blueprint block: a view that draws in the browser, in JavaScript (`frontend.js`, `function render(ctx)`).

| | |
|---|---|
| Block | `ws-model-files-preview-and-download` · version 1 |
| Workspace | Gianluca De Novi's workspace |
| Author | Gianluca De Novi |

## Inputs

| Port | Type | Description |
|---|---|---|
| `stream` Live telemetry | Stream | The files, in pieces, on the channel |

## Settings

| Key | Type | Default | Description |
|---|---|---|---|
| `channel` Channel | text | `"model_files"` |  |
| `show` Show | text | `"patient.glb"` | The glb file to show |

## Use it in EsaSim

In a workspace, open **Projects > Blocks > New block**, choose this kind of block, copy the ports and settings above and paste the code; or ask an AI assistant connected to EsaSim (MCP) to create it from this folder.
