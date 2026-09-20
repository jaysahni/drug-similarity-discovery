# Third-Party Notices

This repository is MIT-licensed (see `LICENSE`). It contains adapted
third-party code, noted here.

## Biomni A1 harness adaptation

The files under `scripts/agent/` are derived in part from
[Biomni](https://github.com/snap-stanford/Biomni), commit
`400c1f366b96a35ca253e13c9b06c5076af41d65`, licensed under the Apache License,
Version 2.0, reaching this repository by way of the NovaKit adaptation in
[jaysahni/cheminformatics-kit](https://github.com/jaysahni/cheminformatics-kit)
(`src/novakit/agent/`).

Retained from that lineage: the `<execute>` / `<solution>` response protocol,
the LangGraph `generate` -> `execute` state machine, the malformed-response
retry, the prompt-based tool retriever, and the local Python/R/Bash executor
with its persistent namespace.

Changed here: the tool surface is this repository's pipeline stages under
`scripts/` rather than a package registry; stage purposes are read from each
script's own module docstring instead of being restated; the NovaKit structure
viewer, credential keychain, and data-lake layers are dropped; the environment
prefix is `AUTOREPURPOSE_AGENT_`; and the system prompt carries this project's
evidence standard and the PROJECT_GOAL.md G8 vocabulary ban.

The Apache License, Version 2.0 is available at
https://www.apache.org/licenses/LICENSE-2.0.
