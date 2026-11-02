# events.tsv files as published in version 1.0.2

Version 1.0.3 rewrote every `sub-*/ieeg/*_events.tsv` so that each file holds one `sequence_start` row per sequence
(epoch) and moved the per-stimulus rows to `derivatives/stimulus_events/` (exact duplicate rows dropped; rows outside
the data listed in `sourcedata/events-removed-2026-10-17/`). So that no published row is lost, the 140 events.tsv files
of version 1.0.2 are kept here byte for byte, at the same relative path (e.g.
`sub-001/ieeg/sub-001_task-fpvsfaces_run-1_events.tsv`). They are not read by BIDS tools.
