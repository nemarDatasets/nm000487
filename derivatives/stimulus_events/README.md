# nm000487 stimulus-level events

Stimulus-level Letswave event rows of nm000487 (Jacques, Jonas, Colnat-Coulbois & Rossion 2025, FPVS faces SEEG),
moved here in version 1.0.3 from the raw `sub-NNN/ieeg/sub-NNN_task-fpvsfaces_run-1_events.tsv`, which now holds one
`sequence_start` row per sequence (epoch) so that `mne_bids.read_epochs_bids` loads each file as its sequences.

- 88 files (the other 52 subjects have no stimulus rows on the file's time axis; see `sourcedata/events-removed-2026-10-17/`
  of the dataset), 80,458 rows. Same file names, columns, line text and order as the raw files before the move;
  `onset` and `sample` are relative to the start of the recording; `sequence` is the epoch (1-based) each row belongs to.
- Only exact duplicate rows (identical in every column; 61,229 rows in 56 files) were dropped.
- `stimulus_events_summary.tsv`: per raw file, number of epochs, rows before, duplicates dropped, rows here.
- Loading example: dataset README, section "Events and epochs".
