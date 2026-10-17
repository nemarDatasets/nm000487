[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000487-blue)](https://doi.org/10.82901/nemar.nm000487)

# Jacques, Jonas, Colnat-Coulbois & Rossion (2025): face-selective SEEG across the ventral occipito-temporal cortex (140 participants)

Intracerebral EEG (SEEG) from 140 patients with refractory epilepsy (University Hospital of Nancy) viewing fast periodic
visual stimulation: natural images of objects at 6 Hz with a face every 5th image (1.2 Hz). 140 files, 9.67 h of
sequences, 22160 SEEG channels in total (20495 with Talairach coordinates). Sampling rate {512.0: 126, 256.0: 3, 500.0: 11} (files per rate).
Released by the authors on OSF (https://osf.io/2qzym/, CC-BY 4.0) with the Letswave loader and participant tables.

THIS IS A DERIVATIVE DATASET: the release holds the authors' minimally preprocessed Letswave 5 files (segmented into
stimulation sequences and merged per participant; no filtering stated), not the raw Micromed recordings.

## Contents
- `sub-NNN/ieeg/sub-NNN_task-fpvsfaces_run-1_ieeg.*`: the participant's FPVS sequences back-to-back (one BrainVision segment per
  sequence; `RecordingType` epoched, `EpochLength` = sequence length incl. the 2 s before and after). Values are the release's
  float32 values, unchanged. The sequence time axis starts at `xstart` (-2 s, or -1.8333 s for participants whose first trigger
  was missing; `participants.tsv` first_trigger_missing, `ieeg.json` SequenceStartTime).
- `channels.tsv`: labels from the Letswave header (a prime ' marks the left hemisphere); all contacts typed SEEG; unit µV
  assumed (Letswave import of Micromed TRC; not stated in the release); reference = the participant's reference contact from
  `140_participants_seeg_info.txt`.
- `sub-NNN_space-Talairach_electrodes.tsv`: Talairach x/y/z per contact from the Letswave header (chanlocs.xyz_coord_Talairach)
  plus the authors' anatomical labels: `anatomy_main_region` (OCC, PTL, ATL, MTL, TP, WM...), `anatomy_sub_region` and `votc`
  (inside the authors' VOTC region or not). Contacts without finite coordinates have n/a.
- `events.tsv`: since version 1.0.3, exactly one `sequence_start` row per sequence (epoch) of the data file, so that
  `mne_bids.read_epochs_bids` loads each file as its sequences (see "Events and epochs").
- `derivatives/stimulus_events/`: the stimulus-level Letswave events (trigger codes, sequence, latency) of the 88 files
  that have them, with onsets in seconds from the start of the recording; the release README says the codes are
  meaningless for most participants (stimulus onsets are periodic). Events that are not on the file's time axis were
  removed in 1.0.3 (see "Events outside the data").
- `participants.tsv`: age at SEEG, sex, handedness (release demographics), reference contact, first-trigger flag, sequence
  length, sampling rate, number of sequences and of channels.
- `sourcedata/osf-2qzym-deidentified/`: the release's Letswave files per participant (headers rewritten, see Privacy), the README,
  loader and the two participant .txt tables (`140_participants_info.xlsx` is not copied: it could not be inspected in the
  conversion environment; the .txt tables carry the same information); `DEIDENTIFICATION_MANIFEST.tsv` lists original and
  new sha-256 per file.

## Privacy
- Letswave event codes with a leading participant-code token were replaced ('seq_...'; other free-text codes 'unk') in
  events.tsv and in the de-identified .lw5 headers (108 distinct tokens; a byte
  search of the output found 0 remaining occurrences). Header `name` fields were
  blanked, history dates reduced to month/year and file paths removed; MAT text-header dates reduced to month/year.
- Participants are identified by number only; demographics are age, sex and handedness as released.

## Events outside the data (corrected 2026-10-17, version 1.0.3)

Up to version 1.0.2, 120 of the 140 `events.tsv` files held Letswave events that are not on the time axis of the data
file: 88 files had onsets after the end of the recording (for example 684 s in a 312 s file), and many more rows fell
inside the file but outside the sequence they are listed in. Each file holds `RecordingDuration` / `EpochLength`
sequences back to back, so the position of a row within its own sequence is `onset` - (`sequence` - 1) x `EpochLength`,
and the sequence covers `latency_in_sequence` from `SequenceStartTime` (-2 s or -1.83 s) to the end of the epoch. Two
rules were applied (35,685 rows removed in total; no row was moved or edited):

- **Outside its own sequence** (68 files, 14,679 rows): a row whose position within its sequence is negative or not
  shorter than `EpochLength` (it lies after the end of its sequence or of the file).
- **Event block on another time axis** (52 files, 21,006 rows; all rows of the file): the file's events form one block
  of about 70 s (`startBlock` ... `endBlock`, ~420 rows), all listed as sequence 1, whose first latency is 9-407 s
  instead of about 0 s and which runs 3-401 s past the end of sequence 1. Because the sequences were cut from 2 s before
  their onset (release README), a block that starts this late and does not fit inside its sequence was recorded on
  another clock (most likely the original continuous recording); its true position in the file cannot be recovered
  from the release. These files are: sub-009, sub-011, sub-012, sub-014, sub-015, sub-018, sub-019, sub-021, sub-023, sub-026, sub-027, sub-028, sub-032, sub-033, sub-035, sub-038, sub-039, sub-040, sub-042, sub-043, sub-044, sub-045, sub-046, sub-047, sub-048, sub-049, sub-050, sub-051, sub-053, sub-054, sub-055, sub-056, sub-057, sub-058, sub-059, sub-061, sub-062, sub-063, sub-064, sub-113, sub-119, sub-123, sub-124, sub-129, sub-131, sub-132, sub-133, sub-135, sub-136, sub-137, sub-138, sub-140.

The 52 files of the second group would otherwise have an empty `events.tsv`. They hold one row per sequence with
`trial_type` `sequence_start` (`onset` = (k - 1) x `EpochLength`, `duration` = `EpochLength`, `latency_in_sequence` =
`SequenceStartTime`, `sample` = onset x sampling rate), which marks the start of each epoch exactly as the file
structure defines it; no stimulus timing is implied. The stimulation is periodic (6 Hz, faces every 5th image, release
README), so stimulus times within a sequence follow from its onset. All removed rows are kept, with the file name and the
rule, in `sourcedata/events-removed-2026-10-17/removed_events.tsv` (per-file counts in `removed_events_summary.json`). The
other 20 files were unchanged by these two rules. The stimulus rows that remained in 88 files were then moved to
`derivatives/stimulus_events/` (next section).

## Events and epochs (version 1.0.3)

Each data file is epoched (`RecordingType` epoched): its sequences are stored back to back, each `EpochLength` long.
`mne_bids.read_epochs_bids` cuts one epoch of `EpochLength` at every row of `events.tsv`, so `events.tsv` now holds
exactly one row per sequence and nothing else: `trial_type` `sequence_start`, `onset` = (k - 1) x `EpochLength`,
`duration` = `EpochLength`, `sequence` = k, `latency_in_sequence` = `SequenceStartTime`, `sample` = (k - 1) x
round(`EpochLength` x sampling rate). The number of sequences of every file was taken from its data header (number of
samples from the BrainVision header and .eeg size, divided by the samples per epoch) and agrees with the .vmrk
`New Segment` markers and `RecordingDuration` / `EpochLength` in `ieeg.json` (452 sequences in 140 files;
`code/fix-20261017_epochs_from_headers.json`).

The stimulus-level Letswave rows of the 88 files that held them are in
`derivatives/stimulus_events/sub-NNN/ieeg/sub-NNN_task-fpvsfaces_run-1_events.tsv` (same file name as the raw
`events.tsv`), unchanged: same columns, same order, `onset` in seconds from the start of the recording, `sequence` =
the epoch the row belongs to (all 80,458 rows lie inside their epoch). Only exact duplicate rows (identical in every
column, 61,229 rows in 56 files; for example every `s1` row of sub-105 appeared twice) were dropped, so no information
is lost. Column descriptions: `derivatives/stimulus_events/task-fpvsfaces_events.json`; per-file counts:
`derivatives/stimulus_events/stimulus_events_summary.tsv`; code: `code/fix-20261017_split_stimulus_events.py`. The
other 52 files have no stimulus rows on the file's time axis (previous section; their original rows are in
`sourcedata/events-removed-2026-10-17/`).

```python
import mne
import pandas as pd
from mne_bids import BIDSPath, read_epochs_bids

root = "nm000487"  # path of the dataset
bp = BIDSPath(root=root, subject="105", task="fpvsfaces", run="1", datatype="ieeg",
              suffix="ieeg", extension=".vhdr")
epochs = read_epochs_bids(bp)  # 4 epochs (one per sequence_start row) x 78.0 s

# stimulus-level events of the same file (onsets in s from the start of the recording)
stim = pd.read_csv(f"{root}/derivatives/stimulus_events/sub-105/ieeg/"
                   "sub-105_task-fpvsfaces_run-1_events.tsv", sep="\t")
n_per = len(epochs.times)                      # samples per epoch
stim["epoch"] = stim["sequence"] - 1           # index into epochs
stim["time_in_epoch"] = (stim["sample"] - stim["epoch"] * n_per) / epochs.info["sfreq"]  # on epochs.times

# or as annotations of the continuous file
raw = mne.io.read_raw_brainvision(bp.fpath)
raw.set_annotations(mne.Annotations(stim["onset"], stim["duration"], stim["trial_type"].astype(str)))
```
