"""nm000487 version 1.0.3 (branch fix-20261017): one events.tsv row per epoch; stimulus rows to a derivative.

`mne_bids.read_epochs_bids` takes every events.tsv row of an epoched file as the start of one epoch of length
`EpochLength`. This script makes every `sub-*/ieeg/*_events.tsv` hold exactly one `sequence_start` row per sequence
(epoch) of the data file, and moves the stimulus-level Letswave rows to `derivatives/stimulus_events/`:

1. Number of epochs per file = n_samples / round(EpochLength x sfreq), with n_samples from the BrainVision header
   (NumberOfChannels, BinaryFormat) and the .eeg size; cross-checked with the .vmrk 'New Segment' markers and
   ieeg.json RecordingDuration (`headers.json`, computed by `headers.py` in a read-only job; all 140 consistent).
2. Files that hold stimulus rows: exact duplicate rows (identical in every column) are dropped, keeping the first
   occurrence; the remaining rows are written unchanged (same line text, same order, onsets relative to the
   recording) to `derivatives/stimulus_events/<same relative path>`. Each row is checked against the epochs
   (`sequence` <= n_epochs and `sample` inside that epoch); rows that do not map are still kept and counted.
3. Every events.tsv is (re)written as n_epochs rows: onset = (k-1) x EpochLength, duration = EpochLength,
   trial_type = sequence_start, sequence = k, latency_in_sequence = SequenceStartTime, sample = (k-1) x n_per
   (same format as the 52 files that already held such rows since the first 1.0.3 change).
Usage: fix-20261017_split_stimulus_events.py <dataset_dir> <headers.json> [--dry]
"""
import csv
import glob
import json
import os
import sys

root, hjson = sys.argv[1], sys.argv[2]
dry = "--dry" in sys.argv
heads = {h["file"]: h for h in json.load(open(hjson))}
HDR = ["onset", "duration", "trial_type", "sequence", "latency_in_sequence", "sample"]
DER = os.path.join(root, "derivatives", "stimulus_events")
summary = []
for ev in sorted(glob.glob(os.path.join(root, "sub-*", "ieeg", "*_events.tsv"))):
    rel = os.path.relpath(ev, root)
    base = rel[: -len("_events.tsv")]
    h = heads[base]
    assert h["ok"] and h["consistent"], base
    j = json.load(open(os.path.join(root, base + "_ieeg.json"), encoding="utf-8"))
    EL, fs, xs = float(j["EpochLength"]), float(j["SamplingFrequency"]), j.get("SequenceStartTime")
    n_ep, n_per = int(h["n_epochs"]), int(h["n_per"])
    n_times = n_ep * n_per
    lines = open(ev, encoding="utf-8", newline="").read().splitlines(True)
    assert lines[0].rstrip("\r\n").split("\t") == HDR, rel
    body = lines[1:]
    tts = [l.rstrip("\r\n").split("\t")[2] for l in body]
    rec = {"file": rel, "n_epochs": n_ep, "samples_per_epoch": n_per, "rows_before": len(body)}
    if body and all(t == "sequence_start" for t in tts):
        rec.update(kind="sequence_start", exact_duplicates_removed=0, rows_to_derivative=0, rows_unmapped=0)
    else:
        assert "sequence_start" not in tts, rel
        seen, keep = set(), []
        for l in body:
            k = l.rstrip("\r\n")
            if k in seen:
                continue
            seen.add(k)
            keep.append(l if l.endswith("\n") else l + "\n")
        unmapped = 0
        for l in keep:
            f = l.rstrip("\r\n").split("\t")
            s, smp, on = int(float(f[3])), int(float(f[5])), float(f[0])
            if not (1 <= s <= n_ep and (s - 1) * n_per <= smp < s * n_per and 0 <= smp < n_times and 0 <= on < n_times / fs):
                unmapped += 1
        rec.update(kind="stimulus", exact_duplicates_removed=len(body) - len(keep), rows_to_derivative=len(keep),
                   rows_unmapped=unmapped)
        if not dry:
            out = os.path.join(DER, rel)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            open(out, "w", encoding="utf-8", newline="").write(lines[0] if lines[0].endswith("\n") else lines[0] + "\n")
            open(out, "a", encoding="utf-8", newline="").write("".join(keep))
    rows = [lines[0] if lines[0].endswith("\n") else lines[0] + "\n"]
    for k in range(1, n_ep + 1):
        on = (k - 1) * EL
        assert round(on * fs) == (k - 1) * n_per, (rel, k)
        row = {"onset": f"{on:.6f}".rstrip("0").rstrip(".") if on else "0.0", "duration": repr(EL),
               "trial_type": "sequence_start", "sequence": str(k),
               "latency_in_sequence": repr(float(xs)) if xs is not None else "n/a", "sample": str((k - 1) * n_per)}
        rows.append("\t".join(row[c] for c in HDR) + "\n")
    rec["rows_after"] = n_ep
    rec["unchanged"] = "".join(rows) == "".join(lines)
    if not dry and not rec["unchanged"]:
        open(ev, "w", encoding="utf-8", newline="").write("".join(rows))
    summary.append(rec)

st = [r for r in summary if r["kind"] == "stimulus"]
tot = {"files": len(summary), "files_sequence_start_already": len(summary) - len(st), "files_with_stimulus_rows": len(st),
       "stimulus_rows_before": sum(r["rows_before"] for r in st),
       "files_with_exact_duplicates": sum(1 for r in st if r["exact_duplicates_removed"]),
       "exact_duplicates_removed": sum(r["exact_duplicates_removed"] for r in st),
       "rows_to_derivative": sum(r["rows_to_derivative"] for r in st),
       "rows_unmapped": sum(r["rows_unmapped"] for r in st),
       "sequence_start_rows_after": sum(r["rows_after"] for r in summary),
       "events_tsv_changed": sum(1 for r in summary if not r["unchanged"])}
print(json.dumps(tot, indent=1))
if not dry:
    cols = ["file", "kind", "n_epochs", "samples_per_epoch", "rows_before", "exact_duplicates_removed",
            "rows_to_derivative", "rows_unmapped", "rows_after"]
    with open(os.path.join(DER, "stimulus_events_summary.tsv"), "w", encoding="utf-8", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in summary:
            w.writerow([r[c] for c in cols])
