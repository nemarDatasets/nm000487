"""nm000487 fix-20261017 (events split): number of epochs per file from the data headers. READ-ONLY.
Runs in a Voyager Job. For each of the 140 runs: fetch the .vhdr (annexed) anonymously from the NEMAR S3 bucket,
check its sha256 against the annex key, parse NumberOfChannels / BinaryFormat / SamplingInterval / DataOrientation,
get the .eeg size with an HTTP HEAD (and compare with the size in the annex key), and compute
  n_samples = eeg_bytes / (n_channels * bytes_per_sample),  n_per = round(EpochLength * sfreq),
  n_epochs  = n_samples / n_per  (must be an integer).
Cross-checks: vmrk 'New Segment' markers (count and positions = 1 + k * n_per) and ieeg.json RecordingDuration.
Usage: headers.py files.json out.json"""
import hashlib, json, re, sys, urllib.request

BASE = "https://nemar.s3.us-east-2.amazonaws.com/nm000487/objects/"
NB = {"IEEE_FLOAT_32": 4, "INT_16": 2, "INT_32": 4, "UINT_16": 2}
files = json.load(open(sys.argv[1]))
out = []
for f in files:
    r = {"file": f["file"]}
    try:
        txt = urllib.request.urlopen(BASE + f["vhdr_key"], timeout=60).read()
        sha = re.search(r"--([0-9a-f]{64})", f["vhdr_key"]).group(1)
        r["vhdr_sha256_ok"] = hashlib.sha256(txt).hexdigest() == sha
        t = txt.decode("utf-8", "replace")
        g = lambda k: (re.search(rf"^{k}=(.*)$", t, re.M) or [None, None])[1]
        nch = int(g("NumberOfChannels")); fmt = g("BinaryFormat").strip(); si = float(g("SamplingInterval"))
        r.update(n_channels=nch, binary_format=fmt, orientation=(g("DataOrientation") or "").strip(),
                 sfreq_vhdr=1e6 / si, datafile=(g("DataFile") or "").strip())
        req = urllib.request.Request(BASE + f["eeg_key"], method="HEAD")
        size = int(urllib.request.urlopen(req, timeout=60).headers["Content-Length"])
        r["eeg_bytes"] = size; r["eeg_size_matches_key"] = size == f["eeg_size_from_key"]
        bps = NB[fmt]
        r["n_samples"] = size // (nch * bps); r["bytes_remainder"] = size % (nch * bps)
        n_per = round(f["EpochLength"] * r["sfreq_vhdr"]); r["n_per"] = n_per
        r["n_epochs"] = r["n_samples"] // n_per; r["samples_remainder"] = r["n_samples"] % n_per
        segs = f["vmrk_segments_1based"]
        r["vmrk_segments"] = len(segs)
        r["vmrk_positions_ok"] = segs == [1 + k * n_per for k in range(len(segs))]
        r["json_epochs"] = f["RecordingDuration"] / f["EpochLength"]
        r["sfreq_matches_json"] = abs(r["sfreq_vhdr"] - f["SamplingFrequency"]) < 1e-6
        r["consistent"] = (r["bytes_remainder"] == 0 and r["samples_remainder"] == 0 and r["vhdr_sha256_ok"]
                           and r["eeg_size_matches_key"] and r["vmrk_segments"] == r["n_epochs"] and r["vmrk_positions_ok"]
                           and abs(r["json_epochs"] - r["n_epochs"]) < 1e-9 and r["sfreq_matches_json"])
        r["ok"] = True
    except Exception as e:  # noqa: BLE001
        r["ok"] = False; r["error"] = f"{type(e).__name__}: {e}"
    print(json.dumps(r), flush=True)
    out.append(r)
json.dump(out, open(sys.argv[2], "w"), indent=0)
print("files", len(out), "ok", sum(r["ok"] for r in out), "consistent", sum(r.get("consistent", False) for r in out),
      "total epochs", sum(r.get("n_epochs", 0) for r in out))
