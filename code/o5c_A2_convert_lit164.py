"""src-literature-164 (OSF 2qzym, CC-BY-4.0; Jacques, Jonas, Colnat-Coulbois, Rossion 2025 eLife reviewed preprint
doi:10.7554/eLife.109640.1 "Concurrent category-selective neural activity across the ventral occipito-temporal cortex supports
a non-hierarchical view of human visual recognition") -> iEEG-BIDS DERIVATIVE (Letswave 5 sequences, minimally preprocessed).

Source: _x/seeg_data/*/P_NNN/P_NNN_ieeg.{lw5,mat}: Letswave 5 pair, header (MAT v5 struct 'header') + v7.3 'data' float32
[sequences x 1 x 1 x 1 x channels x time]. Each sequence = one fast periodic visual stimulation sequence (objects at 6 Hz with a
face every 5th image = 1.2 Hz), cut from 2 s before onset to 2 s after the end, merged per participant (README). chanlocs carry
Talairach coordinates and anatomical labels (main region, sub region, VOTC flag) per contact. Demographics and reference contact
per participant from files/File information/*.txt. Pattern: o5c_C_convert_nancy.py (same lab's Letswave releases, NEMAR nm000403).
Values written unchanged as BrainVision IEEE_FLOAT_32, sequences back-to-back (one New Segment marker per sequence);
Letswave events -> events.tsv (codes with a leading participant-code token are replaced; tokens never printed).
Usage: python o5c_A2_convert_lit164.py <open500 record dir> <bids root>"""
import collections, glob, hashlib, io, json, os, re, shutil, sys
import numpy as np, h5py, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_A2_common import wtsv, wjson, write_vhdr, scrub_mat_header, sha256_file, finalize_common, copy_code

REC, OUT = sys.argv[1], sys.argv[2]
X = os.path.join(REC, "_x", "seeg_data")
INFO = os.path.join(REC, "files", "File information")
SD = os.path.join(OUT, "sourcedata", "osf-2qzym-deidentified")
report = {"files": [], "sourcedata": [], "coords": {}, "privacy": {}}
TOKENS = set()
PREFIX = re.compile(r"^([A-Za-z]{3,}[-_]?\d*)_")
LETTERS = re.compile(r"^[A-Za-z]{3,}$")
SAFE_WORDS = {"end", "start", "stop", "seq", "begin", "pause", "rest", "startblock", "endblock"}
SAFE_PREFIX = {"merged", "seq", "lp", "ds"}
DATE = re.compile(r"^(\d{1,2})-([A-Za-z]{3})-(\d{4})")
PATHRX = re.compile(r"[\\/]|^[A-Za-z]:")


def scrub_code(c):
    m = PREFIX.match(c.strip())
    if m and m.group(1).lower() not in SAFE_PREFIX:
        TOKENS.add(m.group(1))
        return "seq_" + c.strip()[m.end():]
    if LETTERS.match(c) and c.lower() not in SAFE_WORDS:
        TOKENS.add(c)
        return "unk"
    return c


def mask(s):
    return re.sub(r"[A-Za-z]", "a", re.sub(r"\d", "#", s))


def sval(a):
    return str(np.ravel(a)[0]) if isinstance(a, np.ndarray) and a.size else ""


def chtype(l):
    if re.match(r"^ECG", l, re.I): return "ECG"
    if re.match(r"^(PULS|BEAT|SpO2|EMG|EOG)", l, re.I): return "MISC"
    if re.match(r"^MKR", l, re.I): return "TRIG"
    return "SEEG"


# participant tables (cp1252: the reference column uses a typographic apostrophe)
def read_tab(p):
    rows = [l.rstrip("\r\n").split("\t") for l in open(p, encoding="cp1252", errors="replace").read().splitlines() if l.strip()]
    return rows[0], rows[1:]


dh, drows = read_tab(os.path.join(INFO, "140_participants_demographics.txt"))
sh, srows = read_tab(os.path.join(INFO, "140_participants_seeg_info.txt"))
print("demographics header", dh, "n", len(drows)); print("seeg_info header", sh, "n", len(srows))
DEMO = {int(r[0]): r for r in drows}
SINFO = {int(r[0]): r for r in srows}

files = sorted(glob.glob(os.path.join(X, "*", "P_*", "P_*_ieeg.lw5")), key=lambda p: int(re.search(r"P_(\d+)_ieeg", p).group(1)))
assert len(files) == 140, len(files)
scans = {}
TASK = "fpvsfaces"
for n in files:
    pid = int(re.search(r"P_(\d+)_ieeg\.lw5$", n).group(1))
    S = f"sub-{pid:03d}"; D = os.path.join(OUT, S, "ieeg"); os.makedirs(D, exist_ok=True)
    h = sio.loadmat(n, squeeze_me=False, struct_as_record=True)["header"][0, 0]
    ds = [int(v) for v in np.ravel(h["datasize"])]; nch, ns = ds[1], ds[5]
    fs = 1.0 / float(np.ravel(h["xstep"])[0]); x0 = float(np.ravel(h["xstart"])[0])
    C = h["chanlocs"]
    labels = [sval(C.flat[j]["labels"]).replace("’", "'") for j in range(C.size)]
    assert len(labels) == nch, (len(labels), nch)
    with h5py.File(n[:-4] + ".mat", "r") as hf:
        assert not hf["data"].attrs.get("MATLAB_empty", 0), n
        a = hf["data"][()]
    ne = a.shape[-1]
    assert a.shape == (ns, 1, 1, 1, nch, ne), (a.shape, ds)
    assert a.dtype == np.float32
    x = np.concatenate([a[:, 0, 0, 0, :, e] for e in range(ne)], axis=0)
    del a
    stem = f"{S}_task-{TASK}_run-1"
    np.ascontiguousarray(x).astype("<f4").tofile(os.path.join(D, stem + "_ieeg.eeg"))
    ev = []
    E = h["events"]
    if E.dtype.names:
        for j in range(E.size):
            e = E.flat[j]
            ep = int(np.ravel(e["epoch"])[0]); lat = float(np.ravel(e["latency"])[0])
            s0 = int(round((lat - x0) * fs))
            ev.append([round((ep - 1) * ns / fs + (lat - x0), 6), 0, scrub_code(sval(e["code"])), ep, round(lat, 6), (ep - 1) * ns + s0])
    ev.sort(key=lambda r: r[0])
    segs = [("New Segment", f"sequence {k + 1}", k * ns, 1) for k in range(1, ne)]
    write_vhdr(os.path.join(D, stem + "_ieeg"), nch, fs, labels, ["µV"] * nch,
               comment=f"o5c_A2_convert_lit164.py: Letswave {os.path.basename(n)} + .mat, {ne} sequences back-to-back, float32 unchanged", markers=segs)
    wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "sequence", "latency_in_sequence", "sample"], ev)
    typ = [chtype(l) for l in labels]
    # per-contact anatomy from chanlocs
    fn = C.dtype.names
    tal = []; main = []; subr = []; votc = []
    for j in range(C.size):
        c = C.flat[j]
        v = np.ravel(c["xyz_coord_Talairach"]).astype(float) if "xyz_coord_Talairach" in fn and c["xyz_coord_Talairach"].size == 3 else None
        tal.append(v if v is not None and np.all(np.isfinite(v)) else None)
        main.append(sval(c["anatomy_mainRegionLabel"]) if "anatomy_mainRegionLabel" in fn else "")
        subr.append(sval(c["anatomy_subRegionLabel"]) if "anatomy_subRegionLabel" in fn else "")
        votc.append(sval(c["anatomy_VOTC"]) if "anatomy_VOTC" in fn else "")
    si = SINFO.get(pid); de = DEMO.get(pid)
    ref = si[1].replace("’", "'").strip() if si else "n/a"
    first_missing = si[2].strip() if si else "n/a"; seqlen = si[3].strip() if si else "n/a"
    def grp(l):
        return re.sub(r"\d+$", "", l) or "n/a"
    wtsv(os.path.join(D, stem + "_channels.tsv"), ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "reference", "group", "status", "status_description", "description"],
         [[l, t, "µV", "n/a", "n/a", fs, ref if t == "SEEG" else "n/a", grp(l) if t == "SEEG" else "n/a", "good", "n/a",
           (f"intracerebral contact; main region {main[i] or 'n/a'}, sub region {subr[i] or 'n/a'}, {votc[i] or 'n/a'} (Letswave chanlocs)" if t == "SEEG" else "non-SEEG channel (label from the Letswave header)")]
          for i, (l, t) in enumerate(zip(labels, typ))])
    wjson(os.path.join(D, stem + "_ieeg.json"), {
        "TaskName": TASK,
        "TaskDescription": "Fast periodic visual stimulation: natural images of non-face objects at 6 Hz with a face image every 5th item (1.2 Hz); each sequence starts with 4 object images before the first face [O O O O F O O O O F ...]. Sequences of "
                           f"{seqlen} s (release table) cut from 2 s before onset to 2 s after the end and merged per participant (README). Participants: " + ("the first trigger of each sequence is missing and header.xstart was corrected to -1.8333 s" if first_missing == "1" else "first trigger present, xstart -2 s") + ".",
        "Instructions": "n/a (see article)",
        "SamplingFrequency": fs, "PowerLineFrequency": 50,
        "SoftwareFilters": {"Letswave5": {"Description": "Minimal preprocessing: Micromed TRC imported in Letswave 5, segmented into sequences and merged; no filtering stated (README)."}},
        "HardwareFilters": "n/a",
        "iEEGReference": f"{ref} (reference contact of this participant from 140_participants_seeg_info.txt)",
        "Manufacturer": "Micromed",
        "iEEGPlacementScheme": "Stereotactic depth electrodes (SEEG) implanted for presurgical epilepsy evaluation (University Hospital of Nancy); contact labels with a prime (') are in the left hemisphere (README). Talairach coordinates and anatomical labels per contact in electrodes.tsv.",
        "RecordingType": "epoched", "EpochLength": ns / fs, "RecordingDuration": ne * ns / fs,
        "SEEGChannelCount": typ.count("SEEG"), "ECOGChannelCount": 0, "EEGChannelCount": 0, "ECGChannelCount": typ.count("ECG"),
        "MiscChannelCount": typ.count("MISC"), "TriggerChannelCount": typ.count("TRIG"),
        "ElectricalStimulation": False,
        "SequenceStartTime": x0,
    })
    # electrodes.tsv (Talairach from chanlocs) + anatomy columns
    erows = []
    nt = 0
    for i, l in enumerate(labels):
        if typ[i] != "SEEG":
            continue
        v = tal[i]
        if v is not None:
            nt += 1
        erows.append([l, *(["%.3f" % c for c in v] if v is not None else ["n/a"] * 3), "n/a", grp(l), main[i] or "n/a", subr[i] or "n/a", votc[i] or "n/a"])
    wtsv(os.path.join(D, f"{S}_space-Talairach_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "anatomy_main_region", "anatomy_sub_region", "votc"], erows)
    wjson(os.path.join(D, f"{S}_space-Talairach_electrodes.json"), {"SpatialReference": "Talairach atlas space (authors' coordinates)",
        "anatomy_main_region": {"Description": "Authors' anatomical label of the main region of the contact (Letswave chanlocs anatomy_mainRegionLabel; e.g. OCC, PTL, ATL, MTL, TP; WM = white matter; empty in the release -> n/a)"},
        "anatomy_sub_region": {"Description": "Authors' anatomical sub-region label of the contact (chanlocs anatomy_subRegionLabel; 'no_info' as released)"},
        "votc": {"Description": "Whether the contact is within the ventral occipito-temporal cortex region as defined by the authors (chanlocs anatomy_VOTC: VOTC / not_VOTC)"}})
    wjson(os.path.join(D, f"{S}_space-Talairach_coordsystem.json"), {
        "iEEGCoordinateSystem": "Talairach", "iEEGCoordinateUnits": "mm",
        "iEEGCoordinateSystemDescription": "3D Talairach coordinates of each intracerebral contact as stored by the authors in the Letswave header (chanlocs.xyz_coord_Talairach; README).",
        "iEEGCoordinateProcessingDescription": "Copied from the Letswave header; no processing. Contacts without finite coordinates have x/y/z n/a."})
    scans.setdefault(S, []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a"])
    import mne
    rr = mne.io.read_raw_brainvision(os.path.join(D, stem + "_ieeg.vhdr"), preload=False, verbose="error")
    k = min(2000, x.shape[0])
    rt = bool(rr.n_times == x.shape[0] and rr.ch_names == labels and np.allclose(rr.get_data(start=0, stop=k) * 1e6, x[:k].T, rtol=1e-6, atol=1e-6)
              and np.allclose(rr.get_data(start=x.shape[0] - k) * 1e6, x[-k:].T, rtol=1e-6, atol=1e-6))
    report["files"].append(dict(file=os.path.relpath(n, REC), sub=S, pid=pid, sequences=ne, samples_per_sequence=ns, sfreq=fs, xstart=x0, nch=nch, n_events=len(ev),
                                types={t: typ.count(t) for t in set(typ)}, roundtrip_ok=rt, n_tal=nt, reference=ref, first_trigger_missing=first_missing, seq_len=seqlen,
                                main_regions=dict(collections.Counter(m or "n/a" for m in main)), votc=dict(collections.Counter(v or "n/a" for v in votc))))
    report["coords"][S] = {"n_seeg": typ.count("SEEG"), "n_talairach": nt}
    print(stem, ne, nch, ns, fs, x0, "events", len(ev), "tal", nt, "rt", rt, flush=True)
for S, r in scans.items():
    wtsv(os.path.join(OUT, S, f"{S}_scans.tsv"), ["filename", "acq_time"], sorted(r))
subs = sorted(scans)
prow = []
for f in report["files"]:
    de = DEMO.get(f["pid"])
    sex = {"F": "F", "M": "M"}.get(de[1].strip().upper(), "n/a") if de else "n/a"
    hand = {"R": "R", "L": "L", "A": "A"}.get(de[3].strip().upper(), "n/a") if de else "n/a"
    age = de[5].strip() if de and de[5].strip() else "n/a"
    prow.append([f["sub"], age, sex, hand, f"P_{f['pid']:03d}", f["reference"], f["first_trigger_missing"], f["seq_len"], f["sfreq"], f["sequences"], f["types"].get("SEEG", 0), f["n_tal"]])
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "handedness", "source_id", "reference_contact", "first_trigger_missing", "sequence_length", "sampling_frequency", "n_sequences", "n_seeg_channels", "n_contacts_with_coordinates"], prow)
wjson(os.path.join(OUT, "participants.json"), {
    "age": {"Description": "Age at SEEG (140_participants_demographics.txt)", "Units": "year"},
    "sex": {"Description": "Sex (140_participants_demographics.txt)", "Levels": {"F": "female", "M": "male"}},
    "handedness": {"Description": "Handedness (140_participants_demographics.txt)", "Levels": {"R": "right", "L": "left", "A": "ambidextrous"}},
    "source_id": {"Description": "Participant folder in the OSF release"},
    "reference_contact": {"Description": "Reference contact during the recording (140_participants_seeg_info.txt)"},
    "first_trigger_missing": {"Description": "1 if the first trigger of each sequence was missing (header.xstart corrected to -1.8333 s), 0 otherwise (140_participants_seeg_info.txt)"},
    "sequence_length": {"Description": "Length of each stimulation sequence (140_participants_seeg_info.txt)", "Units": "s"},
    "sampling_frequency": {"Description": "Sampling rate of the recording", "Units": "Hz"},
    "n_sequences": {"Description": "Number of FPVS sequences in the participant's file"},
    "n_seeg_channels": {"Description": "Channels typed SEEG"},
    "n_contacts_with_coordinates": {"Description": "SEEG contacts with finite Talairach coordinates in the Letswave header"}})
wjson(os.path.join(OUT, f"task-{TASK}_events.json"), {
    "onset": {"Description": "Event time in this file (sequences back-to-back)", "Units": "s"},
    "duration": {"Description": "0 (instantaneous Letswave event)", "Units": "s"},
    "trial_type": {"Description": "Letswave event code as stored (numeric Micromed trigger codes; a leading participant-code token before '_' is replaced by 'seq'; other free-text codes by 'unk'). The README states the event codes are meaningless for most participants; stimulus onsets are periodic (6 Hz, faces at 1.2 Hz)."},
    "sequence": {"Description": "Sequence (Letswave epoch) number, 1-based"},
    "latency_in_sequence": {"Description": "Letswave event latency on the sequence time axis (header xstart = -2 or -1.8333 s)", "Units": "s"},
    "sample": {"Description": "Sample index in this file (0-based)"}})


def scrub_paths(a):
    n = 0
    if not isinstance(a, np.ndarray):
        return 0
    if a.dtype.names:
        for i in range(a.size):
            for f in a.dtype.names:
                v = a.flat[i][f]
                if isinstance(v, np.ndarray) and v.dtype.kind in "US":
                    s = sval(v)
                    if PATHRX.search(s):
                        a.flat[i][f] = np.array(["path removed"]); n += 1
                else:
                    n += scrub_paths(v)
    elif a.dtype == object:
        for i in range(a.size):
            v = a.flat[i]
            if isinstance(v, np.ndarray) and v.dtype.kind in "US":
                if PATHRX.search(sval(v)):
                    a.flat[i] = np.array(["path removed"]); n += 1
            else:
                n += scrub_paths(v)
    return n


def scrub_lw5(b):
    m = sio.loadmat(io.BytesIO(b), squeeze_me=False, struct_as_record=True)
    H = m["header"][0, 0]
    nchg = 0
    if "events" in H.dtype.names and H["events"].size and H["events"].dtype.names:
        E = H["events"]
        for i in range(E.size):
            c = sval(E.flat[i]["code"])
            nc = scrub_code(c)
            if nc != c:
                E.flat[i]["code"] = np.array([nc]); nchg += 1
    if "history" in H.dtype.names and H["history"].dtype.names:
        Hs = H["history"]
        for i in range(Hs.size):
            if "date" in Hs.dtype.names:
                s = sval(Hs.flat[i]["date"])
                ns_ = DATE.sub(lambda mm: f"01-{mm.group(2)}-{mm.group(3)}", s)
                if ns_ != s:
                    Hs.flat[i]["date"] = np.array([ns_]); nchg += 1
        nchg += scrub_paths(Hs)
    if "name" in H.dtype.names and H["name"].size:
        H["name"] = np.array(["deidentified"]); nchg += 1
    bio = io.BytesIO()
    sio.savemat(bio, {"header": m["header"]}, format="5", do_compression=True)
    nb, _ = scrub_mat_header(bio.getvalue())
    return nb, nchg


os.makedirs(SD, exist_ok=True)
for n in files:
    for p in (n, n[:-4] + ".mat"):
        rel = os.path.relpath(p, X); rel2 = re.sub(r"^[^/]+\.zip\.d/", "", rel)
        b = open(p, "rb").read(); o = hashlib.sha256(b).hexdigest()
        if p.endswith(".lw5"):
            nb, k = scrub_lw5(b); note = f"Letswave header rewritten (scipy savemat): {k} fields changed (event-code tokens, name, history dates/paths)"
        else:
            nb, ch = scrub_mat_header(b); note = "MAT text header date -> Mmm 01 yyyy" if ch else "unchanged"
        t = os.path.join(SD, rel2); os.makedirs(os.path.dirname(t), exist_ok=True)
        open(t, "wb").write(nb)
        report["sourcedata"].append([rel2, len(b), o, hashlib.sha256(nb).hexdigest(), note])
xl_strings = []
for f in sorted(os.listdir(INFO)):
    p = os.path.join(INFO, f)
    if f.lower().endswith(".xlsx"):
        # not copied: its content could not be inspected in the conversion environment (no openpyxl); the two .txt tables carry the same information
        report.setdefault("not_copied", []).append(f)
        continue
    shutil.copy(p, os.path.join(SD, f))
    report["sourcedata"].append([f, os.path.getsize(p), sha256_file(p), sha256_file(os.path.join(SD, f)), "unchanged"])
report["xlsx_strings_sample"] = sorted(set(xl_strings))[:80]
wtsv(os.path.join(SD, "DEIDENTIFICATION_MANIFEST.tsv"), ["path", "bytes_original", "sha256_original", "sha256_here", "change"], report["sourcedata"])
# token search over the output
hits = 0; HITLOG = []
# whole-word search; only code-like tokens (letters+digits, as in the PREFIX rule) are fatal, pure words are reported
CODE_TOKENS = sorted(t for t in TOKENS if re.search(r"\d", t))
WORD_TOKENS = sorted(t for t in TOKENS if not re.search(r"\d", t))
TB = [(t, re.compile(rb"(?<![A-Za-z0-9])" + re.escape(t.encode()) + rb"(?![A-Za-z0-9])")) for t in CODE_TOKENS + WORD_TOKENS]
word_hits = 0
print("tokens (masked):", sorted({mask(t) for t in TOKENS}), "code-like:", len(CODE_TOKENS), "words:", len(WORD_TOKENS))
if TB:
    for root, _, fs_ in os.walk(OUT):
        for fn in fs_:
            p = os.path.join(root, fn)
            if fn.endswith(".lw5"):
                hh = sio.loadmat(p, squeeze_me=False, struct_as_record=True)["header"][0, 0]
                E = hh["events"]
                txt = repr([sval(E.flat[j]["code"]) for j in range(E.size)]).encode() if E.dtype.names else b""
                txt += repr(hh["name"]).encode() + repr(hh["history"]).encode()
            elif fn.endswith((".eeg", ".mat", ".xlsx")):
                continue
            else:
                txt = open(p, "rb").read()
            for t, rx in TB:
                if rx.search(txt):
                    if re.search(r"\d", t):
                        hits += 1
                    else:
                        word_hits += 1
                    if len(HITLOG) < 40:
                        HITLOG.append((os.path.relpath(p, OUT), mask(t), len(t), "code" if re.search(r"\d", t) else "word"))
print("token hits (file, masked token, len):", HITLOG)
cc_ = collections.Counter()
for fe in glob.glob(os.path.join(OUT, "sub-*", "ieeg", "*_events.tsv")):
    for line in open(fe).read().splitlines()[1:]:
        cc_[line.split("\t")[2]] += 1
report["event_codes_after_scrub"] = dict(cc_)
report["privacy"] = {"patient_code_tokens_removed": len(TOKENS), "code_like_tokens": len(CODE_TOKENS), "token_hits_after_scrub": hits, "word_token_hits_after_scrub": word_hits}
copy_code(OUT, [__file__])
F = report["files"]
hrs = sum(f["sequences"] * f["samples_per_sequence"] / f["sfreq"] for f in F) / 3600
ntal = sum(f["n_tal"] for f in F); nseeg = sum(f["types"].get("SEEG", 0) for f in F)
fsc = collections.Counter(f["sfreq"] for f in F)
finalize_common(OUT,
    name="s500-literature-164: Face-selective SEEG across the ventral occipito-temporal cortex in 140 participants, fast periodic visual stimulation (Jacques et al., 2025) (derivative)",
    dtype="derivative", license_="CC-BY-4.0",
    authors=["Corentin Jacques", "Jacques Jonas", "Sophie Colnat-Coulbois", "Bruno Rossion"],
    refs=["https://doi.org/10.7554/eLife.109640.1", "https://doi.org/10.17605/OSF.IO/2QZYM", "https://osf.io/2qzym/", "http://www.nocions.org/letswave/"],
    howto="Cite Jacques C, Jonas J, Colnat-Coulbois S, Rossion B (2025) Concurrent category-selective neural activity across the ventral occipito-temporal cortex supports a non-hierarchical view of human visual recognition. eLife reviewed preprint, doi:10.7554/eLife.109640.1; and the OSF dataset doi:10.17605/OSF.IO/2QZYM.",
    keywords=["SEEG", "intracerebral recordings", "face recognition", "fast periodic visual stimulation", "frequency tagging", "ventral occipito-temporal cortex", "high-frequency broadband", "epilepsy", "Talairach"],
    source=[{"DOI": "doi:10.17605/OSF.IO/2QZYM", "URL": "https://osf.io/2qzym/", "Version": "OSF node 2qzym, files as listed 2026-10-11 (5 zips seeg_data_p*.zip + File information/; sha256 per file in sourcedata/osf-2qzym-deidentified/DEIDENTIFICATION_MANIFEST.tsv)"}],
    generated_by=[{"Name": "Authors' preprocessing in Letswave 5", "Description": "Micromed TRC imported in Letswave 5, each recording segmented into FPVS sequences (2 s before onset to 2 s after the end) and merged per participant; no further transform (README)."},
                  {"Name": "o5c_A2_convert_lit164.py", "Description": "Repackaging: Letswave float32 sequences back-to-back in BrainVision (values unchanged, MNE read-back checked); Letswave events to events.tsv; Talairach coordinates and anatomical labels from the Letswave chanlocs to electrodes.tsv; demographics and reference contact from the release tables. Code in code/."}],
    ethics=["Collected under the source study's approval (University Hospital of Nancy, presurgical SEEG evaluation; see the article doi:10.7554/eLife.109640.1)."],
    readme=f"""# Jacques, Jonas, Colnat-Coulbois & Rossion (2025): face-selective SEEG across the ventral occipito-temporal cortex (140 participants)

Intracerebral EEG (SEEG) from {len(subs)} patients with refractory epilepsy (University Hospital of Nancy) viewing fast periodic
visual stimulation: natural images of objects at 6 Hz with a face every 5th image (1.2 Hz). {len(F)} files, {hrs:.2f} h of
sequences, {nseeg} SEEG channels in total ({ntal} with Talairach coordinates). Sampling rate {dict(fsc)} (files per rate).
Released by the authors on OSF (https://osf.io/2qzym/, CC-BY 4.0) with the Letswave loader and participant tables.

THIS IS A DERIVATIVE DATASET: the release holds the authors' minimally preprocessed Letswave 5 files (segmented into
stimulation sequences and merged per participant; no filtering stated), not the raw Micromed recordings.

## Contents
- `sub-NNN/ieeg/sub-NNN_task-{TASK}_run-1_ieeg.*`: the participant's FPVS sequences back-to-back (one BrainVision segment per
  sequence; `RecordingType` epoched, `EpochLength` = sequence length incl. the 2 s before and after). Values are the release's
  float32 values, unchanged. The sequence time axis starts at `xstart` (-2 s, or -1.8333 s for participants whose first trigger
  was missing; `participants.tsv` first_trigger_missing, `ieeg.json` SequenceStartTime).
- `channels.tsv`: labels from the Letswave header (a prime ' marks the left hemisphere); all contacts typed SEEG; unit µV
  assumed (Letswave import of Micromed TRC; not stated in the release); reference = the participant's reference contact from
  `140_participants_seeg_info.txt`.
- `sub-NNN_space-Talairach_electrodes.tsv`: Talairach x/y/z per contact from the Letswave header (chanlocs.xyz_coord_Talairach)
  plus the authors' anatomical labels: `anatomy_main_region` (OCC, PTL, ATL, MTL, TP, WM...), `anatomy_sub_region` and `votc`
  (inside the authors' VOTC region or not). Contacts without finite coordinates have n/a.
- `events.tsv`: Letswave events (trigger codes, sequence, latency); the README says the codes are meaningless for most
  participants (stimulus onsets are periodic).
- `participants.tsv`: age at SEEG, sex, handedness (release demographics), reference contact, first-trigger flag, sequence
  length, sampling rate, number of sequences and of channels.
- `sourcedata/osf-2qzym-deidentified/`: the release's Letswave files per participant (headers rewritten, see Privacy), the README,
  loader and the two participant .txt tables (`140_participants_info.xlsx` is not copied: it could not be inspected in the
  conversion environment; the .txt tables carry the same information); `DEIDENTIFICATION_MANIFEST.tsv` lists original and
  new sha-256 per file.

## Privacy
- Letswave event codes with a leading participant-code token were replaced ('seq_...'; other free-text codes 'unk') in
  events.tsv and in the de-identified .lw5 headers ({report['privacy']['patient_code_tokens_removed']} distinct tokens; a byte
  search of the output found {report['privacy']['token_hits_after_scrub']} remaining occurrences). Header `name` fields were
  blanked, history dates reduced to month/year and file paths removed; MAT text-header dates reduced to month/year.
- Participants are identified by number only; demographics are age, sex and handedness as released.
""",
    changes_note="BIDS (derivative) repackaging of the OSF 2qzym release (Jacques et al. 2025) for NEMAR (o5c lane A2).")
json.dump(report, open(os.path.join(OUT, "code", "conversion_report.json"), "w"), indent=1, default=str)
print(json.dumps({"files": len(F), "subjects": len(subs), "hours": hrs, "rt_all": all(f["roundtrip_ok"] for f in F), "privacy": report["privacy"], "n_tal": ntal, "n_seeg": nseeg,
                  "xlsx_strings_sample": report["xlsx_strings_sample"][:40]}, indent=1))
assert hits == 0, "patient-code tokens remain"
assert all(f["roundtrip_ok"] for f in F)
print("CONVERT_DONE")
