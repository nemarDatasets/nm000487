"""o5c lane A2 (OSF): thin wrapper over o5c_A1_common (same helpers); copy_code also ships both common modules."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_A1_common import *  # noqa: F401,F403
from o5c_A1_common import write_epochs, events_json, h5str, h5cellstr, wtsv, wjson, finalize_common, copy_sourcedata, sha256_file, add_electrodes_na, write_run, mat_strings, flag_strings, CH_COLS, scrub_mat_header, R, OPEN500  # noqa: F401


def copy_code(OUT, files):
    import shutil
    os.makedirs(os.path.join(OUT, "code"), exist_ok=True)
    H = os.path.dirname(os.path.abspath(__file__))
    for fn in files + ["o5c_A2_common.py", "o5c_A1_common.py", "b3w10_common.py", "b3w3_common.py"]:
        shutil.copy(os.path.join(H, os.path.basename(fn)), os.path.join(OUT, "code", os.path.basename(fn)))
