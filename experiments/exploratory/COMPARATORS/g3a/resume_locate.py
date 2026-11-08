"""Resume the locate stage of repro_quakemigrate.py after the process was killed at event 48454 of 49249.

The triggered events from FIRST_ID onward (in the order locate reads them: the day files concatenated) are written to
a trigger file and located with QuakeScan.locate(trigger_file=...), the scan configured exactly as
repro_quakemigrate.locate does (same onset, picker, marginal window, threads). Each event is located independently,
so this gives the events the uninterrupted run would have given. Events already written are not touched.
"""
import sys
from pathlib import Path

import pandas as pd

REPO = r"C:\Users\juanj\OneDrive\Escritorio\blog\timing-tolerant-mt-detection\experiments\exploratory\COMPARATORS"
WORK = r"C:\Users\juanj\ttmt\qm"
FIRST_ID = 20090121231648784          # "EVENT - 48454 of 49249" in logs\locate.log, interrupted during picking
sys.path.insert(0, REPO)
import repro_quakemigrate as RQ  # noqa: E402


def main():
    RQ._single_thread()
    run = RQ._paths(WORK)["runs"] / RQ.RUN
    files = sorted((run / "trigger" / "events").glob("*_TriggeredEvents.csv"))
    trig = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    start = int(trig.index[trig["EventID"] == FIRST_ID][0])
    rest = trig.iloc[start:]
    out = Path(WORK) / "resume" / "remaining_TriggeredEvents.csv"
    rest.to_csv(out, index=False)
    print(f"{len(trig)} triggered; resuming at row {start + 1} ({FIRST_ID}); {len(rest)} events -> {out}", flush=True)
    s = RQ._scan(WORK, "centred", decimate=False, picker=True)
    s.marginal_window = 0.06
    s.threads = 4
    s.plot_event_summary = False
    s.write_cut_waveforms = False
    s.locate(trigger_file=str(out))


if __name__ == "__main__":
    main()
