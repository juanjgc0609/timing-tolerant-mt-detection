import hashlib, os, subprocess, re
HERE = "/home/claude/work"
reg = open(os.path.join(HERE, "PREREGISTRO_1j_1k.txt"), "rb").read()
lines = open(os.path.join(HERE, "PREREGISTRO_hashes.txt")).read().split("\n")
entries = []
for i, ln in enumerate(lines):
    if ln.strip().endswith("PREREGISTRO_1j_1k.txt"):
        entries.append((ln.split()[0], lines[i + 1].strip() if i + 1 < len(lines) else ""))
h = hashlib.sha256()
pref = {}
want = {e[0] for e in entries}
for off in range(len(reg) + 1):
    if off:
        h.update(reg[off - 1:off])
    d = h.hexdigest()
    if d in want:
        pref[d] = off
print(f"registry size {len(reg)} bytes; {len(entries)} hash entries")
prev = -1
ok_chain = True
for d, ts in entries:
    off = pref.get(d)
    tail = reg[:off][-160:].decode("ascii", "replace").replace("\n", " | ") if off is not None else ""
    mono = off is not None and off > prev
    ok_chain &= mono
    print(f"{d[:16]}... {ts}: prefix length {off} ({'monotone' if mono else 'NOT FOUND/NOT MONOTONE'})")
    print(f"     ...ends with: {tail[-150:]}")
    prev = off if off is not None else prev
print("chain:", "OK" if ok_chain else "FAIL")
# prefix at calibration note hash contains note, not result
txt = reg.decode("latin-1")
cal = [d for d, ts in entries if "15:38:20" in ts]
if cal:
    p = reg[:pref[cal[0]]].decode("latin-1")
    print("calibration-note prefix contains 'Nota de calibracion de 2E':", "Nota de calibracion de 2E" in p,
          "; contains 'Resultado de 2E':", "Resultado de 2E" in p)
    idx = txt.find("Tarea 2E (")
    p0 = [d for d, ts in entries if "14:44:14" in ts][0]
    print("2E-section prefix (14:44:14) contains 'Tarea 2E (':", "Tarea 2E (" in reg[:pref[p0]].decode("latin-1"),
          "; contains calibration note:", "Nota de calibracion de 2E" in reg[:pref[p0]].decode("latin-1"))
    print("bytes after calibration-note prefix and before 2E result prefix differ only by appended text:",
          txt[pref[cal[0]]:pref[entries[-1][0]]][:120].replace("\n", " | "))
