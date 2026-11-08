"""compare_tables.py - mean change in recall per method between two outputs of analyze_2c.py or analyze_2d.py (E52).
    python posthoc_e52/compare_tables.py FLOOR.txt NEAREST.txt TABLE_NUMBER"""
import re
import sys

import numpy as np


def table(path, num):
    out, reg, cols = {}, None, None
    for line in open(path):
        m = re.match(r"### Tabla (\d+) \((\w+)[,)]", line)
        if m:
            reg = m.group(2) if m.group(1) == num else None
            cols = None
            continue
        if line.startswith("### "):
            reg = None
        if reg and line.startswith("| Celda"):
            cols = [c.strip() for c in line.strip().strip("|").split("|")]
            continue
        if reg and cols and re.match(r"\| .*\(\d{4}\)", line):
            c = [x.strip() for x in line.strip().strip("|").split("|")]
            for name, v in zip(cols, c):
                if re.fullmatch(r"[-+]?\d\.\d+/[-+]?\d\.\d+", v):
                    out[(reg, c[0], name)] = tuple(map(float, v.split("/")))
    return out


def main(a, b, num):
    A, B = table(a, num), table(b, num)
    for name in dict.fromkeys(k[2] for k in A):
        ks = [k for k in A if k[2] == name and k in B]
        d = np.array([[B[k][j] - A[k][j] for j in (0, 1)] for k in ks])
        print(f"{name:16s} cells {len(ks):2d}  mean change @0.1 {d[:, 0].mean():+.3f}  @0.25 {d[:, 1].mean():+.3f}  "
              f"range @0.25 [{d[:, 1].min():+.3f}, {d[:, 1].max():+.3f}]")


if __name__ == "__main__":
    main(*sys.argv[1:4])
