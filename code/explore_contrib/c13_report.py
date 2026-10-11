"""C13: assembles REPORT_C13.txt from the prose in c13_report_text.py and the tables of out/c13_report_tables.txt."""
import os
import re

import c13_lib as L
import c13_report_text as T

TAB = open(os.path.join(L.HERE, "out", "c13_report_tables.txt")).read().split("\n")
HEAD = re.compile(r"^T[A-Z]\d?\. ")
SUBS = ("  rule ", "  decoy ", "  STD both ALL, rule", "  DENSE both ALL, rule", "  at t", "  ratio ", "  per stratum",
        "  clutter-only", "  noise-only")


def section(name):
    i0 = next(i for i, x in enumerate(TAB) if x.startswith(name + ". "))
    i1 = next((i for i in range(i0 + 1, len(TAB)) if HEAD.match(TAB[i])), len(TAB))
    out = TAB[i0:i1]
    while out and not out[-1].strip():
        out.pop()
    return out


def sub(name, header, keep=None):
    """the sub-table of section `name` that starts with the line `header` (up to the next sub-table header)."""
    sec = section(name)
    i0 = next(i for i, x in enumerate(sec) if x.startswith(header))
    i1 = next((i for i in range(i0 + 1, len(sec)) if sec[i].startswith(SUBS) or not sec[i].strip()), len(sec))
    out = sec[i0:i1]
    if keep is not None:
        out = [x for k, x in enumerate(out) if k < 2 or keep(x)]
    return out


def head(name, n):
    return section(name)[:n]


# captions and headers of the included tables, in Spanish (whole-line replacements by prefix; then word replacements
# on the short header lines)
LINES_ES = [
    ("TA. Sample and benchmark", "TA. Muestra y puntos de operacion del benchmark (oraculo) de MT-FW: recall con <= 0.1 / <= 0.25 eventos falsos"),
    ("    window, and the false discovery proportion", "    por ventana realizados, y proporcion de falsos (FDP) entre los atomos aceptados en esos puntos"),
    ("TC. Windows with events:", "TC. Ventanas con eventos: falsos del blanco frente a senuelos. Atomos por ventana con eventos por encima de t80 | t95"),
    ("    of the target scores of the clutter-only", "    (cuantiles 0.80 y 0.95 de las puntuaciones del blanco en las ventanas CLUT del estrato). false = falsos del"),
    ("    tmiss (within 5 km", "    blanco bajo LOC+T, de los cuales tmiss (a <= 5 km, error de origen > 0.5 s), near (5-10 km), far (> 10 km);"),
    ("    labelled as hits;", "    fLOC = falsos bajo LOC; hit = atomos etiquetados como acierto; columnas D* = atomos senuelo sobre el umbral"),
    ("  at t80", "  en t80"), ("  at t95", "  en t95"),
    ("  ratio false targets (LOC+T) / decoys at t80", "  cociente falsos del blanco (LOC+T) / senuelos en t80 [95%]; un senuelo valido necesita <= 1"),
    ("  ratio false targets under LOC / decoys", "  cociente falsos bajo LOC / senuelos en t80 y t95 [95%], y distancia KS de las puntuaciones de los falsos (LOC+T)"),
    ("  clutter-only windows (CLUT + empty STD)", "  ventanas de solo clutter (CLUT + STD vacias); KS = distancia de Kolmogorov-Smirnov [p de permutacion] del atomo\n  maximo (top) y de todos los atomos (all); r80, r95 = #blancos / #senuelos sobre los cuantiles 0.80 y 0.95 [95%]"),
]
WORDS_ES = [("  rule ", "  regla "), ("  decoy ", "  senuelo "), (", rule ", ", regla "), (", batch ", ", lote "),
            (" windows, q", " ventanas, q"), ("  line  ", "  linea "), ("regime cell", "regimen celda"),
            ("      ratio ", "      cociente "), (" win ", " ven "), ("meth ", "regla")]


def spanish(lines):
    out = []
    for x in lines:
        for pre, rep in LINES_ES:
            if x.startswith(pre):
                x = rep
                break
        else:
            if x.startswith(("  rule ", "  decoy ", "  line", "  STD both ALL, rule", "  DENSE both ALL, rule", "  regime",
                             "      ratio", "    ")) and ("[" not in x or x.startswith("      ratio")):
                for a, b in WORDS_ES:
                    x = x.replace(a, b)
        out.extend(x.split("\n"))
    return out


def build():
    out = []
    for item in T.PARTS:
        if isinstance(item, str):
            out.extend(item.strip("\n").split("\n"))
        else:
            kind = item[0]
            if kind == "section":
                out.extend(spanish(section(item[1])))
            elif kind == "head":
                out.extend(spanish(head(item[1], item[2])))
            elif kind == "td":
                sec = section("TD")
                i0 = next(i for i, x in enumerate(sec) if x.strip() == item[1])
                i1 = next((i for i in range(i0 + 1, len(sec)) if sec[i][:2] == "  " and sec[i][2:3] != " "), len(sec))
                out.extend(spanish(sec[i0:i1]))
            elif kind == "sub":
                out.extend(spanish(sub(item[1], item[2], item[3] if len(item) > 3 else None)))
        out.append("")
    txt = "\n".join(out).rstrip("\n") + "\n"
    bad = sorted({c for c in txt if ord(c) > 127})
    assert not bad, bad
    wide = max(len(x) for x in txt.split("\n"))
    open(os.path.join(L.HERE, "REPORT_C13.txt"), "w").write(txt)
    print("REPORT_C13.txt:", len(txt.split("\n")), "lines; widest line", wide)


if __name__ == "__main__":
    build()
