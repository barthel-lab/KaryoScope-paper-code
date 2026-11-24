#!/usr/bin/env python3
"""Write an RGB copy of an Illustrator file saved in CMYK, with the source RGB colours restored.

ED Figs. 2 and 3 and Supplementary Fig. 1 were composed in Illustrator documents in CMYK
(U.S. Web Coated SWOP v2). Illustrator converted every placed RGB colour to CMYK, so the exported
figures show shifted colours (#fa0000 prints as #ed1c24, #0081fa as #4179bd) that no longer match
the tracks drawn from the colour files. This script finds, for every CMYK colour operator in the
.ai, the RGB colour it was converted from, and writes a copy in which that operator sets the RGB
colour instead. Nothing else in the file changes. The .ai itself is never modified.

    python ai_cmyk_to_rgb.py supfig2.ai supfig2.rgb.pdf --palettes palettes.json \\
        --panel A=karyo_centromere.region --panel B=tree_region --panel C=tree_repeat \\
        --panel D=tree_gene --panel E=tree_acrocentric

The source colours of each panel (palettes.json: panel name -> {hex: count}) are the colours of
the RGB files the panel was made from: the component SVGs of the karyotype panels and the trees
as plot_featureset_tree_bars.draw_panel draws them. For each operator:
  1. its position on the page is found (every operator gets a unique placeholder colour in a
     scratch copy, which is then looked up among the drawn paths and text),
  2. the panel is the panel letter above and to the left of it,
  3. each colour of that panel's palette is converted to CMYK with the document's own profile
     (relative colorimetric, black-point compensation: Illustrator's default), and the operator's
     CMYK value is matched to the nearest; the match must be close (--max-dist) and clear of the
     second nearest (--min-gap), or the script stops.
An operator that draws nothing visible takes the colour its CMYK value resolves to elsewhere on
the page. A report (one line per distinct CMYK value: panels, RGB, distances) is printed.
"""

import argparse
import io
import json
import re
import sys
from collections import Counter, defaultdict

import pymupdf as fitz
from PIL import Image, ImageCms

NUM = r"-?(?:\d+\.?\d*|\.\d+)"
OP = re.compile(rf"({NUM})\s+({NUM})\s+({NUM})\s+({NUM})\s+(scn|SCN|k|K)(?![A-Za-z])")


def content_streams(doc):
    page = doc[0]
    xs = set(page.get_contents())
    for x in range(1, doc.xref_length()):
        try:
            if doc.xref_is_stream(x) and (x in xs or "/Form" in doc.xref_object(x)):
                yield x
        except Exception:
            pass


def cmyk_profile(doc):
    for x in range(1, doc.xref_length()):
        try:
            o = doc.xref_object(x)
        except Exception:
            continue
        m = re.search(r"/ICCBased\s+(\d+)\s+0\s+R", o)
        if m:
            data = doc.xref_stream(int(m.group(1)))
            if data[16:20] == b"CMYK":
                return ImageCms.ImageCmsProfile(io.BytesIO(data))
    sys.exit("no CMYK ICC profile in the file: nothing to convert")


def rewrite(doc, colour_of):
    """Replace every CMYK operator; colour_of(stream xref, occurrence index, cmyk, op) -> op text."""
    for x in content_streams(doc):
        s = doc.xref_stream(x).decode("latin1")
        k = iter(range(10 ** 9))
        new = OP.sub(lambda m: colour_of(x, next(k), tuple(float(v) for v in m.group(1, 2, 3, 4)), m.group(5)), s)
        if new != s:
            doc.update_stream(x, new.encode("latin1"))


def drop_cmyk_space(doc):
    """Remove the now unused CMYK ICC colour space: its names in /ColorSpace resource dicts and the
    'cs' / 'CS' operators that select it (every colour after them is set with rg / RG)."""
    names = set()
    for x in range(1, doc.xref_length()):
        try:
            o = doc.xref_object(x)
        except Exception:
            continue
        for name, ref in re.findall(r"/(\w+)\s*\[\s*/ICCBased\s+(\d+)\s+0\s+R\s*\]", o):
            if doc.xref_stream(int(ref))[16:20] == b"CMYK":
                names.add(name)
        for name, ref in re.findall(r"/(\w+)\s+(\d+)\s+0\s+R", o):
            try:
                t = doc.xref_object(int(ref))
            except Exception:
                continue
            m = re.fullmatch(r"\s*\[\s*/ICCBased\s+(\d+)\s+0\s+R\s*\]\s*", t)
            if m and doc.xref_stream(int(m.group(1)))[16:20] == b"CMYK":
                names.add(name)
    n = 0
    pat = re.compile(r"/(%s)\s+(cs|CS)(?![A-Za-z])" % "|".join(map(re.escape, names))) if names else None
    for x in content_streams(doc):
        s = doc.xref_stream(x).decode("latin1")
        if pat:
            s2, k = pat.subn("", s)
            if k:
                n += k
                doc.update_stream(x, s2.encode("latin1"))
        if OP.search(s):
            sys.exit("a CMYK colour operator survived")
    for x in range(1, doc.xref_length()):
        try:
            kind, val = doc.xref_get_key(x, "Resources/ColorSpace")
        except Exception:
            continue
        targets = [(x, "Resources/ColorSpace")] if kind == "dict" else []
        try:
            k2, v2 = doc.xref_get_key(x, "ColorSpace")
            if k2 == "dict":
                targets.append((x, "ColorSpace"))
        except Exception:
            pass
        for xx, key in targets:
            for name in names:
                if doc.xref_get_key(xx, f"{key}/{name}")[0] != "null":
                    doc.xref_set_key(xx, f"{key}/{name}", "null")
    return n


def rgb_op(rgb, op):
    r, g, b = (v / 255 for v in rgb)
    return f"{r:.6f} {g:.6f} {b:.6f} " + ("rg" if op in ("scn", "k") else "RG")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ai")
    ap.add_argument("out")
    ap.add_argument("--palettes", required=True)
    ap.add_argument("--panel", action="append", required=True, metavar="LETTER=PALETTE")
    ap.add_argument("--max-dist", type=float, default=0.05, help="largest CMYK distance (0-1, max over inks) accepted")
    ap.add_argument("--min-gap", type=float, default=0.03, help="smallest margin to the second nearest colour")
    a = ap.parse_args()
    palettes = json.load(open(a.palettes))
    panel_pal = {}
    for p in a.panel:
        letter, name = p.split("=")
        panel_pal[letter] = sorted(palettes[name])

    src = fitz.open(a.ai)
    page = src[0]
    letters = [(s["text"].strip(), fitz.Rect(s["bbox"]))
               for b in page.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]
               if re.fullmatch(r"[A-H]", s["text"].strip()) and "Bold" in s["font"] and s["size"] >= 15]
    if sorted(t for t, _ in letters) != sorted(panel_pal):
        sys.exit(f"panel letters in the file {sorted(t for t, _ in letters)} differ from --panel {sorted(panel_pal)}")

    def panel_at(x, y):
        cand = [(round(r.x0 / 10), r.y0, t) for t, r in letters if r.x0 - 5 <= x and r.y0 - 5 <= y]
        return max(cand)[2] if cand else None

    # forward conversion of every palette colour with the document's profile
    prof = cmyk_profile(src)
    srgb = ImageCms.createProfile("sRGB")
    t = ImageCms.buildTransform(srgb, prof, "RGB", "CMYK", renderingIntent=ImageCms.Intent.RELATIVE_COLORIMETRIC,
                                flags=ImageCms.Flags.BLACKPOINTCOMPENSATION)
    allhex = sorted({h for v in panel_pal.values() for h in v})
    im = Image.new("RGB", (len(allhex), 1))
    im.putdata([tuple(int(h[i:i + 2], 16) for i in (1, 3, 5)) for h in allhex])
    fwd = {h: tuple(v / 255 for v in px) for h, px in zip(allhex, ImageCms.applyTransform(im, t).getdata())}

    def resolve(cmyk, letter):
        ds = sorted((max(abs(p - q) for p, q in zip(fwd[h], cmyk)), h) for h in panel_pal[letter])
        (d1, h1), (d2, h2) = ds[0], ds[1] if len(ds) > 1 else (9, None)
        return h1, d1, d2, h2

    # 1. placeholders: occurrence n -> a unique RGB value
    scratch = fitz.open(a.ai)
    occ = {}

    def place(x, k, cmyk, op):
        n = len(occ) + 1
        occ[n] = (x, k, cmyk, op)
        return rgb_op(((n >> 16) & 255, (n >> 8) & 255, n & 255), op)
    rewrite(scratch, place)
    sp = scratch[0]
    where = defaultdict(set)
    for d in sp.get_drawings():
        for key in ("fill", "color"):
            c = d.get(key)
            if c:
                n = sum(round(v * 255) << s for v, s in zip(c, (16, 8, 0)))
                if n in occ:
                    r = d["rect"]
                    where[n].add(panel_at((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
    for b in sp.get_text("dict")["blocks"]:
        for l in b.get("lines", []):
            for s in l["spans"]:
                if s["color"] in occ:
                    r = fitz.Rect(s["bbox"])
                    where[s["color"]].add(panel_at((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))

    # 2. resolve every occurrence by its panel(s); unplaced ones by their CMYK value elsewhere
    by_value = defaultdict(set)
    result = {}
    report = defaultdict(lambda: defaultdict(Counter))
    bad = []
    for n, (x, k, cmyk, op) in occ.items():
        panels = where.get(n, set()) - {None}
        hexes = set()
        for p in panels:
            h, d1, d2, h2 = resolve(cmyk, p)
            if d1 > a.max_dist or d2 - d1 < a.min_gap:
                bad.append((cmyk, p, h, round(d1, 3), h2, round(d2, 3)))
            hexes.add(h)
            report[cmyk][p][(h, round(d1, 3), round(d2, 3))] += 1
        if len(hexes) > 1:
            bad.append((cmyk, sorted(panels), sorted(hexes), "one operator, two panels, two colours"))
        if hexes:
            result[n] = hexes.pop()
            by_value[cmyk].add(result[n])
    for n, (x, k, cmyk, op) in occ.items():
        if n not in result:
            vals = by_value.get(cmyk, set())
            if len(vals) == 1:
                result[n] = next(iter(vals))
                report[cmyk]["(not drawn)"][(result[n], None, None)] += 1
            else:                       # drawn nowhere: resolve against all palettes together
                ds = sorted((max(abs(p - q) for p, q in zip(fwd[h], cmyk)), h) for h in allhex)
                if ds[0][0] > a.max_dist or (len(ds) > 1 and ds[1][0] - ds[0][0] < a.min_gap and ds[1][1] != ds[0][1]):
                    bad.append((cmyk, "not drawn", ds[:2]))
                result[n] = ds[0][1]
                report[cmyk]["(not drawn, all palettes)"][(ds[0][1], round(ds[0][0], 3), None)] += 1
    if bad:
        for b in bad[:20]:
            print("UNRESOLVED", b)
        sys.exit(f"{len(bad)} operators could not be resolved unambiguously")

    # 3. the RGB copy
    out = fitz.open(a.ai)
    by_pos = {(x, k): n for n, (x, k, _, _) in occ.items()}
    rewrite(out, lambda x, k, cmyk, op: rgb_op(tuple(int(result[by_pos[(x, k)]][i:i + 2], 16) for i in (1, 3, 5)), op))
    dropped = drop_cmyk_space(out)
    out.save(a.out, garbage=3, deflate=True)
    print(f"  CMYK colour space removed: {dropped} 'cs'/'CS' operators and its resource entries")
    print(f"{a.ai}: {len(occ)} CMYK operators, {len(report)} distinct values, all resolved -> {a.out}")
    for cmyk, panels in sorted(report.items(), key=lambda kv: -sum(sum(c.values()) for c in kv[1].values())):
        desc = "; ".join(f"{p}: " + ", ".join(f"{h} x{c} (d {d1}, next {d2})" for (h, d1, d2), c in cs.items())
                         for p, cs in sorted(panels.items()))
        print(f"  cmyk {cmyk}: {desc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
