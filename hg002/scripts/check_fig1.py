#!/usr/bin/env python3
"""Check a composed Fig. 1 (compose_fig1.py) against the Nature limits and its sources.

    python check_fig1.py Fig1.pdf [--reference Fig1_previous.pdf] [--ed-tree ED_Fig2.pdf]

1. page 180 mm wide and at most 170 mm tall; all text Arial, embedded; text 5 to 7 pt; panel letters
   8 pt bold lowercase, in reading order (a, b, ... left to right, then top to bottom);
2. with --reference: the karyotype labels (chromosome names, h1/h2, Mb values, the scale-legend
   strings) are the same multisets as in the reference figure;
3. every tree in the figure (found by its bold root label) equals, row by row, the tree
   plot_featureset_tree_bars draws from the v2 hierarchy: label, k-mer count string (with
   "(multigroup)" where drawn) and bar colour; each bar colour equals the colour of that label in the
   colour file the karyotype tracks are drawn with;
4. with --ed-tree: the region tree's labels (in order) and count strings equal those of the given
   Extended Data PDF (ED Fig. 2b), whose colours went through other software and are not compared;
5. no two text spans overlap in the trees or in the karyotype panel.
Exits 1 if a check fails. KS_RESOURCES must point to vendored/KaryoScope-archive-resources/resources.
"""

import argparse
import collections
import re
import sys
from pathlib import Path

import fitz  # PyMuPDF

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "vendored" / "KaryoScope-ISCN"))
import plot_featureset_tree_bars as tp  # noqa: E402

MM = 72 / 25.4
TREES = ("chromosome", "region", "subtelomeric", "repeat", "gene", "acrocentric")
failed = []


def note(ok, msg):
    print(("   " if ok else " X ") + msg)
    if not ok:
        failed.append(msg)


def hexc(c):
    return "#" + "".join(f"{round(v * 255):02x}" for v in c)


def spans(page, clip=None):
    return [dict(s, rect=fitz.Rect(s["bbox"]), t=s["text"].replace("\xa0", " ").replace("\xad", "-").strip())
            for b in page.get_text("dict", clip=clip)["blocks"] for l in b.get("lines", []) for s in l["spans"]
            if s["text"].strip()]


def near(c, ref):
    return all(abs(((c >> k) & 255) - ((ref >> k) & 255)) <= 3 for k in (0, 8, 16))


def tree_rows(page, clip):
    """{label: (count string, bar colour)} of the tree drawn inside clip, in drawing order. A bar must be
    visible: get_drawings also returns paths that a clip hides (compose_fig1.py places panel a in pieces,
    each the whole panel clipped), so the rendered page must show the bar's colour at its centre."""
    sp = spans(page, clip)
    area = clip & page.rect                  # the render covers only the page
    pix = page.get_pixmap(dpi=600, clip=area, alpha=False)
    z = 600 / 72

    def visible(d):
        r = d["rect"]
        x, y = int(((r.x0 + r.x1) / 2 - area.x0) * z), int(((r.y0 + r.y1) / 2 - area.y0) * z)
        if not (0 <= x < pix.width and 0 <= y < pix.height):
            return False
        return max(abs(a - round(b * 255)) for a, b in zip(pix.pixel(x, y), d["fill"])) <= 2
    labels = [s for s in sp if near(s["color"], 0x333333)]
    counts = [s for s in sp if near(s["color"], 0x555555)]
    bars = [d for d in page.get_drawings() if d.get("fill") and d["type"] == "f" and clip.contains(d["rect"])
            and d["rect"].height > 0.5 and d["rect"].width > 0.3
            and hexc(d["fill"]) not in ("#000000", "#333333", "#555555", "#999999", "#ffffff") and visible(d)]
    rows = {}
    for lb in sorted(labels, key=lambda s: s["rect"].y0):
        cy, hh = (lb["rect"].y0 + lb["rect"].y1) / 2, lb["rect"].height
        same = lambda r: abs((r.y0 + r.y1) / 2 - cy) < 0.4 * hh and r.x0 >= lb["rect"].x1 - 0.5
        bar = min((d for d in bars if same(d["rect"]) and d["rect"].height < hh), key=lambda d: d["rect"].x0)
        cnt = min((c for c in counts if same(c["rect"])), key=lambda c: c["rect"].x0)
        rows[lb["t"]] = (re.sub(r"(\d) (\d)", r"\1.\2", cnt["t"]), hexc(bar["fill"]))
    return rows


def tree_clips(page, letters=()):
    """Each tree's area: from its bold root label to the next tree to its right, and down to the next
    tree or panel letter below it in its column (or the page bottom)."""
    roots = {s["t"]: s["rect"] for s in spans(page) if "Bold" in s["font"] and s["t"] in TREES}
    out = {}
    for k, r in roots.items():
        right = min([o.x0 for o in roots.values() if o.x0 > r.x0 + 5] + [page.rect.width])
        below = min([o.y0 for o in roots.values() if o.y0 > r.y1 and abs(o.x0 - r.x0) < 30]
                    + [l["rect"].y0 for l in letters if l["rect"].y0 > r.y1 + 5 and l["rect"].x0 > r.x0 - 30]
                    + [page.rect.height])
        out[k] = fitz.Rect(r.x0 - 3, r.y0 - 14, right - 1, below - 1)
    return out


def generator_rows(fs, trees):
    rows = {}
    for r in tp._flatten(trees[fs])[1:]:
        n = r["node"]
        cnt = tp.format_count(n.kmer_count) if n.kmer_count > 0 else "0"
        suffix = " (multigroup)" if not r["is_leaf"] and n.kmer_count > 0 else ""
        rows[r["name"].replace("_", " ")] = (cnt + suffix, n.color.lower())
    return rows


def label_groups(sp):
    t = [s["t"] for s in sp]
    return {"chromosome names": sorted(x for x in t if re.fullmatch(r"chr(\d+|X|Y)", x)),
            "h1/h2": sorted(x for x in t if x in ("h1", "h2")),
            "Mb values": sorted(x for x in t if re.fullmatch(r"[\d.]+ (Mb|kb)", x)),
            # exact legend words: the colour key's header "Centromere zoom" must not count
            "scale legend": sorted(x for x in t if re.search(r"^(Genome|Centromere|Telomere)$|zoom\)$", x))}


def overlaps(sp):
    bad = []
    for i, a in enumerate(sp):
        for b in sp[i + 1:]:
            r = a["rect"] & b["rect"]
            if not r.is_empty and abs(r) > 0.15 * min(abs(a["rect"]), abs(b["rect"])):
                bad.append((a["t"], b["t"]))
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdf")
    ap.add_argument("--reference", help="a previous Fig. 1 whose karyotype labels must be the same")
    ap.add_argument("--ed-tree", help="ED Fig. 2 PDF: compare the region tree's labels and counts")
    ap.add_argument("--colors", default=str(tp.V2_DB_DIR / "KS_human_CHM13_v2.{fs}.colors.txt"))
    a = ap.parse_args()
    page = fitz.open(a.pdf)[0]
    sp = spans(page)

    # 1. page, text, fonts, letters
    w, h = page.rect.width / MM, page.rect.height / MM
    note(abs(w - 180) < 0.05 and h <= 170 + 1e-6, f"1. page {w:.1f} x {h:.1f} mm")
    letters = [s for s in sp if re.fullmatch(r"[a-h]", s["t"]) and abs(s["size"] - 8) < 0.01 and "Bold" in s["font"]]
    other = [s for s in sp if s not in letters]
    lo, hi = min(s["size"] for s in other), max(s["size"] for s in other)
    note(lo >= 5 - 1e-6 and hi <= 7 + 1e-6, f"   text {lo:.2f}-{hi:.2f} pt over {len(other)} spans")
    order = [s["t"] for s in sorted(letters, key=lambda s: (round(s["rect"].y0 / 30), s["rect"].x0))]
    note(order == list("abcdefgh"[:len(order)]), f"   panel letters {''.join(order)} (8 pt bold), in reading order")
    fonts = {re.sub(r"^[A-Z]{6}\+", "", f[3]): f[1] for f in page.get_fonts(full=True)}
    doc = page.parent
    embedded = all(doc.extract_font(x)[3] for x in {f[0] for f in page.get_fonts(full=True)})
    note(all(f.startswith("Arial") for f in fonts) and embedded, f"   fonts {sorted(fonts)}, all embedded: {embedded}")

    # karyotype panel: from its letter (the first panel letter below every tree's root label, or the
    # last letter if there is no tree) to the page bottom
    clips = tree_clips(page, letters)
    roots_bottom = max([c.y0 + 14 for c in clips.values()] + [0])
    below = sorted((s for s in letters if s["rect"].y0 > roots_bottom + 5), key=lambda s: s["rect"].y0)
    kletter = below[0] if below else max(letters, key=lambda s: s["rect"].y0)
    kclip = fitz.Rect(0, kletter["rect"].y0, page.rect.width, page.rect.height)
    for fs in clips:                       # a tree ends where the karyotype starts
        clips[fs].y1 = min(clips[fs].y1, kclip.y0 - 1)

    # 2. karyotype labels vs the reference
    if a.reference:
        g, r = label_groups(spans(page, kclip)), label_groups(spans(fitz.open(a.reference)[0]))
        for k in r:
            note(g[k] == r[k], f"2. {k}: {len(g[k])}, same multiset as the reference: {g[k] == r[k]}")

    # 3. trees vs the generator and the colour files
    trees = tp.load_all_featuresets_v2()
    for root in trees.values():
        tp._apply_color_overrides(root)
    for fs, clip in clips.items():
        rows, gen = tree_rows(page, clip), generator_rows(fs, trees)
        same = list(rows.items()) == list(gen.items())
        diff = [(k, rows.get(k), gen.get(k)) for k in dict.fromkeys(list(rows) + list(gen)) if rows.get(k) != gen.get(k)]
        note(same, f"3. {fs} tree: {len(rows)} rows; label, count and colour equal to the generator's tree, row by row: "
                   f"{same}" + (f" (differ: {diff[:4]})" if diff else ""))
        cf = {}
        for line in Path(a.colors.format(fs=fs)).read_text().splitlines()[1:]:
            p = line.split()
            if len(p) >= 2:
                cf[p[0].replace("_", " ")] = p[1].lower()
        off = {k: (v[1], cf[k]) for k, v in rows.items() if k in cf and v[1] != cf[k]}
        note(not off, f"   {fs} tree: bar colour equals the track colour (colour file) for the "
                      f"{sum(k in cf for k in rows)} labels in the file: {not off}" + (f" (differ: {off})" if off else ""))
        if a.ed_tree and fs == "region":
            ed = fitz.open(a.ed_tree)[0]
            ed_letters = [s for s in spans(ed) if re.fullmatch(r"[a-h]", s["t"]) and abs(s["size"] - 8) < 0.01]
            ed_clip = tree_clips(ed, ed_letters).get("region")
            er = tree_rows(ed, ed_clip) if ed_clip else {}
            ok = list(er) == list(rows) and all(er[k][0] == rows[k][0] for k in rows)
            note(ok, f"4. region tree: labels (in order) and counts equal to {Path(a.ed_tree).name}: {ok}")
        ov = overlaps(spans(page, clip))
        note(not ov, f"5. {fs} tree: overlapping text spans: {len(ov)}")
    ov = overlaps(spans(page, kclip))
    note(not ov, f"5. karyotype panel: overlapping text spans: {len(ov)}")
    print("ALL CHECKS PASS" if not failed else f"{len(failed)} CHECKS FAILED")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
