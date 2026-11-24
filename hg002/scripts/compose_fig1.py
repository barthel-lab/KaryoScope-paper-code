#!/usr/bin/env python3
"""Compose Fig. 1 at the Nature main-figure size: panel a from the Illustrator file, panel b
from the karyotype render.

Panel a is the FISH / KaryoScope schematic, cropped from the published Illustrator file
(PDF-compatible, read with PyMuPDF and never modified) and placed as vector artwork, scaled so
that its smallest text prints at --panel-a-min-pt. Its two row labels (FISH, KaryoScope) are set
again at --row-label-pt in Arial Bold, at the same centre, because at that scale they would
print above 7 pt. Panel b is the SVG written by
KaryoScope_assembly_figure.py --chroms-per-row 12 --no-trees --width-mm 180, converted to PDF
with cairosvg and placed at its own size. Panel letters a and b are 8 pt bold lowercase Arial.

    python compose_fig1.py --ai HG002_figure_v2_v1cen_grid.ai --panel-b fig1b.svg --out Fig1

writes Fig1.pdf (vector, fonts embedded) and Fig1.png (RGB, 300 ppi), and prints the page size
and the smallest and largest text size at final size. Arial must be the font fontconfig returns
for "sans-serif" when the SVG is converted (see hg002/README.md); the script stops if the PDF
holds any other font.
"""

import argparse
import re
import sys
from pathlib import Path

import cairosvg
import fitz  # PyMuPDF

MM = 72 / 25.4


def spans(page, clip=None):
    for b in page.get_text("dict", clip=clip)["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                if s["text"].strip():
                    yield dict(s, dir=tuple(line["dir"]))


def letter_span(page, letter):
    """The panel letter of the Illustrator file: a bold single capital, 20 pt."""
    found = [s for s in spans(page) if s["text"].strip() == letter and "Bold" in s["font"]
             and s["size"] >= 15]
    if len(found) != 1:
        sys.exit(f"panel letter {letter}: found {len(found)} candidates in the .ai, expected 1")
    return found[0]


def panel_a_rect(page, letter_b):
    """Everything drawn above panel b's letter (shapes and text), with a 1 pt margin."""
    limit = letter_b["bbox"][1] - 2
    r = fitz.Rect()
    for _, bbox in page.get_bboxlog():
        b = fitz.Rect(bbox)
        if b.is_empty or b.width >= page.rect.width * 0.99 and b.height >= page.rect.height * 0.99:
            continue
        if b.y1 <= limit:
            r |= b
    for s in spans(page):
        if s["bbox"][3] <= limit:
            r |= fitz.Rect(s["bbox"])
    return r + (-1, -1, 1, 1)


def compact_pieces(src, clip, scale, max_w, gap=None):
    """Panel a with short arrows: the strip left of the first column (row labels, letter A), the four
    columns (header box, FISH box, KaryoScope box, each column's text), and between two columns an
    arrow cut from the .ai: its head and the end of its shaft. Columns are `gap` .ai units apart
    (arrow included), or as far apart as fills max_w page units, at most as far as in the .ai.
    Returns ([(source rect, page x offset)], width on the page, gap)."""
    grey = lambda d: d.get("fill") and all(abs(c - 0.749) < 0.01 for c in d["fill"])
    dr = [d for d in src.get_drawings() if clip.contains(d["rect"])]
    heads = sorted({(round(d["rect"].x0, 1), round(d["rect"].x1, 1)) for d in dr
                    if grey(d) and len(d["items"]) == 3 and 18 < d["rect"].width < 30})
    shafts = [d for d in dr if d["type"] == "s" and d["rect"].height < 0.5 and d.get("width", 0) >= 2]
    if len(heads) != 3:
        sys.exit(f"--compact-a: expected 3 arrow heads between the columns, found {len(heads)}")
    # column extents: everything that is not an arrow, split at the arrows
    arrow_x = [(min(d["rect"].x0 for d in shafts if d["rect"].x1 <= hx0 + 1 and d["rect"].x1 >= hx0 - 1), hx1)
               for hx0, hx1 in heads]
    items = [d["rect"] for d in dr if not (grey(d) and len(d["items"]) == 3 and 18 < d["rect"].width < 30)
             and d not in shafts and d["rect"].width < clip.width * 0.9]
    items += [fitz.Rect(s["bbox"]) for s in spans(src, clip)]
    bounds = [clip.x0] + [v for ax in arrow_x for v in ax] + [clip.x1]
    cols = []
    for k in range(4):
        lo, hi = bounds[2 * k], bounds[2 * k + 1]
        inside = [r for r in items if r.x0 >= lo - 0.5 and r.x1 <= hi + 0.5]
        if k == 0:                         # column 1 starts at its first box, the strip is left of it
            first = min(r.x0 for r in inside if r.width > 100)
            cols.append((first - 1.5, max(r.x1 for r in inside) + 1.5))
            strip = (clip.x0, first - 1.5)
        else:
            cols.append((min(r.x0 for r in inside) - 1.5, max(r.x1 for r in inside) + 1.5))
    content = (strip[1] - strip[0]) + sum(c1 - c0 for c0, c1 in cols)
    orig_gap = [cols[k + 1][0] - cols[k][1] for k in range(3)]
    if gap is None:
        gap = (max_w / scale - content) / 3
    gap = min(gap, min(orig_gap))
    head_w = heads[0][1] - heads[0][0]
    if gap < head_w + 12:
        sys.exit(f"--compact-a: {gap * scale / MM:.1f} mm between columns leaves no room for an arrow "
                 f"(needs {(head_w + 12) * scale / MM:.1f} mm)")
    shaft = min(arrow_x[0][1] - arrow_x[0][0] - head_w, gap - head_w - 8.0)   # at most as in the .ai
    margin = (gap - head_w - shaft) / 2
    pieces, x = [], 0.0
    pieces.append((fitz.Rect(strip[0], clip.y0, strip[1], clip.y1), x))
    x += strip[1] - strip[0]
    for k, (c0, c1) in enumerate(cols):
        pieces.append((fitz.Rect(c0, clip.y0, c1, clip.y1), x * scale))
        x += c1 - c0
        if k < 3:
            hx1 = heads[k][1] + 1.0
            alen = gap - 2 * margin       # arrow length kept: the head and the end of the shaft
            pieces.append((fitz.Rect(hx1 - alen, clip.y0, hx1, clip.y1), (x + margin) * scale))
            x += gap
    return pieces, x * scale, gap


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ai", required=True, help="published Illustrator file (HG002_figure_v2_v1cen_grid.ai)")
    ap.add_argument("--panel-b", required=True, help="panel b SVG from KaryoScope_assembly_figure.py")
    ap.add_argument("--out", required=True, help="output prefix (writes .pdf and .png)")
    ap.add_argument("--arial-bold", required=True, help="Arial Bold .ttf (panel letters, row labels)")
    ap.add_argument("--width-mm", type=float, default=180)
    ap.add_argument("--max-height-mm", type=float, default=170)
    ap.add_argument("--panel-a-min-pt", type=float, default=5.05,
                    help="scale panel a so its smallest text prints at this size")
    ap.add_argument("--row-label-pt", type=float, default=7.0, help="FISH / KaryoScope row labels")
    ap.add_argument("--letter-pt", type=float, default=8.0)
    ap.add_argument("--panel-a-align", choices=["center", "left"], default="center",
                    help="horizontal position of panel a on the page (letter a stays at the left edge)")
    ap.add_argument("--gap-mm", type=float, default=1.5, help="space between panel a and letter b")
    ap.add_argument("--compact-a", action="store_true",
                    help="panel a with short arrows: its four columns placed individually, each arrow "
                         "cut from the .ai (head and the end of its shaft)")
    ap.add_argument("--arrow-gap-mm", type=float, default=None,
                    help="with --compact-a: space between two columns, arrow included (default: "
                         "what fills the width left of the trees)")
    ap.add_argument("--tree", action="append", default=[], metavar="SVG",
                    help="a tree (draw_fig1_trees.py) in the top band right of panel a, with its own "
                         "panel letter; repeatable, left to right")
    ap.add_argument("--letter-slot-mm", type=float, default=4.0,
                    help="room left of each tree for its panel letter")
    ap.add_argument("--ppi", type=int, default=300)
    a = ap.parse_args()
    out = Path(a.out)

    # -- panel a from the Illustrator file -------------------------------------------------
    ai = fitz.open(a.ai)
    src = ai[0]
    letter_a, letter_b = letter_span(src, "A"), letter_span(src, "B")
    clip = panel_a_rect(src, letter_b)
    text_a = [s for s in spans(src, clip) if s is not letter_a]
    body = [s for s in text_a if s["size"] < 18]
    scale_a = a.panel_a_min_pt / min(s["size"] for s in body)
    # Remove the letter A and the two rotated row labels; they are set again on the page below.
    rows = [s for s in text_a if s["dir"] == (0.0, -1.0) and s["size"] >= 18]
    if sorted(s["text"].strip() for s in rows) != ["FISH", "KaryoScope"]:
        sys.exit(f"expected the row labels FISH and KaryoScope, found {[s['text'] for s in rows]}")
    for s in [letter_a] + rows:
        src.add_redact_annot(fitz.Rect(s["bbox"]))
    src.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                         text=fitz.PDF_REDACT_TEXT_REMOVE)
    # Delete the old panels b-e from the copy, so the placed page holds only panel a (a clip
    # alone would keep them, invisible, in the PDF).
    src.add_redact_annot(fitz.Rect(src.rect.x0, clip.y1 + 0.5, src.rect.x1, src.rect.y1))
    src.apply_redactions(images=fitz.PDF_REDACT_IMAGE_REMOVE,
                         graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_TOUCHED,
                         text=fitz.PDF_REDACT_TEXT_REMOVE)
    src.clean_contents(sanitize=True)  # rewrites the content stream and drops unused resources
    kept = [s for s in spans(src) if s["bbox"][1] > clip.y1]
    if kept:
        sys.exit(f"text below panel a survived the redaction: {[s['text'] for s in kept][:5]}")

    # -- panel b ------------------------------------------------------------------------
    b_pdf = out.with_name(out.name + "_panel_b.pdf")
    cairosvg.svg2pdf(url=str(a.panel_b), write_to=str(b_pdf))
    pb = fitz.open(b_pdf)
    b_rect = pb[0].rect

    # -- page ---------------------------------------------------------------------------
    W = a.width_mm * MM
    if abs(b_rect.width - W) > 0.5:
        sys.exit(f"panel b is {b_rect.width / MM:.1f} mm wide; render it with --width-mm {a.width_mm}")
    lp = a.letter_pt
    a_w, a_h = clip.width * scale_a, clip.height * scale_a
    a_top = 0.0
    a_left = (W - a_w) / 2 if a.panel_a_align == "center" else 0.0
    pieces = [(clip, 0.0)]                 # (source rect in the .ai, x offset on the page from a_left)
    trees = []
    for t in a.tree:
        tpdf = out.with_name(out.name + f"_tree{len(trees)}.pdf")
        cairosvg.svg2pdf(url=str(t), write_to=str(tpdf))
        trees.append((tpdf, fitz.open(tpdf)))
    tree_w = sum(td[0].rect.width + a.letter_slot_mm * MM for _, td in trees)
    if a.compact_a:
        a_left = 0.0
        pieces, a_w, gap_units = compact_pieces(src, clip, scale_a, W - tree_w - 1.0 * MM,
                                                a.arrow_gap_mm * MM / scale_a if a.arrow_gap_mm else None)
        print(f"  compact panel a: {a_w / MM:.1f} mm wide, columns {gap_units * scale_a / MM:.1f} mm apart "
              f"(arrow included; {168.0 * scale_a / MM:.1f} mm before)")
    band_h = max([a_h] + [td[0].rect.height for _, td in trees])
    letter_b_base = a_top + band_h + a.gap_mm * MM + 0.72 * lp
    b_top = letter_b_base + 0.6 * MM
    H = b_top + b_rect.height
    if H > a.max_height_mm * MM + 1e-6:
        sys.exit(f"page would be {H / MM:.1f} mm tall (limit {a.max_height_mm} mm)")

    doc = fitz.open()
    page = doc.new_page(width=W, height=H)
    for rect, dx in pieces:
        x0 = a_left + dx
        page.show_pdf_page(fitz.Rect(x0, a_top, x0 + rect.width * scale_a, a_top + rect.height * scale_a),
                           ai, 0, clip=rect)
    page.show_pdf_page(fitz.Rect(0, b_top, W, b_top + b_rect.height), pb, 0)
    page.insert_font(fontname="ArialBold", fontfile=a.arial_bold)
    font = fitz.Font(fontfile=a.arial_bold)

    def to_page(x, y):                   # the strip with the row labels and letter A is piece 0
        return a_left + (x - clip.x0) * scale_a, a_top + (y - clip.y0) * scale_a

    # Row labels: same centre, reading upwards, at --row-label-pt
    for s in rows:
        x0, y0, x1, y1 = s["bbox"]
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        n = font.text_length(s["text"].strip(), fontsize=a.row_label_pt)
        px, py = to_page(cx, cy)
        # origin of upward text: the baseline is the left side of the glyphs, the text starts at the bottom
        page.insert_text(fitz.Point(px + 0.36 * a.row_label_pt, py + n / 2), s["text"].strip(),
                         fontname="ArialBold", fontsize=a.row_label_pt, rotate=90)
    # Panel letters at the page's left edge: a at the .ai letter's height, b above panel b's
    # first header
    ax, ay = to_page(*letter_a["origin"])
    ax -= a_left
    page.insert_text(fitz.Point(max(ax, 0.5), max(ay, 0.75 * lp)), "a", fontname="ArialBold", fontsize=lp)
    letters = iter("bcdefg")
    tx = W
    placed = []
    for tpdf, td in reversed(trees):      # flush right, right to left
        tw, th = td[0].rect.width, td[0].rect.height
        tx -= tw
        placed.append((tx, tw, th, td))
        tx -= a.letter_slot_mm * MM
    for tx, tw, th, td in reversed(placed):
        page.show_pdf_page(fitz.Rect(tx, a_top, tx + tw, a_top + th), td, 0)
        page.insert_text(fitz.Point(tx - a.letter_slot_mm * MM + 0.5, max(ay, 0.75 * lp)), next(letters),
                         fontname="ArialBold", fontsize=lp)
        print(f"  tree at x {tx / MM:.1f}-{(tx + tw) / MM:.1f} mm, {tw / MM:.1f} x {th / MM:.1f} mm")
    if placed and a_left + a_w > placed[-1][0] - a.letter_slot_mm * MM + 1e-6:
        sys.exit(f"panel a ({(a_left + a_w) / MM:.1f} mm) runs into the trees "
                 f"({(placed[-1][0] - a.letter_slot_mm * MM) / MM:.1f} mm)")
    page.insert_text(fitz.Point(max(ax, 0.5), letter_b_base), next(letters), fontname="ArialBold", fontsize=lp)
    doc.subset_fonts()
    doc.save(str(out.with_suffix(".pdf")), garbage=3, deflate=True)
    b_pdf.unlink()
    for tpdf, _ in trees:
        tpdf.unlink()

    # -- checks + PNG -------------------------------------------------------------------
    final = fitz.open(str(out.with_suffix(".pdf")))
    fp = final[0]
    fonts = {f[3] for f in fp.get_fonts(full=True)}
    base = {re.sub(r"^[A-Z]{6}\+", "", f) for f in fonts}
    if not all(f.startswith("Arial") for f in base):
        sys.exit(f"non-Arial fonts in the PDF: {sorted(base)}")
    sizes = [(round(s["size"], 2), s["text"]) for s in spans(fp)]
    letters = [s for s in sizes if re.fullmatch(r"[a-g]", s[1].strip()) and s[0] == a.letter_pt]
    other = [s for s in sizes if s not in letters]
    pix = fp.get_pixmap(dpi=a.ppi, alpha=False, colorspace=fitz.csRGB)
    pix.set_dpi(a.ppi, a.ppi)
    pix.save(str(out.with_suffix(".png")))
    print(f"{out.with_suffix('.pdf')}: {W / MM:.1f} x {H / MM:.1f} mm; panel a {a_w / MM:.1f} x "
          f"{a_h / MM:.1f} mm at x {a_left / MM:.1f} mm (scale {scale_a:.4f} of the .ai), panel b {b_rect.width / MM:.1f} x "
          f"{b_rect.height / MM:.1f} mm")
    print(f"  text: {min(s[0] for s in other)}-{max(s[0] for s in other)} pt over {len(other)} spans; "
          f"panel letters {sorted(set(s[0] for s in letters))} pt; fonts {sorted(base)}")
    print(f"  {out.with_suffix('.png')}: {pix.width} x {pix.height} px at {a.ppi} ppi")
    return 0


if __name__ == "__main__":
    sys.exit(main())
