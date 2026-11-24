#!/usr/bin/env python3
"""Rearrange the FISH panels of Extended Data Fig. 9 so that they read in the order of the text.

The Illustrator file draws one panel c: the chr9 parental FISH (four image columns) on the left
and the chr12 and chr19 FISH (two image columns, one sample per row) on the right. The text
discusses chr12 and chr19 first and the parents second. This script takes the PDF written by
nature_ai_export.py and swaps the two blocks, without redrawing anything:

    c  chr12 (top left)        e  chr9 parents (right, father above mother)
    d  chr19 (bottom left)

Everything above panel c (panels a and b) is copied unchanged. The blocks are found from the
positions of the FISH images, so nothing is hard-coded except the panel letter to split.

    python illustrator/nature_ai_export.py HPRC_cent_analysis.suppl_tpo.YC.ai ED_Fig9_export \
        --max-pt 7 --anchor center --drop-small --arial
    python illustrator/ed_fig9_panels.py ED_Fig9_export.pdf ED_Fig9

writes ED_Fig9.pdf (vector, images untouched) and ED_Fig9.tif (RGB, 300 ppi).
"""

import argparse
import os
import sys
from pathlib import Path

import fitz  # PyMuPDF
from PIL import Image

PPI = 300
_FONTS = os.environ.get("ARIAL_DIR", "/System/Library/Fonts/Supplemental")
ARIAL_BOLD = os.environ.get("ARIAL_BOLD", f"{_FONTS}/Arial Bold.ttf")


def spans(page):
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                if s["text"].strip():
                    yield s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdf", help="ED Fig. 9 as written by nature_ai_export.py")
    ap.add_argument("out", help="output path without extension")
    a = ap.parse_args()

    src = fitz.open(a.pdf)
    page = src[0]
    W, H = page.rect.width, page.rect.height

    letters = {s["text"].strip(): s for s in spans(page)
               if len(s["text"].strip()) == 1 and s["size"] > 7.5 and "Bold" in s["font"]}
    if sorted(letters) != ["a", "b", "c"]:
        sys.exit(f"expected panel letters a, b, c and found {sorted(letters)}: is this the exporter's PDF?")
    c = fitz.Rect(letters["c"]["bbox"])
    size = letters["c"]["size"]
    base = letters["c"]["origin"][1]

    imgs = [fitz.Rect(i["bbox"]) for i in page.get_image_info() if i["bbox"][1] > c.y0]
    if len(imgs) != 12:
        sys.exit(f"expected 12 FISH images below panel letter c and found {len(imgs)}")
    xs = sorted({round(r.x0, 1) for r in imgs})
    gaps = [(xs[i + 1] - xs[i], i) for i in range(len(xs) - 1)]
    cut_i = max(gaps)[1]                       # the widest step between image columns
    left = [r for r in imgs if r.x0 <= xs[cut_i] + 0.1]      # chr9 parents
    right = [r for r in imgs if r.x0 > xs[cut_i] + 0.1]      # chr12 and chr19
    if (len(left), len(right)) != (8, 4):
        sys.exit(f"expected 8 + 4 images and found {len(left)} + {len(right)}")
    l0, l1 = min(r.x0 for r in left), max(r.x1 for r in left)
    r0, r1 = min(r.x0 for r in right), max(r.x1 for r in right)
    mid = (l1 + r0) / 2
    ysplit = c.y0 - 3                          # just above the letter and the sample titles
    rows = sorted({round(r.y0, 1) for r in right})
    if len(rows) != 2:
        sys.exit(f"expected two rows of images and found {len(rows)}")

    out = fitz.open()
    pg = out.new_page(width=W, height=H)
    pg.show_pdf_page(fitz.Rect(0, 0, W, ysplit), src, 0, clip=fitz.Rect(0, 0, W, ysplit))
    dx_right = l0 - r0                         # chr12 / chr19 block to the left edge
    dx_left = r1 - l1                          # parents block to the right edge
    # The parents block is copied in two pieces that leave out the old letter c (the strip to its
    # right, and the strip below it), so the letter is not in the output at all.
    clips = [fitz.Rect(mid, ysplit, W, H),                 # chr12 and chr19
             fitz.Rect(c.x1 + 1, ysplit, mid, H),          # parents, right of the old letter
             fitz.Rect(0, c.y1 + 1, c.x1 + 1, H)]          # parents, below the old letter
    for clip, dx in zip(clips, (dx_right, dx_left, dx_left)):
        pg.show_pdf_page(clip + (dx, 0, dx, 0), src, 0, clip=clip)

    x = letters["c"]["origin"][0]
    for ch, px, py in (("c", x, base), ("d", x, base + rows[1] - rows[0]), ("e", x + dx_left - 1, base)):
        pg.insert_text(fitz.Point(px, py), ch, fontname="ArialBold", fontfile=ARIAL_BOLD, fontsize=size)

    o = Path(a.out)
    out.save(str(o.with_suffix(".pdf")), garbage=3, deflate=True)
    pix = pg.get_pixmap(matrix=fitz.Matrix(PPI / 72, PPI / 72), colorspace=fitz.csRGB, alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        str(o.with_suffix(".tif")), compression="tiff_lzw", dpi=(PPI, PPI))

    final = fitz.open(str(o.with_suffix(".pdf")))[0]
    got = sorted(s["text"].strip() for s in spans(final)
                 if len(s["text"].strip()) == 1 and s["size"] > 7.5 and "Bold" in s["font"])
    sizes = [s["size"] for s in spans(final)]
    print(f"{o.name}: {W / 72 * 25.4:.1f} x {H / 72 * 25.4:.1f} mm, {pix.width} x {pix.height} px; "
          f"panel letters {''.join(got)}; text {min(sizes):.2f}-{max(sizes):.2f} pt")
    if got != ["a", "b", "c", "d", "e"]:
        sys.exit(f"panel letters are {got}, expected a to e")


if __name__ == "__main__":
    main()
