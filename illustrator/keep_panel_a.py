#!/usr/bin/env python3
"""Keep only panel a of an exported figure, without its panel letter.

Supplementary Fig. 1 drew the region hierarchy again as panel b, in a hand-compacted layout. The
same tree is Fig. 1b and Extended Data Fig. 2b, so the figure now shows panel a alone and its legend
points to those. This script takes the PDF written by nature_ai_export.py, cuts the page just above
the panel letter b, and leaves out the letter a. Nothing is redrawn.

    python illustrator/nature_ai_export.py supfig4.rgb.pdf SuppFig1_export --arial --display-names v2 \
        --rename "p arm telomere (satellite featureset)=p arm telomere (region feature set)" \
        --rename "q arm telomere (satellite featureset)=q arm telomere (region feature set)"
    python illustrator/keep_panel_a.py SuppFig1_export.pdf SuppFig1

writes SuppFig1.pdf (vector) and SuppFig1.tif (RGB, 300 ppi).
"""

import argparse
import sys
from pathlib import Path

import fitz  # PyMuPDF
from PIL import Image

PPI = 300


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdf")
    ap.add_argument("out", help="output path without extension")
    ap.add_argument("--margin-mm", type=float, default=0.3, help="white space kept below panel a")
    a = ap.parse_args()

    src = fitz.open(a.pdf)
    page = src[0]
    W = page.rect.width
    letters = {}
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                t = s["text"].strip()
                if len(t) == 1 and s["size"] > 7.5 and "Bold" in s["font"]:
                    letters[t] = fitz.Rect(s["bbox"])
    if sorted(letters) != ["a", "b"]:
        sys.exit(f"expected panel letters a and b, found {sorted(letters)}")
    la, lb = letters["a"], letters["b"]
    # Panel a ends with its row of chromosome names. (The header of the tree in panel b can sit
    # higher than the letter b, so the letter alone does not mark where panel a stops.)
    import re
    names = [fitz.Rect(s["bbox"]).y1 for blk in page.get_text("dict")["blocks"] for ln in blk.get("lines", [])
             for s in ln["spans"] if re.fullmatch(r"chr(\d+|X|Y|X\s*chrY)", s["text"].replace("\xa0", " ").strip())]
    if not names:
        sys.exit("no chromosome names found: is this a karyotype view?")
    cut = max(names) + a.margin_mm * 72 / 25.4
    if cut > lb.y1 + 2:
        sys.exit("the lowest chromosome names are below the letter b: panel a is not on top")

    out = fitz.open()
    pg = out.new_page(width=W, height=cut)
    # two pieces that leave out the letter a: the strip to its right, and the strip below it
    for clip in (fitz.Rect(la.x1 + 1, 0, W, cut), fitz.Rect(0, la.y1 + 1, la.x1 + 1, cut)):
        pg.show_pdf_page(clip, src, 0, clip=clip)
    o = Path(a.out)
    out.save(str(o.with_suffix(".pdf")), garbage=3, deflate=True)
    pix = pg.get_pixmap(matrix=fitz.Matrix(PPI / 72, PPI / 72), colorspace=fitz.csRGB, alpha=False)
    Image.frombytes("RGB", (pix.width, pix.height), pix.samples).save(
        str(o.with_suffix(".tif")), compression="tiff_lzw", dpi=(PPI, PPI))
    final = fitz.open(str(o.with_suffix(".pdf")))[0]
    left = [s["text"].strip() for b in final.get_text("dict")["blocks"] for ln in b.get("lines", [])
            for s in ln["spans"] if len(s["text"].strip()) == 1 and s["size"] > 7.5 and "Bold" in s["font"]]
    if left:
        sys.exit(f"panel letters left in the output: {left}")
    print(f"{o.name}: {W / 72 * 25.4:.1f} x {cut / 72 * 25.4:.1f} mm, {pix.width} x {pix.height} px, no panel letters")


if __name__ == "__main__":
    main()
