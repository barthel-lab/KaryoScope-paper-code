#!/usr/bin/env python3
"""Export an Illustrator (.ai, PDF-compatible) figure to the Nature figure specification.

For a figure that exists only as an Illustrator file, this script:
  1. crops the page to the artwork,
  2. scales it to the Nature width (183 mm; or to 247 mm tall if that binds first),
  3. replaces each panel letter (a standalone bold capital A-H) with an 8 pt bold lowercase
     Arial letter at the same baseline position (Nature: 8 pt bold upright lowercase),
  4. optionally shrinks text that would print above --max-pt at final size,
  5. writes a vector PDF and an RGB TIFF at 300 ppi, and
  6. reports the text sizes at final size, so too-small text can be found.

The .ai file itself is never modified.

Usage: nature_ai_export.py FIGURE.ai OUT_PREFIX [--max-pt 7] [--min-pt 5] [--width-mm 180]
       [--height-mm 170] [--letter-min-size 9] [--relabel C=E,E=C]
"""

import argparse
import os
import re
import sys
from pathlib import Path

import fitz  # PyMuPDF
from PIL import Image

MM = 72 / 25.4
# Arial for the text this script sets; macOS locations by default, each with an environment override.
_FONTS = os.environ.get("ARIAL_DIR", "/System/Library/Fonts/Supplemental")
ARIAL_BOLD = os.environ.get("ARIAL_BOLD", f"{_FONTS}/Arial Bold.ttf")
ARIAL = os.environ.get("ARIAL", f"{_FONTS}/Arial.ttf")
ARIAL_ITALIC = os.environ.get("ARIAL_ITALIC", f"{_FONTS}/Arial Italic.ttf")
ARIAL_BOLD_ITALIC = os.environ.get("ARIAL_BOLD_ITALIC", f"{_FONTS}/Arial Bold Italic.ttf")
WIDTH_MM, HEIGHT_MM, PPI = 183, 247, 300


def content_rect(page):
    """Union of everything drawn on the page except a full-page background."""
    r = fitz.Rect()
    for kind, bbox in page.get_bboxlog():
        b = fitz.Rect(bbox)
        if b.is_empty or b.width >= page.rect.width * 0.99 and b.height >= page.rect.height * 0.99:
            continue
        r |= b
    for s in spans(page):  # glyph boxes can reach past the drawn shapes (letters, top labels)
        r |= fitz.Rect(s["bbox"])
    return (r + (-2, -2, 2, 2)) & page.rect  # 2 pt margin


def strip_hidden_layers(doc, page):
    """Delete the content of layers that are switched off in Illustrator.

    Rendering skips a hidden layer, but show_pdf_page() copies it without its visibility setting, so
    without this the vector PDF draws hidden layers too.
    """
    off = {x for x, ocg in doc.get_ocgs().items() if not ocg["on"]}
    kind, props = doc.xref_get_key(page.xref, "Resources/Properties")
    if not off or kind != "dict":
        return []
    hidden = {name for name, x in re.findall(r"/(\S+?)\s+(\d+)\s+0\s+R", props) if int(x) in off}
    stream = b"".join(doc.xref_stream(x) for x in page.get_contents())
    out, depth, skip, pos = [], 0, None, 0
    for m in re.finditer(rb"/OC\s*/(\S+)\s+BDC|\bBDC\b|\bBMC\b|\bEMC\b", stream):
        if m.group(0).endswith((b"BDC", b"BMC")):
            if skip is None and m.group(1) and m.group(1).decode() in hidden:
                out.append(stream[pos:m.start()])
                skip = depth
            depth += 1
        else:
            depth -= 1
            if skip is not None and depth == skip:
                pos, skip = m.end(), None
    out.append(stream[pos:])
    page.set_contents(page.get_contents()[0])
    doc.update_stream(page.get_contents()[0], b"".join(out))
    return sorted(doc.get_ocgs()[x]["name"] for x in off)


def spans(page):
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line["spans"]:
                if s["text"].strip():
                    yield dict(s, dir=tuple(line["dir"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ai")
    ap.add_argument("out")
    ap.add_argument("--max-pt", type=float, default=None,
                    help="shrink text larger than this at final size (e.g. 7)")
    ap.add_argument("--anchor", choices=["left", "center"], default="left",
                    help="keep each shrunk label's left edge or its centre (centred titles)")
    ap.add_argument("--drop-small", action="store_true",
                    help="remove text that would print under 5 pt (only where it repeats other content)")
    ap.add_argument("--min-pt", type=float, default=None,
                    help="enlarge text that would print under this size (e.g. 5); a label that would "
                         "then overlap its neighbour in a stack of tick labels is removed")
    ap.add_argument("--width-mm", type=float, default=WIDTH_MM,
                    help="final width (183 for Extended Data; main figures are at most 180)")
    ap.add_argument("--height-mm", type=float, default=HEIGHT_MM,
                    help="final height limit (247 for Extended Data; main figures are at most 170)")
    ap.add_argument("--letter-min-size", type=float, default=9,
                    help="nominal size above which a bold single capital counts as a panel letter")
    ap.add_argument("--relabel", default="",
                    help="re-letter panels, e.g. C=E,E=C (the .ai keeps its original letters)")
    ap.add_argument("--arial", action="store_true",
                    help="set all remaining text that is not Arial again in embedded Arial "
                         "(same origin, size, colour and direction)")
    a = ap.parse_args()
    relabel = dict(pair.split("=") for pair in a.relabel.split(",")) if a.relabel else {}

    doc = fitz.open(a.ai)
    page = doc[0]
    hidden = strip_hidden_layers(doc, page)
    clip = content_rect(page)
    scale = min(a.width_mm * MM / clip.width, a.height_mm * MM / clip.height)

    # Panel letters: bold, single capital A-H, large.
    letters = [s for s in spans(page)
               if re.fullmatch(r"[A-H]", s["text"].strip()) and "Bold" in s["font"]
               and s["size"] >= a.letter_min_size]
    for s in letters:
        # Only the capital's ink (baseline up to cap height): the full glyph box reaches into
        # neighbouring labels, e.g. a vertical scale-bar label just below the letter.
        x0, _, x1, _ = s["bbox"]
        y = s["origin"][1]
        page.add_redact_annot(fitz.Rect(x0, y - 0.75 * s["size"], x1, y + 0.02 * s["size"]))
    # Oversized text, if asked: remove and re-insert smaller at the same baseline (left-aligned).
    shrink = []
    if a.max_pt:
        for s in spans(page):
            if s in letters or s["size"] * scale <= a.max_pt:
                continue
            shrink.append(s)
            page.add_redact_annot(fitz.Rect(s["bbox"]))
    dropped = []
    if a.drop_small:
        for s in spans(page):
            if s not in letters and s["size"] * scale < 5:
                dropped.append(s)
                # Only the middle of the line: glyph boxes of stacked labels overlap vertically.
                x0, y0, x1, y1 = s["bbox"]
                h = y1 - y0
                page.add_redact_annot(fitz.Rect(x0, y0 + 0.35 * h, x1, y1 - 0.35 * h))
    # Undersized text, if asked: remove and re-insert at --min-pt, keeping the label's alignment.
    grow, thinned = [], []
    if a.min_pt:
        small = [s for s in spans(page) if s not in letters and s not in dropped
                 and s["size"] * scale < a.min_pt]
        size = a.min_pt / scale
        horiz = lambda s: s["dir"] == (1.0, 0.0)
        for s in small:
            # Alignment: labels that share a left (or right) edge with two others are a
            # left- (right-) aligned column of tick labels; anything else is centred.
            lo, hi = (0, 2) if horiz(s) else (1, 3)
            same = lambda k: sum(abs(o["bbox"][k] - s["bbox"][k]) < 0.3 and o["dir"] == s["dir"]
                                 for o in small)
            s["anchor"] = ("center" if not horiz(s) else
                           "start" if same(lo) >= 3 else "end" if same(hi) >= 3 else "center")
            font = "hebo" if "Bold" in s["font"] else "helv"
            s["newlen"] = fitz.get_text_length(s["text"], fontname=font, fontsize=size)
        # A stack of tick labels that would overlap at the new size: keep every other one.
        keep = []
        for s in sorted(small, key=lambda s: (round(s["bbox"][0]), -s["origin"][1])):
            # (only text well under the limit: a label within 10% of it grows too little to collide)
            if s["size"] * scale < 0.9 * a.min_pt and horiz(s) and s["anchor"] != "center" \
                    and keep and horiz(keep[-1]) \
                    and abs(keep[-1]["bbox"][0] - s["bbox"][0]) < 0.3 \
                    and keep[-1]["origin"][1] - s["origin"][1] < size * 1.2:
                thinned.append(s)
            else:
                keep.append(s)
        grow, grow_obstacles = keep, []
        # New origin of each enlarged label, then two nudges so enlarged labels stay clear:
        # a centred horizontal label moves down if it would touch the label above it, and a
        # vertical label moves right if it would touch horizontal labels on its left.
        for s in grow:
            x, y = s["origin"]
            old = (s["bbox"][2] - s["bbox"][0]) if horiz(s) else (s["bbox"][3] - s["bbox"][1])
            shift = {"start": 0, "center": (s["newlen"] - old) / 2, "end": s["newlen"] - old}[s["anchor"]]
            s["new"] = (x - shift, y) if horiz(s) else (x, y + shift)
        for o in spans(page):  # labels that keep their size are obstacles too
            if not any(o["bbox"] == g["bbox"] and o["text"] == g["text"] for g in small) and horiz(o):
                grow_obstacles.append(dict(o, new=o["origin"], newlen=o["bbox"][2] - o["bbox"][0]))
        for s in grow:
            x, y = s["new"]
            if horiz(s) and s["anchor"] == "center":
                above = [o["new"][1] for o in grow + grow_obstacles if o is not s and horiz(o)
                         and 0 < y - o["new"][1] < size * 1.15
                         and o["new"][0] < x + s["newlen"] and x < o["new"][0] + o["newlen"]]
                if above:
                    s["new"] = (x, max(above) + size * 1.15)
            elif not horiz(s):
                left = [o["new"][0] + o["newlen"] for o in grow if horiz(o)
                        and o["new"][0] < x and o["new"][0] + o["newlen"] > x - 0.75 * size
                        and y - s["newlen"] < o["new"][1] and o["new"][1] - 0.72 * size < y]
                if left:
                    s["new"] = (max(left) + 0.75 * size + 0.1 * size, y)
        for s in small:
            x0, y0, x1, y1 = s["bbox"]
            if horiz(s):
                h = y1 - y0
                page.add_redact_annot(fitz.Rect(x0, y0 + 0.35 * h, x1, y1 - 0.35 * h))
            else:
                w = x1 - x0
                page.add_redact_annot(fitz.Rect(x0 + 0.35 * w, y0, x1 - 0.35 * w, y1))
    # Text that a redaction box touches is removed too (stacked labels overlap); put it back.
    rects = [annot.rect for annot in page.annots(types=[fitz.PDF_ANNOT_REDACT])]
    collateral = [s for s in spans(page)
                  if not any(s["bbox"] == t["bbox"] and s["text"] == t["text"]
                             for t in letters + shrink + dropped + grow + thinned)
                  and any(fitz.Rect(s["bbox"]).intersects(r) for r in rects)]
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                          graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                          text=fitz.PDF_REDACT_TEXT_REMOVE)
    # A redaction box can touch a glyph box without removing the text; re-insert only what is gone.
    left = {(s["text"], tuple(round(v) for v in s["bbox"])) for s in spans(page)}
    collateral = [s for s in collateral
                  if (s["text"], tuple(round(v) for v in s["bbox"])) not in left]
    page.insert_font(fontname="AB", fontfile=ARIAL_BOLD)
    # Text set again by this script: embedded Arial with --min-pt. Without it, the base-14
    # Helvetica names as before, so earlier exports stay reproducible pixel for pixel. Arial has
    # the same advance widths as Helvetica, so the lengths measured above hold for both.
    if a.min_pt:
        page.insert_font(fontname="AR", fontfile=ARIAL)
    def face(s):
        bold = "Bold" in s["font"]
        return ("AB" if bold else "AR") if a.min_pt else ("hebo" if bold else "helv")
    for s in letters:
        letter = s["text"].strip()
        page.insert_text(fitz.Point(s["origin"]), relabel.get(letter, letter).lower(),
                         fontname="AB", fontsize=8 / scale, color=(0, 0, 0))
    for s in collateral:
        x, y = s["origin"]
        if s["dir"] == (0.0, -1.0):  # reading upwards
            page.insert_text(fitz.Point(x, y), s["text"],
                             fontname=face(s),
                             fontsize=s["size"], color=fitz.sRGB_to_pdf(s["color"]), rotate=90)
            continue
        if s["dir"] != (1.0, 0.0):
            sys.exit(f"rotated text {s['text']!r} touches a redaction; re-inserting it is not supported")
        # If the line stacked under this one was dropped, recentre this one in the freed space.
        below = [d for d in dropped if d["origin"][1] > y and
                 min(d["bbox"][2], s["bbox"][2]) > max(d["bbox"][0], s["bbox"][0])]
        if below:
            y += (min(d["origin"][1] for d in below) - y) / 2
        page.insert_text(fitz.Point(x, y), s["text"],
                         fontname=face(s),
                         fontsize=s["size"], color=fitz.sRGB_to_pdf(s["color"]))
    for s in shrink:
        font = "hebo" if "Bold" in s["font"] else "helv"
        size = a.max_pt / scale
        x, y = s["origin"]
        if a.anchor == "center":
            centre = (s["bbox"][0] + s["bbox"][2]) / 2
            x = centre - fitz.get_text_length(s["text"], fontname=font, fontsize=size) / 2
        page.insert_text(fitz.Point(x, y), s["text"], fontname=face(s), fontsize=size,
                         color=fitz.sRGB_to_pdf(s["color"]))

    for s in grow:
        font = "hebo" if "Bold" in s["font"] else "helv"
        size = a.min_pt / scale
        x, y = s["new"]
        if s["dir"] == (1.0, 0.0):
            page.insert_text(fitz.Point(x, y), s["text"], fontname=face(s), fontsize=size,
                             color=fitz.sRGB_to_pdf(s["color"]))
        elif s["dir"] == (0.0, -1.0):  # reading upwards
            page.insert_text(fitz.Point(x, y), s["text"], fontname=face(s), fontsize=size,
                             color=fitz.sRGB_to_pdf(s["color"]), rotate=90)
        else:
            sys.exit(f"text direction {s['dir']} of {s['text']!r} is not supported by --min-pt")
    # Any text still in another typeface (Helvetica from Illustrator): the same string in Arial at
    # the same origin. Arial has Helvetica's advance widths, so labels keep their extent.
    converted = []
    if a.arial:
        # All text comes off and goes back, including text that is Arial already: a redaction box
        # also removes single characters of neighbouring labels, which would then be left in part.
        before = list(spans(page))
        converted = [s for s in before if "Arial" not in s["font"]]
        for s in before:
            page.add_redact_annot(fitz.Rect(s["bbox"]))
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                              graphics=fitz.PDF_REDACT_LINE_ART_NONE,
                              text=fitz.PDF_REDACT_TEXT_REMOVE)
        if list(spans(page)):
            sys.exit(f"--arial: text left after removal: {[s['text'] for s in spans(page)][:5]}")
        faces = {"AR": ARIAL, "AB": ARIAL_BOLD, "AI": ARIAL_ITALIC, "ABI": ARIAL_BOLD_ITALIC}
        for name, path in faces.items():
            page.insert_font(fontname=name, fontfile=path)
        for s in before:
            bold = "Bold" in s["font"]
            ital = "Italic" in s["font"] or "Oblique" in s["font"]
            name = {(0, 0): "AR", (1, 0): "AB", (0, 1): "AI", (1, 1): "ABI"}[(bold, ital)]
            turn = {(1.0, 0.0): 0, (0.0, -1.0): 90, (-1.0, 0.0): 180, (0.0, 1.0): 270}.get(s["dir"])
            if turn is None:
                sys.exit(f"--arial: text direction {s['dir']} of {s['text']!r} is not supported")
            page.insert_text(fitz.Point(s["origin"]), s["text"], fontname=name, fontsize=s["size"],
                             color=fitz.sRGB_to_pdf(s["color"]), rotate=turn)
    if grow:  # enlarged labels can reach past the first crop
        clip = content_rect(page)
        scale = min(a.width_mm * MM / clip.width, a.height_mm * MM / clip.height)

    # RGB TIFF at 300 ppi, alpha flattened onto white. Render before set_cropbox(): after it, page
    # coordinates are relative to the crop box, and `clip` would cut off the top and left edges.
    out = Path(a.out)
    pix = page.get_pixmap(matrix=fitz.Matrix(scale * PPI / 72, scale * PPI / 72), clip=clip,
                          colorspace=fitz.csRGB, alpha=False)
    img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
    img.save(str(out.with_suffix(".tif")), compression="tiff_lzw", dpi=(PPI, PPI))

    # Vector PDF cropped to the artwork, then scaled to final size.
    vec = fitz.open()
    vp = vec.new_page(width=clip.width * scale, height=clip.height * scale)
    vp.show_pdf_page(vp.rect, doc, 0, clip=clip)
    vec.save(str(out.with_suffix(".pdf")), garbage=3, deflate=True)

    sizes = [round(s["size"] * scale, 2) for s in spans(page)
             if not re.fullmatch(r"[a-h]", s["text"].strip())]
    small = sum(x < 5 for x in sizes)
    large = sum(x > 7 for x in sizes)
    print(f"{Path(a.ai).name}: {clip.width * scale / MM:.0f} x {clip.height * scale / MM:.0f} mm, "
          f"{img.size[0]} x {img.size[1]} px; panel letters replaced: "
          f"{''.join(sorted(s['text'].strip() for s in letters))}"
          + (f" (relabelled {a.relabel})" if relabel else "") + f"; text shrunk: {len(shrink)}"
          + (f"; dropped under 5 pt: {len(dropped)}" if a.drop_small else "")
          + (f"; enlarged to {a.min_pt} pt: {len(grow)}, tick labels removed: "
             f"{[s['text'] for s in thinned]}" if a.min_pt else "")
          + (f"; re-inserted (collateral): {len(collateral)}" if collateral else "")
          + (f"; set again in Arial: {len(converted)}" if a.arial else "")
          + (f"; hidden layers removed: {', '.join(hidden)}" if hidden else ""))
    if sizes:
        print(f"  text at final size: {min(sizes)}-{max(sizes)} pt; "
              f"{small}/{len(sizes)} spans < 5 pt, {large}/{len(sizes)} > 7 pt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
