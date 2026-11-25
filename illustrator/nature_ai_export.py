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


def rename_spans(spans_, a):
    """Apply --rename and --display-names to the text of the spans, in place. -> [(old, new)]"""
    if not (a.rename or a.display_names):
        return []
    explicit = dict(r.split("=", 1) for r in a.rename)
    if a.display_names:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "common"))
        import display_names as dn
        listed = dn.REGION_V1 if a.display_names == "v1" else dn.REGION_V2
    done = []
    for s in spans_:
        old = s["text"]
        key = old.replace("\xa0", " ").strip()
        new = key
        if key in explicit:
            new = explicit[key]
        elif a.display_names:
            if a.no_multigroup and key.endswith(" " + dn.MULTIGROUP):
                new = key[: -len(dn.MULTIGROUP) - 1]
            else:
                for label in (key, key.replace(" ", "_")):
                    base = label[: -len("_multigroup1")] if label.endswith("_multigroup1") else \
                        label[: -len("_specific")] if label.endswith("_specific") else label
                    if base in listed or base in dn.OTHER:
                        new = dn.display_name(label, a.display_names)
                        break
        if new != key:
            s["renamed"] = True
            s["text"] = old.replace(old.strip(), new.replace(" ", "\xa0") if "\xa0" in old else new)
            done.append((key, new))
    return done


def _label(s):
    return s["text"].replace("\xa0", " ").strip()


def _text_end(s, faces):
    bold = "Bold" in s["font"]
    ital = "Italic" in s["font"] or "Oblique" in s["font"]
    face = faces[{(0, 0): "AR", (1, 0): "AB", (0, 1): "AI", (1, 1): "ABI"}[(bold, ital)]]
    return s["origin"][0] + fitz.Font(fontfile=face).text_length(s["text"], s["size"])


def _redraw(page, moves):
    """moves: [(drawing, dx, dy)] of single-rectangle drawings. All are taken off first, then all
    are drawn again, so that one entry can move into the place another one leaves."""
    for d, _, _ in moves:
        page.add_redact_annot(fitz.Rect(d["rect"]) + (-0.3, -0.3, 0.3, 0.3))
    if moves:
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                              graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                              text=fitz.PDF_REDACT_TEXT_NONE)
    for d, dx, dy in moves:
        page.draw_rect(fitz.Rect(d["rect"]) + (dx, dy, dx, dy), color=d.get("color"), fill=d.get("fill"),
                       width=d.get("width") or 0)


def arrange_key(page, spans_, sort_spec, shift_spec, scale):
    """--sort-key and --shift-key: move legend entries (the label and the swatch at its left). The
    text is off the page at this point (it is set again afterwards), so only swatches are redrawn."""
    rects = [d for d in page.get_drawings() if len(d["items"]) == 1 and d["items"][0][0] == "re"]

    def entry(label):
        found = []
        for s in spans_:
            if _label(s) != label or s["dir"] != (1.0, 0.0):
                continue
            x0, y0, x1, y1 = s["bbox"]
            h = y1 - y0
            sw = [d for d in rects if d.get("fill") is not None and d["rect"].x1 <= x0 + 0.5
                  and x0 - d["rect"].x1 < 1.5 * h and d["rect"].y0 < (y0 + y1) / 2 < d["rect"].y1
                  and d["rect"].width < 2.5 * h]
            if len(sw) == 1:
                box = fitz.Rect(sw[0]["rect"]) + (-1.5, -1.5, 1.5, 1.5)
                found.append((s, [d for d in rects if box.contains(d["rect"])]))
        if len(found) != 1:
            sys.exit(f"legend entry {label!r}: found {len(found)} (need exactly one label with a swatch)")
        return found[0]

    delta = {}                                   # id(span) -> [span, parts, dx_text, dy, dx_swatch]
    if sort_spec:
        # The lists are the columns as they should read. The places are a grid: the x of each
        # existing column (its swatch and its text, taken from the column's first entry, so that
        # every entry of a column lines up exactly) and the rows of the longest existing column.
        lists = [[n.strip() for n in col.split(",")] for col in sort_spec.split(";")]
        ents = {n: entry(n) for col in lists for n in col}

        def sw_x(e):
            return min(d["rect"].x0 for d in e[1])
        cols = []
        for e in sorted(ents.values(), key=lambda e: e[0]["origin"][0]):
            if cols and e[0]["origin"][0] - cols[-1][0][0]["origin"][0] < 3:
                cols[-1].append(e)
            else:
                cols.append([e])
        if len(cols) != len(lists):
            sys.exit(f"--sort-key: {len(lists)} lists for {len(cols)} legend columns")
        rows = sorted(e[0]["origin"][1] for e in max(cols, key=len))
        pitch = (rows[-1] - rows[0]) / (len(rows) - 1)
        first = min(cols[0], key=lambda e: e[0]["origin"][1])
        gap = first[0]["origin"][0] - sw_x(first)          # swatch to text, the same in every column
        for want, col in zip(lists, cols):
            ref = min(col, key=lambda e: e[0]["origin"][1])
            tx = ref[0]["origin"][0]
            sx = tx - gap
            for i, n in enumerate(want):
                sp, parts = ents[n]
                y = rows[i] if i < len(rows) else rows[-1] + pitch * (i - len(rows) + 1)
                delta[id(sp)] = [sp, parts, tx - sp["origin"][0], y - sp["origin"][1], sx - sw_x((sp, parts))]
    if shift_spec:
        names, dx_mm = shift_spec.rsplit(":", 1)
        dx = float(dx_mm) * MM / scale
        for n in names.split(","):
            sp, parts = entry(n.strip())
            e = delta.setdefault(id(sp), [sp, parts, 0.0, 0.0, 0.0])
            e[2] += dx
            e[4] += dx
    _redraw(page, [(d, dxs, dy) for _, parts, _, dy, dxs in delta.values() for d in parts])
    for sp, _, dx, dy, _ in delta.values():
        x0, y0, x1, y1 = sp["bbox"]
        sp["origin"] = (sp["origin"][0] + dx, sp["origin"][1] + dy)
        sp["bbox"] = (x0 + dx, y0 + dy, x1 + dx, y1 + dy)


def tidy_leaders(page, spans_, faces):
    """--tidy-leaders: in the hierarchy trees, start every dotted leader the same distance after its
    label (the generator placed them by an estimated text width, and renamed labels differ in length).
    Where a label would reach its bar, the bars and counts of that tree column move right together."""
    drawings = page.get_drawings()
    leaders = [d for d in drawings if d.get("dashes") and not d["dashes"].startswith("[]")
               and len(d["items"]) == 1 and d["items"][0][0] == "l" and d["rect"].height < 0.5
               and d["rect"].width > 1]
    bars = [d for d in drawings if len(d["items"]) == 1 and d["items"][0][0] == "re"
            and d.get("fill") is not None]
    flat = [s for s in spans_ if s["dir"] == (1.0, 0.0)]
    rows = []
    for d in leaders:
        y, x0, x1 = d["rect"].y0, d["rect"].x0, d["rect"].x1
        on_row = [s for s in flat if s["bbox"][1] < y < s["bbox"][3]]
        left = [s for s in on_row if s["origin"][0] < x0 + 1]
        if not left:
            continue
        lab = max(left, key=lambda s: s["origin"][0])
        bar = [b for b in bars if b["rect"].y0 < y < b["rect"].y1 and -1 < b["rect"].x0 - x1 < 6]
        bar = min(bar, key=lambda b: b["rect"].x0) if bar else None
        count = [s for s in on_row if bar and s["origin"][0] >= bar["rect"].x1 - 0.5]
        count = min(count, key=lambda s: s["origin"][0]) if count else None
        gap = 0.4 * lab["size"]
        rows.append({"d": d, "lab": lab, "bar": bar, "count": count, "gap": gap,
                     "start": _text_end(lab, faces) + gap, "end": x1})
    # a column of a tree = rows whose bars start at the same x; it moves as one
    cols = {}
    for r in rows:
        if r["bar"]:
            cols.setdefault(round(r["bar"]["rect"].x0, 1), []).append(r)
    moves = []
    for col in cols.values():
        need = max(r["start"] - r["end"] for r in col)      # only where a label reaches its bar
        if need <= 0:
            continue
        for r in col:
            r["end"] += need
            moves.append((r["bar"], need, 0.0))
            if r["count"]:
                x0, y0, x1, y1 = r["count"]["bbox"]
                r["count"]["origin"] = (r["count"]["origin"][0] + need, r["count"]["origin"][1])
                r["count"]["bbox"] = (x0 + need, y0, x1 + need, y1)
    for r in rows:
        page.add_redact_annot(fitz.Rect(r["d"]["rect"]) + (-0.5, -1.2, 0.5, 1.2))
    _redraw(page, moves)          # also applies the leader redactions (with any others)
    if rows and not moves:
        page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE,
                              graphics=fitz.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                              text=fitz.PDF_REDACT_TEXT_NONE)
    for r in rows:
        if r["end"] - r["start"] > 1:
            d = r["d"]
            y = d["rect"].y0
            page.draw_line(fitz.Point(r["start"], y), fitz.Point(r["end"], y), color=d.get("color"),
                           width=d.get("width") or 0.4, dashes=d["dashes"])
    return len(rows), len(moves)


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
    ap.add_argument("--display-names", choices=["v1", "v2"], default=None,
                    help="print feature labels by their display names (common/display_names.py); the "
                         "value is the database vocabulary of the figure. Needs --arial")
    ap.add_argument("--no-multigroup", action="store_true",
                    help="with --display-names: drop the '(multigroup)' suffix after parent-row counts")
    ap.add_argument("--rename", action="append", default=[], metavar="OLD=NEW",
                    help="print NEW where a text span reads exactly OLD (repeatable; applied before "
                         "--display-names). Needs --arial")
    ap.add_argument("--sort-key", default=None, metavar="LABEL,LABEL,...[;LABEL,...]",
                    help="put the entries of a legend column in this order, top to bottom (one "
                         "comma-separated list per column, columns separated by ';'). Needs --arial")
    ap.add_argument("--tidy-leaders", action="store_true",
                    help="start every dotted leader of the hierarchy trees the same distance after "
                         "its label. Needs --arial")
    ap.add_argument("--shift-key", default=None, metavar="LABEL,LABEL,...:DX_MM",
                    help="move these legend entries (the label as printed, and the swatch at its left) "
                         "DX_MM to the right at final size, to make room for longer names in the "
                         "column before them. Needs --arial")
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
    renamed = []
    tidied = (0, 0)
    if (a.rename or a.display_names) and not a.arial:
        sys.exit("--rename and --display-names need --arial (the text is set again there)")
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
        renamed = rename_spans(before, a)
        faces = {"AR": ARIAL, "AB": ARIAL_BOLD, "AI": ARIAL_ITALIC, "ABI": ARIAL_BOLD_ITALIC}
        # These redraw line art; they run before the fonts are added, because a redaction drops
        # font resources that are not in use yet.
        if a.sort_key or a.shift_key:
            arrange_key(page, before, a.sort_key, a.shift_key, scale)
        if a.tidy_leaders:
            tidied = tidy_leaders(page, before, faces)
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
          + (f"; labels renamed: {len(renamed)}" if a.arial and (a.rename or a.display_names) else "")
          + (f"; leaders set again: {tidied[0]}, bars moved: {tidied[1]}" if a.arial and a.tidy_leaders else "")
          + (f"; hidden layers removed: {', '.join(hidden)}" if hidden else ""))
    if sizes:
        print(f"  text at final size: {min(sizes)}-{max(sizes)} pt; "
              f"{small}/{len(sizes)} spans < 5 pt, {large}/{len(sizes)} > 7 pt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
