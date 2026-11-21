# Figures composed in Illustrator

Figs. 1 and 5, Extended Data (ED) Figs. 1-3 and 9, and Supplementary Fig. 1 were
composed in Adobe Illustrator. The `.ai` files are kept in the KaryoScope-paper
repository, not here. `nature_ai_export.py` (from `scripts/` in KaryoScope-paper
`main`, including `--relabel` from PR #42)
makes the submitted PDF and TIFF from each `.ai` file (PDF-compatible) without
modifying it. It:

1. crops the page to the artwork,
2. scales it to 183 mm wide (or 247 mm tall, if that binds first); `--width-mm` and
   `--height-mm` set other limits (main figures are at most 180 x 170 mm),
3. drops hidden layers,
4. replaces each panel letter with an 8 pt bold lowercase Arial letter at the
   same baseline position, redacting only the letter,
5. optionally shrinks text that would print above `--max-pt` at final size, and enlarges text
   that would print below `--min-pt` (a tick label that would then touch its neighbour is removed),
6. writes a vector PDF and an RGB TIFF at 300 ppi, and reports text sizes.

```bash
python3 nature_ai_export.py FIGURE.ai OUT_PREFIX [flags]
```

| figure | flags |
|---|---|
| Fig. 1 | `--relabel C=E,E=C` (the chromosome tree becomes c and the subtelomeric tree e, so panels are cited in order, while the `.ai` keeps its letters) |
| Fig. 5 | `--min-pt 5 --width-mm 180 --height-mm 170 --arial` |
| ED Fig. 1 | `--min-pt 5 --arial` (the scale-bar label "10 Mbp") |
| ED Fig. 3, Supp. Fig. 1 | `--arial` |
| ED Fig. 2 | `--max-pt 7 --arial` |
| ED Fig. 9 | `--max-pt 7 --anchor center --drop-small --arial` |

`--arial` removes all text and sets the same strings again in embedded Arial
(the Illustrator files carry Helvetica labels, and Nature asks for Arial).

The script requires PyMuPDF and Pillow. The Arial files are read from the macOS
font folder by default. On other systems, set `ARIAL_DIR` to the folder holding
`Arial.ttf`, `Arial Bold.ttf`, `Arial Italic.ttf`, and `Arial Bold Italic.ttf`
(or set `ARIAL`, `ARIAL_BOLD`, `ARIAL_ITALIC`, and `ARIAL_BOLD_ITALIC` one by one).
These overrides are the only change from the copy in the paper repository.
