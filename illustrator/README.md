# Illustrator-composed figures

Figs. 1 and 5, Extended Data Figs. 1-3 and 9, and Supplementary Fig. 1 were
composed in Adobe Illustrator. The `.ai` files are kept in the KaryoScope-paper
repository, not here. `nature_ai_export.py` (from KaryoScope-paper PR #42,
branch `r1i/fig1-panel-order` at `f8f87db`, `scripts/`; adds `--relabel`) makes the submitted PDF and TIFF from each `.ai` file
(PDF-compatible) without modifying it. It:

1. crops the page to the artwork,
2. scales it to 183 mm wide (or 247 mm tall, if that binds first),
3. drops hidden layers,
4. replaces each panel letter with an 8 pt bold lowercase Arial letter at the
   same baseline position, redacting only the letter,
5. optionally shrinks text that would print above `--max-pt` at final size,
6. writes a vector PDF and an RGB TIFF at 300 ppi, and reports text sizes.

```bash
python3 nature_ai_export.py FIGURE.ai OUT_PREFIX [flags]
```

| figure | flags |
|---|---|
| Fig. 1 | `--relabel C=E,E=C` (the chromosome tree becomes c and the subtelomeric tree e, so panels are cited in order; the `.ai` keeps its letters) |
| Fig. 5, ED Fig. 1, ED Fig. 3, Supp. Fig. 1 | none |
| ED Fig. 2 | `--max-pt 7` |
| ED Fig. 9 | `--max-pt 7 --anchor center --drop-small` |

Requires PyMuPDF and Pillow. The replacement letters use Arial Bold from the
macOS location by default; set `ARIAL_BOLD` to the path of `Arial Bold.ttf`
elsewhere (the only change from the paper-repository copy).
