# Figures composed in Illustrator

Figs. 1 and 5, Extended Data (ED) Figs. 1-3 and 9, and Supplementary Fig. 1 were
composed in Adobe Illustrator. Fig. 1 is now made by `hg002/scripts/compose_fig1.py`,
which takes only panel a from its `.ai` file (see `hg002/README.md`). The `.ai` files are kept in the KaryoScope-paper
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
| Fig. 1 (first submission, panels a-e) | `--relabel C=E,E=C` (the chromosome tree becomes c and the subtelomeric tree e, so panels are cited in order, while the `.ai` keeps its letters). The revised Fig. 1 is not exported with this script. |
| Fig. 5 | `--min-pt 5 --width-mm 180 --height-mm 170 --arial`, plus the display-name options (below) |
| ED Fig. 1 | `--min-pt 5 --arial` (the scale-bar label "10 Mbp"), plus `--display-names v2 --no-multigroup` |
| ED Figs. 2 and 3 | `--max-pt 7 --arial`, run on the RGB copy (below). In ED Fig. 3 the cap brings the "subtelomeric" tree header from 7.03 to 7.0 pt. |
| Supp. Fig. 1 | `--arial`, run on the RGB copy (below) |
| ED Fig. 9 | `--max-pt 7 --anchor center --drop-small --arial`, then `ed_fig9_panels.py` on the exported PDF (below) |

`--arial` removes all text and sets the same strings again in embedded Arial
(the Illustrator files carry Helvetica labels, and Nature asks for Arial).

The script requires PyMuPDF and Pillow. The Arial files are read from the macOS
font folder by default. On other systems, set `ARIAL_DIR` to the folder holding
`Arial.ttf`, `Arial Bold.ttf`, `Arial Italic.ttf`, and `Arial Bold Italic.ttf`
(or set `ARIAL`, `ARIAL_BOLD`, `ARIAL_ITALIC`, and `ARIAL_BOLD_ITALIC` one by one).
These overrides are the only change from the copy in the paper repository.

## CMYK documents: ED Figs. 2 and 3, Supplementary Fig. 1

`supfig2.ai`, `supfig3.ai` and `supfig4.ai` (ED Figs. 2 and 3, Supp. Fig. 1) are Illustrator documents
in CMYK (U.S. Web Coated SWOP v2). Illustrator converted every placed RGB colour to CMYK, so a
direct export prints shifted colours (`#fa0000` as `#ed1c24`, `#0081fa` as `#4179bd`) that do not
match the tracks of the other figures. Switching the document to RGB in Illustrator does not undo
this; it converts the CMYK values once more. `ai_cmyk_to_rgb.py` writes an RGB copy of the `.ai`
with the source colours restored: each CMYK colour operator is matched, panel by panel, to the RGB
colour of that panel's source (the component SVG of a karyotype panel, the generator's tree for a
tree panel) that converts to it with the document's own profile, and is replaced by it. Nothing else
changes; the script stops if a match is not close and unambiguous. `source_palettes.py` lists the
source colours. Export the RGB copy with `nature_ai_export.py` as in the table.

```bash
export KS_RESOURCES=$PWD/vendored/KaryoScope-archive-resources/resources
C=KaryoScope-paper/figures/fig1/supplement/components/HG002.KS_human_CHM13_v2.presmoothed
python illustrator/source_palettes.py --out palettes.json \
    --svg karyo_centromere.region=$C.centromere.region.KaryoScope.svg \
    --svg karyo_subtelomere.subtelomeric=$C.subtelomere.subtelomeric.KaryoScope.svg \
    --svg karyo_subtelomere.region=$C.subtelomere.region.KaryoScope.svg \
    --tree region --tree repeat --tree gene --tree acrocentric --tree subtelomeric
python illustrator/ai_cmyk_to_rgb.py supfig2.ai supfig2.rgb.pdf --palettes palettes.json \
    --panel A=karyo_centromere.region --panel B=tree_region --panel C=tree_repeat \
    --panel D=tree_gene --panel E=tree_acrocentric
python illustrator/ai_cmyk_to_rgb.py supfig3.ai supfig3.rgb.pdf --palettes palettes.json \
    --panel A=karyo_subtelomere.subtelomeric --panel B=tree_subtelomeric
python illustrator/ai_cmyk_to_rgb.py supfig4.ai supfig4.rgb.pdf --palettes palettes.json \
    --panel A=karyo_subtelomere.region --panel B=tree_region
python illustrator/nature_ai_export.py supfig2.rgb.pdf ED_Fig2 --max-pt 7 --arial
python illustrator/nature_ai_export.py supfig3.rgb.pdf ED_Fig3 --max-pt 7 --arial
python illustrator/nature_ai_export.py supfig4.rgb.pdf SuppFig1 --arial
```

The other `.ai` files are RGB (DeviceRGB or sRGB) and are exported directly.

## ED Fig. 9: FISH panels in the order of the text

The Illustrator file draws one panel c holding the chr9 parental FISH (left) and the chr12 and
chr19 FISH (right). The text discusses chr12 and chr19 first, so `ed_fig9_panels.py` swaps the two
blocks in the exported PDF and letters them c (chr12), d (chr19), and e (chr9 parents). Nothing is
redrawn: panels a and b are copied unchanged and the FISH images and labels are moved as they are.

```bash
python illustrator/nature_ai_export.py HPRC_cent_analysis.suppl_tpo.YC.ai ED_Fig9_export \
    --max-pt 7 --anchor center --drop-small --arial
python illustrator/ed_fig9_panels.py ED_Fig9_export.pdf ED_Fig9
```

## Display names in the Illustrator figures

The figure scripts print feature labels by the names in `common/display_names.py` when given
`--display-names`. The Illustrator files hold the old labels as live text, and `--arial` sets all
text again, so `nature_ai_export.py` can print the same names without the `.ai` being edited:

- `--display-names v1|v2`: each text span that is a feature label is printed by its display name.
  The value is the database vocabulary of the figure (v1 for Fig. 5, v2 for the HG002 figures).
- `--no-multigroup`: drops "(multigroup)" after the k-mer count of parent rows in the trees.
- `--rename OLD=NEW` (repeatable): for text that is not a bare label (callouts, axis titles).
- `--tidy-leaders`: sets every dotted leader of the hierarchy trees again, starting 0.4 of the text
  size after its label. The trees were drawn with an estimated text width, so the gap varied, and
  renamed labels differ in length.
- `--sort-key LABEL,...[;LABEL,...]`: puts the legend entries (label and swatch) in the given order,
  one list per column. Every entry of a column is set at the same x (swatch and text), on the rows
  of the existing legend. For Fig. 5 this is hierarchy order (`display_sort_key`): each parent directly
  above its children, as in the trees.
- `--shift-key LABEL,...:DX_MM`: moves legend entries to the right, where a longer name in the column
  before them needs the room.

```bash
python illustrator/nature_ai_export.py supfig1.ai ED_Fig1 --min-pt 5 --arial \
    --display-names v2 --no-multigroup --tidy-leaders
python illustrator/nature_ai_export.py supfig2.rgb.pdf ED_Fig2 --max-pt 7 --arial \
    --display-names v2 --no-multigroup --tidy-leaders
python illustrator/nature_ai_export.py supfig3.rgb.pdf ED_Fig3 --max-pt 7 --arial \
    --display-names v2 --no-multigroup --tidy-leaders \
    --rename "p arm telomere (subtelomere featureset)=p arm telomere (subtelomeric feature set)" \
    --rename "q arm telomere (subtelomere featureset)=q arm telomere (subtelomeric feature set)"
python illustrator/nature_ai_export.py HPRC_cent_analysis.v2_tpo.YC.ai Fig5 \
    --min-pt 5 --width-mm 180 --height-mm 170 --arial --display-names v1 \
    --rename "chr3: del Hsat1A=chr3: del HSat1A" --rename "chr5: del Hsat3=chr5: del HSat3" \
    --rename "chr12: dup asat1=chr12: dup aSat" --rename "Hsat1A=HSat1A" --rename "Hsat3=HSat3" \
    --sort-key "centromeric,aSat,aSat HOR,active HOR,divergent HOR,inactive HOR,monomeric,bSat,cenSat,ct,gSat;HSat,HSat1,HSat1A,HSat1B,HSat2,HSat3,rDNA,arm,p arm,q arm" \
    --shift-key "HSat,HSat1,HSat1A,HSat1B,HSat2,HSat3,rDNA,arm,p arm,q arm:3"
python illustrator/nature_ai_export.py supfig4.rgb.pdf SuppFig1_export --arial --display-names v2 \
    --rename "p arm telomere (satellite featureset)=p arm telomere (region feature set)" \
    --rename "q arm telomere (satellite featureset)=q arm telomere (region feature set)"
python illustrator/keep_panel_a.py SuppFig1_export.pdf SuppFig1
```

The `.rgb.pdf` inputs are the RGB copies written by `ai_cmyk_to_rgb.py` (above). Without these
options the exports are unchanged.

Supp. Fig. 1 keeps panel a only (`keep_panel_a.py`). Its panel b drew the region hierarchy again in
a hand-compacted layout in which the longer names do not fit. The same tree is Fig. 1b and ED Fig. 2b.
