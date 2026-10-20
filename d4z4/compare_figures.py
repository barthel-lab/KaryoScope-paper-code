#!/usr/bin/env python3
"""
Compare regenerated figures with the published PNGs.

    python compare_figures.py --published <dir with D4Z4-*.png> --regenerated results/figures/

For each figure writes <name>.contact.png (published | regenerated | overlay) and
<name>.diff.png (absolute difference), and prints the mean absolute difference.
The overlay shows the published figure in red and the regenerated one in cyan: anything
that lines up is grey/black, anything that moved shows as a red/cyan pair.
"""

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageOps

PUBLISHED = {'fig4': 'D4Z4-fig-1-v2.png', 'figS9': 'D4Z4-sup-fig-1-v2.png', 'figS10': 'D4Z4-sup-fig-2.png'}


def load(path, size=None):
    im = Image.open(path).convert('RGB')
    if size and im.size != size:
        im = im.resize(size, Image.LANCZOS)
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--published', required=True)
    ap.add_argument('--regenerated', default='results/figures')
    ap.add_argument('--out', default=None, help='default: the regenerated directory')
    ap.add_argument('--width', type=int, default=1400, help='width of each image in the contact sheet')
    args = ap.parse_args()
    out = Path(args.out or args.regenerated)
    out.mkdir(parents=True, exist_ok=True)
    for name, pub_name in PUBLISHED.items():
        reg_path = Path(args.regenerated) / f'{name}.png'
        if not reg_path.exists():
            continue
        pub = load(Path(args.published) / pub_name)
        reg = load(reg_path, pub.size)
        a = np.asarray(pub, dtype=np.int16)
        b = np.asarray(reg, dtype=np.int16)
        diff = np.abs(a - b).max(axis=2).astype(np.uint8)
        Image.fromarray(255 - diff).save(out / f'{name}.diff.png')
        ga = np.asarray(ImageOps.grayscale(pub), dtype=np.uint8)
        gb = np.asarray(ImageOps.grayscale(reg), dtype=np.uint8)
        overlay = Image.fromarray(np.stack([gb, ga, ga], axis=2))  # published red, regenerated cyan
        h = round(pub.size[1] * args.width / pub.size[0])
        sheet = Image.new('RGB', (3 * args.width + 40, h + 40), 'white')
        draw = ImageDraw.Draw(sheet)
        for k, (label, im) in enumerate([('published', pub), ('regenerated', reg),
                                         ('overlay: published red, regenerated cyan', overlay)]):
            sheet.paste(im.resize((args.width, h), Image.LANCZOS), (k * (args.width + 20), 40))
            draw.text((k * (args.width + 20) + 5, 10), f'{name}: {label}', fill='black')
        sheet.save(out / f'{name}.contact.png')
        print(f'{name}: mean abs diff {np.abs(a - b).mean():.2f}/255; '
              f'{(diff > 64).mean() * 100:.2f}% of pixels differ by > 64')


if __name__ == '__main__':
    main()
