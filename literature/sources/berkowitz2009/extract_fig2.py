"""Recover Figure 2 marker centers from the supplied vector PDF (pdfplumber)."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
from statistics import mean

import pdfplumber


HERE = Path(__file__).resolve().parent
DEFAULT_PDF = HERE.parents[1] / "inbox" / "Berkowitz2009Laboratory.pdf"
# Observed fill paths. Each marker also has a stroke path. Reading only its fill
# recovers both the white (open) and shaded (filled) series when they overlap.
# (object, width in points, height, number of drawing commands, fill gray)
SERIES = {
    ("rect", 3.5, 3.5, 5, 1.0): ("FC_square", 0.3, 49),
    ("rect", 3.5, 3.5, 5, 0.0): ("CF_square", 0.3, 63),
    ("curve", 3.9, 3.9, 5, 1.0): ("FC_circle", 0.3, 54),
    ("curve", 3.9, 3.9, 5, 0.492): ("CF_circle", 0.3, 63),
    ("curve", 4.0, 4.7, 4, 1.0): ("FC_triangle", 0.4, 81),
    ("curve", 4.0, 4.7, 4, 0.241): ("CF_triangle", 0.4, 92),
    ("curve", 3.8, 4.3, 7, 1.0): ("FC_hexagon", 0.4, 82),
    ("curve", 5.0, 5.0, 5, 0.747): ("CF_diamond", 0.4, 90),
    ("curve", 4.2, 4.0, 6, 1.0): ("FC_pentagon", 1.0, 64),
    ("curve", 4.2, 4.0, 6, 0.0): ("CF_pentagon", 1.0, 62),
}


def marker_center(obj):
    """Circle and rectangle centers, or the vertex mean of a regular polygon."""
    path = obj["path"]
    if obj["object_type"] == "rect" or any(op[0] == "c" for op in path):
        return (obj["x0"] + obj["x1"]) / 2, (obj["top"] + obj["bottom"]) / 2
    vertices = [op[-1] for op in path if op[0] in ("m", "l")]
    if vertices[-1] != vertices[0]:
        raise ValueError("Expected a closed regular marker polygon")
    return tuple(mean(axis) for axis in zip(*vertices[:-1]))


def extract(pdf):
    expected_hash = json.loads((HERE / "metadata.json").read_text())["pdf_sha256"]
    if hashlib.sha256(pdf.read_bytes()).hexdigest() != expected_hash:
        raise ValueError("The source differs from the reviewed PDF")
    groups = defaultdict(list)
    with pdfplumber.open(pdf) as document:
        page = document.pages[2]
        time_ticks = sorted(line["x0"] for line in page.lines
                            if 3.67 < line["height"] < 3.70 and abs(line["top"] - 241.316) < 0.01)
        ec_ticks = sorted((line["top"] for line in page.lines
                           if 3.67 < line["width"] < 3.70 and abs(line["x0"] - 198.465) < 0.01),
                          reverse=True)
        if (len(time_ticks), len(ec_ticks)) != (13, 11):
            raise ValueError("Tick count differs from the printed axes")
        for obj in [*page.rects, *page.curves]:
            if not obj["fill"]:
                continue
            signature = (obj["object_type"], round(obj["width"], 1), round(obj["height"], 1),
                         len(obj["path"]), obj["non_stroking_color"])
            if signature in SERIES:
                series, _, _ = SERIES[signature]
                groups[series].append(marker_center(obj))
    for name, _, count in SERIES.values():
        if len(groups[name]) != count:
            raise ValueError(f"Expected {count} plotted markers for {name}, found {len(groups[name])}")
    return time_ticks, ec_ticks, groups


def write_markers(output, time_ticks, ec_ticks, groups):
    lines = [
        '# Figure 2, printed p. 3. PDF points with origin at the top left of the page.',
        '# Points are [x, y]. Same-time markers and the separately plotted symbols are retained.',
        'time_tick_min = ' + str(list(range(50, 651, 50))),
        'time_tick_x_pt = ' + str([round(x, 3) for x in time_ticks]),
        'ec_tick_uS = ' + str(list(range(0, 101, 10))),
        'ec_tick_y_pt = ' + str([round(y, 3) for y in ec_ticks]),
    ]
    for name, flow, _ in sorted(SERIES.values(), key=lambda row: (row[1], row[0])):
        group = sorted(groups[name])
        lines.extend(['', '[[series]]', f'id = "{name}"', f'flow_ml_min = {flow}',
                      f'direction = "{name[:2]}"', 'points_pt = ['])
        for first in range(0, len(group), 5):
            lines.append('  ' + ' '.join(f'[{x:.4f}, {y:.4f}],' for x, y in group[first:first + 5]))
        lines.append(']')
    output.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def draw_overlay(pdf, output, groups):
    """Visual QA on the source. The overlay is kept outside the repository."""
    import pypdfium2
    from PIL import ImageDraw

    scale = 4
    with pypdfium2.PdfDocument(pdf) as document:
        page = document[2].render(scale=scale).to_pil().convert("RGB")
    pen = ImageDraw.Draw(page)
    for name, group in groups.items():
        color = "#e52521" if name.startswith("FC") else "#1574c5"
        for x, y in group:
            x, y = scale * x, scale * y
            pen.ellipse((x - 1.6, y - 1.6, x + 1.6, y + 1.6), fill=color)
    page.crop((165 * scale, 46 * scale, 460 * scale, 267 * scale)).save(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overlay", type=Path)
    args = parser.parse_args()
    ticks_x, ticks_y, markers = extract(args.pdf)
    write_markers(args.output, ticks_x, ticks_y, markers)
    if args.overlay:
        draw_overlay(args.pdf, args.overlay, markers)
    print(f"Extracted {sum(map(len, markers.values()))} markers from {len(markers)} symbol series.")
