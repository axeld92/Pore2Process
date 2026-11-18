"""Extract the public Figure 2 data and sphere coordinates to an ignored directory."""

import argparse
import csv
import json
from math import hypot
from pathlib import Path
from tempfile import TemporaryFile
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET
from zipfile import ZipFile


ARCHIVE = (
    "https://data.4tu.nl/file/4b8d0613-5ee7-4971-b695-835f287336e8/"
    "6d43bf0a-e0f1-460f-b2e0-3d09ef50a3de"
)
ARCHIVE_SIZE = 315_044_865
# Figure 2 entries, the coordinate sheet, and the archive central directory.
RANGES = ((13_469, 24_525), (295_071_667, 295_123_976), (314_913_793, 315_044_864))
FIGURE_ENTRIES = {
    f"SupplementaryMaterial/Figures/Figure2/Figure_2_{panel}.csv": crc
    for panel, crc in zip("abcdef", (0x8468EC4F, 0xBEABB884, 0x5724B88C,
                                      0x1A101D85, 0xAE26E9FB, 0x869DC44B), strict=True)
}
COORDINATE_ENTRY = "SupplementaryMaterial/Reconstructed packing/Particle positions.xlsx"
OUTPUT = Path("results/adsorption/breakthrough_hypotheses/benchmark")


def read_entries(archive, output):
    for entry, expected_crc in {**FIGURE_ENTRIES, COORDINATE_ENTRY: 0x1148FF19}.items():
        if archive.getinfo(entry).CRC != expected_crc:
            raise ValueError(f"Unexpected archive revision for {entry}")
        name = Path(entry).name if entry != COORDINATE_ENTRY else "particle_positions.xlsx"
        (output / name).write_bytes(archive.read(entry))  # Also verifies CRC after decompression.


def download_entries(output):
    with TemporaryFile() as sparse_archive:
        sparse_archive.truncate(ARCHIVE_SIZE)
        for first, last in RANGES:
            request = Request(ARCHIVE, headers={"Range": f"bytes={first}-{last}"})
            with urlopen(request, timeout=30) as response:
                expected_range = f"bytes {first}-{last}/{ARCHIVE_SIZE}"
                if response.status != 206 or response.headers.get("Content-Range") != expected_range:
                    raise ValueError(f"Expected HTTP 206 and Content-Range {expected_range}")
                data = response.read(last - first + 2)
            if len(data) != last - first + 1:
                raise ValueError(f"Incomplete archive range {first}–{last}")
            sparse_archive.seek(first)
            sparse_archive.write(data)
        with ZipFile(sparse_archive) as archive:
            read_entries(archive, output)


def convert_coordinates(output):
    namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with ZipFile(output / "particle_positions.xlsx") as sheet:
        table = ET.fromstring(sheet.read("xl/worksheets/sheet1.xml"))
    data = [
        {cell.attrib["r"][0]: float(cell.find("x:v", namespace).text) for cell in row}
        for row in table.findall(".//x:row", namespace)[1:]
    ]
    # Table 2. All three directions have the same voxel size.
    pixel_m = .025 / 128
    diameter_px, row_center, column_center = (data[0][key] for key in "FGH")
    tube_radius_m = diameter_px / 2 * pixel_m
    # x is axial in the project. y and z are measured from the tube center.
    spheres = [
        (row["C"] * pixel_m, (row["A"] - row_center) * pixel_m,
         (row["B"] - column_center) * pixel_m, row["D"] * pixel_m)
        for row in data
    ]
    with (output / "sphere_coordinates.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("x_axial_m", "y_from_axis_m", "z_from_axis_m", "radius_m"))
        writer.writerows(spheres)
    wall_overlaps = [hypot(y, z) + radius - tube_radius_m for _, y, z, radius in spheres]
    geometry = {
        "data_doi": "10.4121/4b8d0613-5ee7-4971-b695-835f287336e8",
        "entry_crc32": "1148ff19",
        "voxel_m": pixel_m,
        "nominal_tube_length_m": .061,
        "nominal_particle_diameter_m": .003,
        "tube_radius_m": tube_radius_m,
        "sphere_count": len(spheres),
        "axial_sphere_centers_m": [min(x for x, *_ in spheres), max(x for x, *_ in spheres)],
        "axial_sphere_extents_m": [min(x - r for x, _, _, r in spheres),
                                   max(x + r for x, _, _, r in spheres)],
        "figure_2_axial_average_window_m": [.015875, .040875],
        "sphere_count_intersecting_tube_wall": sum(depth > 0 for depth in wall_overlaps),
        "maximum_tube_wall_overlap_m": max(0, *wall_overlaps),
    }
    (output / "geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--archive", type=Path, help="Use a locally downloaded SupplementaryMaterial.zip")
    mode.add_argument("--local", action="store_true", help="Convert coordinates from already extracted files")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.archive:
        with ZipFile(args.archive) as archive:
            read_entries(archive, args.output)
    elif not args.local:
        download_entries(args.output)
    convert_coordinates(args.output)


if __name__ == "__main__":
    main()
