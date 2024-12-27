#!/usr/bin/env python3

import sys
import csv
import Bio.SeqIO
import argparse

parser = argparse.ArgumentParser(description="extract crispr regions from fasta file")
parser.add_argument(
    "fna",
    type=str,
    help="Path to contigs FASTA file",
)
parser.add_argument(
    "crispr",
    type=str,
    help="Path to crispr array TSV file",
)
parser.add_argument(
    "--flanking",
    type=int,
    help="Number of bp flanking each CRISPR (20000)",
    default=20000,
)
parser.add_argument(
    "--min_spacers",
    type=int,
    help="Min # spacers per array (3)",
    default=3,
)
args = vars(parser.parse_args())

seqs = dict([[_.id, str(_.seq)] for _ in Bio.SeqIO.parse(args["fna"], "fasta")])

contig_coords = {}
for r in csv.DictReader(open(args["crispr"]), delimiter="\t"):
    if int(r["num_spacers"]) < args["min_spacers"]:
        continue
    contig, start, end = r["contig_id"], int(r["start_pos"]), int(r["end_pos"])
    start = max([1, start - args["flanking"]])
    end = min([len(seqs[contig]), end + args["flanking"]])
    if contig not in contig_coords:
        contig_coords[contig] = [[start, end]]
    else:
        contig_coords[contig].append([start, end])

for contig, coords in contig_coords.items():
    coords = sorted(coords)
    updated = [coords[0]]
    if len(coords) > 1:
        for start, end in coords[1:]:
            if start <= updated[-1][-1]:
                updated[-1][-1] = max([updated[-1][-1], end])
            else:
                updated.append([start, end])
    contig_coords[contig] = updated

for contig, coords in contig_coords.items():
    for index, coord in enumerate(coords):
        start_pos, end_pos = coord
        contig_length = len(seqs[contig])
        sub_seq = seqs[contig][start_pos - 1 : end_pos]
        sys.stdout.write(
            f">{contig}:{start_pos}-{end_pos}/{contig_length}\n{sub_seq}\n"
        )
