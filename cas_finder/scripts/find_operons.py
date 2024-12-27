#!/usr/bin/env python3

import time
import os
import sys
import csv
import json
import argparse
import Bio.SeqIO
from collections import defaultdict
from operator import itemgetter


class Operon:
    def __init__(self, id, contig, assembly):
        self.id = id
        self.contig_id = contig
        self.contig_alias = contigs[contig]["alias"]
        self.assembly_id = assembly
        self.coords = []
        self.unannotated = 0

    def update(self, coord):
        self.coords.append(coord)
        if coord[-2] == "protein" and proteins[coord[-1]]["hmm_aln"] is None:
            self.unannotated += 1
        else:
            self.unannotated = 0

    def end_operon(self, coord, max_unknown_genes=3, max_bp=5000):
        if len(self.coords) == 0:
            return False
        if self.unannotated >= max_unknown_genes:
            return True
        if (coord[0] - self.coords[-1][1]) >= max_bp:
            return True
        return False

    def trim_unannotated(self):

        index = None
        for index in range(len(self.coords)):
            coord = self.coords[index]
            if coord[-2] == "crispr":
                index += 1
                break
            elif coord[-2] == "protein" and proteins[coord[-1]]["hmm_aln"] is not None:
                index += 1
                break
        self.coords = self.coords[index - 1 :]

        index = None
        for index in range(len(self.coords))[::-1]:
            coord = self.coords[index]
            if coord[-2] == "crispr":
                index += 1
                break
            elif coord[-2] == "protein" and proteins[coord[-1]]["hmm_aln"] is not None:
                index += 1
                break
        self.coords = self.coords[:index]

    def add_flanking(self, contigs, n_bp):

        self.flanking_start = max([1, self.start_pos - n_bp])
        self.flanking_end = min(
            [self.end_pos + n_bp, len(contigs[self.contig_id]["sequence"])]
        )

        self.flanking_seq = contigs[self.contig_id]["sequence"][
            self.flanking_start - 1 : self.flanking_end
        ]
        self.flanking_masked = contigs[self.contig_id]["masked"][
            self.flanking_start - 1 : self.flanking_end
        ]

        self.max_intergenic = max([len(_) for _ in self.flanking_masked.split("N")])
        self.num_intergenic = len(
            [_ for _ in self.flanking_masked.split("N") if len(_) > 0]
        )

    def summarize(self):

        self.start_pos = min([_[0] for _ in self.coords])  # self.coords[0][0]
        self.end_pos = max([_[1] for _ in self.coords])  # self.coords[-1][1]
        self.length = self.end_pos - self.start_pos + 1

        self.genes = [proteins[_[-1]] for _ in self.coords if _[-2] == "protein"]
        self.protein_ids = [_[-1] for _ in self.coords if _[-2] == "protein"]
        self.n_proteins = len(self.protein_ids)

        self.crispr_ids = [_[-1] for _ in self.coords if _[-2] == "crispr"]
        self.n_crisprs = len(self.crispr_ids)
        #self.crispr_id = (
        #    self.crispr_ids[0] if self.n_crisprs > 0 else None
        #)  # need better way of choosing here

        self.cas_genes = []
        self.effectors = []
        self.cas_genes = []
        for protein_id in self.protein_ids:
            hmm_aln = proteins[protein_id]["hmm_aln"]
            if hmm_aln is None:
                continue
            self.cas_genes.append(hmm_aln["gene_name"])
            if (
                hmm_aln["gene_name"].startswith("Cas9")
                or hmm_aln["gene_name"].startswith("Cas12")
                or hmm_aln["gene_name"].startswith("Cas13")
            ):
                self.effectors.append(protein_id)
        self.n_cas_genes = len(self.cas_genes)
        self.n_effectors = len(self.effectors)        
        self.effector_id = (
            self.effectors[0] if self.n_effectors > 0 else None
        )  # need better way of choosing here

    def print_summary(self):
        print(f"operon: {self.id}")
        print(f"length: {self.length}")
        print(f"n_proteins: {self.n_proteins}")
        print(f"n_cas_genes: {self.n_cas_genes}")
        print(f"n_effectors: {self.n_effectors}")
        print(f"n_crisprs: {self.n_crisprs}")
        print(f"genes: {self.cas_genes}")
        print()

    def keep(self):
        if self.n_crisprs >= 1 and self.n_cas_genes >= 1:
            return True
        else:
            return False

    def make_json(self):
        self.json = {
            "metadata": {
                "source_db": None,
                "assembly_id": self.assembly_id,
                "assembly_type": None,
                "contig_id": self.contig_id,
                "contig_alias": self.contig_alias,
                "contig_length": len(contigs[self.contig_id]["sequence"]),
            },
            "operon": {
                "n_proteins": self.n_proteins,
                "n_crisprs": self.n_crisprs,
                "n_effectors": self.n_effectors,
                "n_cas_genes": self.n_cas_genes,
                "cas_genes": self.cas_genes,
                "length": self.length,
                "start_pos": self.start_pos,
                "end_pos": self.end_pos,                
            },
            "tracr": None,
            "crisprs": [crisprs[_].json for _ in self.crispr_ids],
            "genes": [protein_to_json(_) for _ in self.genes],
            "effector": protein_to_json(proteins[self.effector_id])
                        if self.effector_id
                        else None,
            "genomic": {
                "sequence": self.flanking_seq,
                "masked": self.flanking_masked,
                "start_pos": self.flanking_start,
                "end_pos": self.flanking_end,
                "length": len(self.flanking_seq),
                "max_intergenic": self.max_intergenic,
                "num_intergenic": self.num_intergenic,
            },
        }

    def sanity_check(self):

        # make sure flanking coords reflect true length
        flanking_length = (
            self.json["genomic"]["end_pos"] - self.json["genomic"]["start_pos"] + 1
        )
        assert len(self.json["genomic"]["sequence"]) == flanking_length

        # make sure extracted nucleotide seqs match gene coordinates
        for gene in self.json["genes"]:

            start = gene["start_pos"] - self.json["genomic"]["start_pos"]
            end = gene["end_pos"] - self.json["genomic"]["start_pos"]
            masked = self.json["genomic"]["masked"][start : end + 1]
            assert masked == len(masked) * "N"

            # make sure gene coords match up with ATG start
            if (
                gene["strand"] == 1
                and gene["partial"] == "00"
                and gene["start_type"] == "ATG"
            ):

                atg = self.json["genomic"]["sequence"][start : start + 3]
                assert atg == "ATG"

            # make sure gene length is divisible by 3
            assert int((end - start + 1) / 3) == gene["length"]

            # make sure the gene never extends past the contig
            assert end < len(self.json["genomic"]["sequence"])


def protein_to_json(r):
    x = {
        "gene_id": r["gene_id"],
        "start_pos": r["start_pos"],
        "end_pos": r["end_pos"],
        "strand": r["strand"],
        "genetic_code": r["genetic_code"],
        "length": r["length"],
        "partial": r["info"]["partial"],
        "start_type": r["info"]["start_type"],
        "hmm_align": r["hmm_aln"],
        "protein": r["protein"].rstrip("*"),
        "dna": r["dna"],
    }
    if x["hmm_align"]:
        if "seq_name" in x["hmm_align"]:
            del x["hmm_align"]["seq_name"]
        if "seq_name" in x["hmm_align"]:
            del x["hmm_align"]["tacc"]
        # x["hmm_align"]["hmm_id"] = x["hmm_align"]["hmm_name"]
        if "seq_name" in x["hmm_align"]:
            del x["hmm_align"]["hmm_name"]
    return x


class CRISPR:
    def __init__(self, r):
        self.id = f"{r['contig_id']}_{r['array_num']}"
        self.contig_id = r["contig_id"]
        self.start_pos = int(r["start_pos"])
        self.end_pos = int(r["end_pos"])
        self.n_spacers = int(r["num_spacers"])
        self.repeat_identity = float(r["percent_id"])
        self.spacer_length = float(r["spacer_length"])
        self.repeat_length = int(r["repeat_length"])
        self.consensus_repeat = r["consensus_repeat"]
        self.partial = (
            "00"
            if r["truncated"] == "neither"
            else "01"
            if r["truncated"] == "right"
            else "10"
        )
        self.spacers = []
        self.subtype = None
        self.subtype_probability = None

    def make_json(self):
        self.json = {
            "crispr_id": self.id,
            "start_pos": self.start_pos,
            "end_pos": self.end_pos,
            "n_spacers": self.n_spacers,
            "repeat_identity": self.repeat_identity,
            "spacer_length": self.spacer_length,
            "repeat_length": self.repeat_length,
            "consensus_repeat": self.consensus_repeat,
            "partial": self.partial,
            "subtype": self.subtype,
            "subtype_probability": self.subtype_probability,
            "spacers": self.spacers,
        }


def parse_prodigal(path):
    for r in Bio.SeqIO.parse(path, "fasta"):
        x = {}
        x["id"] = r.id
        x["contig_id"] = r.id.rsplit(":", 1)[0]
        x["offset"] = int(r.id.rsplit(":", 1)[1].split("-")[0]) - 1
        x["seq"] = str(r.seq)
        x["start_pos"] = int(r.description.split()[2]) + x["offset"]
        x["end_pos"] = int(r.description.split()[4]) + x["offset"]
        x["strand"] = int(r.description.split()[6])
        x["info"] = dict([_.split("=") for _ in r.description.split()[-1].split(";")])
        x["genetic_code"] = x["info"]["genetic_code"]
        x["length"] = len(x["seq"])
        x["gene_id"] = f"{x['contig_id']}_{x['start_pos']}"
        x["hmm_aln"] = None
        yield x


def parse_args():
    parser = argparse.ArgumentParser(
        description="identify crispr operons from search results"
    )
    parser.add_argument(
        "dir",
        type=str,
        help="Path to CRISPR directory containing temporary files",
    )
    parser.add_argument(
        "sample",
        type=str,
        help="Sample identifier",
    )
    parser.add_argument(
        "--db",
        type=str,
        help="Path to database directory",
        required=True,
    )    
    parser.add_argument(
        "--max_unknown",
        type=int,
        help="Max allowed number of unknown genes between elements in CRISPR operon (3)",
        default=3,
    )
    parser.add_argument(
        "--max_bp",
        type=int,
        help="Max allowed number of bp between elements in CRISPR operon (5000)",
        default=5000,
    )
    parser.add_argument(
        "--flanking",
        type=int,
        help="Number of bp flanking the operon (500)",
        default=500,
    )
    args = vars(parser.parse_args())
    args["dir"] = args["dir"].rstrip("/")
    return args

if __name__ == "__main__":
    
    # parse args
    args = parse_args()

    # initialize contigs
    contigs = {}
    p = f"{args['dir']}/tmp/fna/pruned_contigs.fna"
    for r in Bio.SeqIO.parse(p, "fasta"):
        contigs[r.id] = {
            "sequence": str(r.seq).upper(),
            "masked": str(r.seq).upper(),
            "alias": str(r.description.split()[-1]),
        }

    # initialize proteins
    proteins = {}
    p = f"{args['dir']}/tmp/faa/crispr_contigs.faa"
    for r in parse_prodigal(p):
        r["protein"] = r["seq"]
        del r["seq"]
        proteins[r["id"]] = r
    p = f"{args['dir']}/tmp/ffn/crispr_contigs.ffn"
    for r in parse_prodigal(p):
        proteins[r["id"]]["dna"] = r["seq"]

    # annotate proteins
    p = f"{args['dir']}/tmp/hmm/hmm_full.tsv"
    for r in csv.DictReader(open(p), delimiter="\t"):
        # format record
        r["evalue"] = float(r["evalue"])
        r["score"] = float(r["score"])
        r["seq_cov"] = float(r["seq_cov"])
        r["hmm_cov"] = float(r["hmm_cov"])
        r["gene_name"] = r["hmm_name"].split("_")[0]
        # filter hmm results
        if float(r["evalue"]) > 1e-5:
            continue
        elif float(r["seq_cov"]) < 0.25 and float(r["hmm_cov"]) < 0.25:
            continue
        # store best hits
        if proteins[r["seq_name"]]["hmm_aln"] is None:
            proteins[r["seq_name"]]["hmm_aln"] = r
        elif r["score"] > proteins[r["seq_name"]]["hmm_aln"]["score"]:
            proteins[r["seq_name"]]["hmm_aln"] = r

    # store crisprs
    crisprs = {}
    p = f"{args['dir']}/tmp/crispr/crispr_arrays.tsv"
    for r in csv.DictReader(open(p), delimiter="\t"):
        crispr = CRISPR(r)
        crisprs[crispr.id] = crispr

    # add spacers to crispr arrays
    p = f"{args['dir']}/tmp/crispr/crispr_spacers.tsv"
    for r in csv.DictReader(open(p), delimiter="\t"):
        crispr_id = f"{r['contig_id']}_{r['array_num']}"
        for field in ["contig_id", "array_num"]:
            del r[field]
        r["start_pos"] = int(r["start_pos"])
        r["repeat_num"] = int(r["repeat_num"])
        crisprs[crispr_id].spacers.append(r)


    # make crispr json
    for crispr in crisprs.values():
        crispr.make_json()

    # create sorted list of features (proteins and crisprs)
    feature_coords = defaultdict(list)
    for r in proteins.values():
        coord = int(r["start_pos"]), int(r["end_pos"]), "protein", r["id"]
        feature_coords[r["contig_id"]].append(coord)
    for r in crisprs.values():
        coord = r.start_pos, r.end_pos, "crispr", r.id
        feature_coords[r.contig_id].append(coord)
    for contig in feature_coords:
        feature_coords[contig] = sorted(feature_coords[contig], key=itemgetter(0, 1))

    # mask all features from contigs with Ns
    for contig in feature_coords:
        for coords in feature_coords[contig]:
            start_pos, end_pos, feat_type, protein_id = coords
            masked = contigs[contig]["masked"]
            contigs[contig]["masked"] = (
                masked[0 : start_pos - 1]
                + "N" * (end_pos - start_pos + 1)
                + masked[end_pos:]
            )

    # group features into operons
    assembly = args["sample"]  # os.path.basename(args["dir"]).rsplit(".fna", 1)[0]
    operons = {}
    operon_num = 0
    for contig, coords in feature_coords.items():
        operon_num += 1
        operon_id = f"{assembly}@{operon_num}"
        operons[operon_id] = Operon(operon_id, contig, assembly)
        for coord in coords:
            if operons[operon_id].end_operon(
                coord, args["max_unknown"], args["max_bp"]
            ):
                operon_num += 1
                operon_id = f"{assembly}@{operon_num}"
                operons[operon_id] = Operon(operon_id, contig, assembly)
            operons[operon_id].update(coord)

    # trim operons
    for id, operon in operons.items():
        operon.trim_unannotated()
    operons = dict(
        [(id, operon) for id, operon in operons.items() if len(operon.coords) > 0]
    )

    # finalize
    for id, operon in operons.items():
        operon.summarize()
        operon.add_flanking(contigs, n_bp=args["flanking"])
        operon.make_json()
        operon.sanity_check()

    # write json formatted operons as output
    document = dict(
        [(id, operon.json) for id, operon in operons.items() if operon.keep()]
    )
    p = os.path.join(args["dir"], "crispr-cas.json")
    with open(p, "w") as output:
        json.dump(document, output, indent=6)
