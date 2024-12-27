#!/usr/bin/env python3

import sys
import re


def parse_domtblout(path):
    # http://eddylab.org/software/hmmer/Userguide.pdf (Page 70: The domain hits table)
    for line in open(path):
        if line.startswith("#"):
            continue
        r = re.split(r"\s+", line)
        d = {}
        d["seq_name"] = r[0]
        d["seq_acc"] = r[1]
        d["seq_len"] = int(r[2])
        d["hmm_name"] = r[3]
        d["hmm_acc"] = r[4]
        d["hmm_len"] = int(r[5])
        d["evalue_full"] = float(r[6])
        d["score_full"] = float(r[7])
        d["bias_full"] = r[8]
        d["dom_num"] = r[9]
        d["dom_total"] = r[10]
        d["evalue_domain"] = float(r[11])
        # d["evalue_i"] = float(r[12])
        d["score_domain"] = float(r[13])
        d["bias_domain"] = float(r[14])
        d["hmm_start"], d["hmm_end"] = int(r[15]), int(r[16])
        d["seq_start"], d["seq_end"] = int(r[17]), int(r[18])
        d["env_start"], d["env_end"] = int(r[19]), int(r[20])
        d["seq_aln"] = d["seq_end"] - d["seq_start"] + 1
        d["hmm_aln"] = d["hmm_end"] - d["hmm_start"] + 1
        d["seq_cov"] = round(1.0 * d["seq_aln"] / d["seq_len"], 3)
        d["hmm_cov"] = round(1.0 * d["hmm_aln"] / d["hmm_len"], 3)
        yield d


fields = ["seq_name", "hmm_name", "hmm_acc", "seq_len", "hmm_len", "seq_cov", "hmm_cov", "evalue", "score"]
sys.stdout.write("\t".join(fields) + "\n")
for r in parse_domtblout(sys.argv[1]):
    r["evalue"] = r["evalue_domain"]
    r["score"] = r["score_domain"]
    row = "\t".join([str(r[f]) for f in fields]) + "\n"
    sys.stdout.write(row)
