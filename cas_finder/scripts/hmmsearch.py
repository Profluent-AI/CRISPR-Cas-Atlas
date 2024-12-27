#!/usr/bin/env python3

import pyhmmer
import sys

def run_hmmsearch(seq_path, hmm_path, out_path, size, cpus=0):
        
    # inputs
    hmm_file = pyhmmer.plan7.HMMFile(hmm_path)    
    seq_file = pyhmmer.easel.SequenceFile(seq_path, digital=True)
    
    # run search
    hits = pyhmmer.hmmer.hmmsearch(hmm_file, seq_file, cpus=cpus)
    
    # write results
    with open(out_path, "wb") as out:
        for hit in hits:
            hit.write(out, header=False, format="domains")

if __name__ == "__main__":
    # out_path = "test/output/faa/crispr_contigs.faa.domtblout"
    # hmm_path = "cas-finder-db/profiles.hmm"
    # seq_path = "test/output/faa/crispr_contigs.faa"
    out_path = sys.argv[1]
    hmm_path = sys.argv[2]
    seq_path = sys.argv[3]
    size = 1000
    run_hmmsearch(seq_path, hmm_path, out_path, size)