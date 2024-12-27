# CRISPR-Cas Atlas

This repository contains documentation regarding the CRISPR-Cas Atlas. The database is described in detail in **Ruffolo, Nayfach, Gallagher, and Bhatnagar et al.,** (2024). *Design of highly functional genome editors by modeling the universe of CRISPR-Cas sequences.* bioRxiv. [https://doi.org/2024.04.22.590591v1](https://www.biorxiv.org/content/10.1101/2024.04.22.590591v1).

More information about OpenCRISPR can be found at https://github.com/Profluent-AI/OpenCRISPR.

## Quickstart

Download the database from Google Cloud Storage:
```bash
wget https://storage.googleapis.com/profluent-public/crispr-cas-atlas/crispr-cas-atlas-v1.0.json
```

## Data format
The CRISPR-Cas Atlas is formatted as a JSON document. 

An example record can be found in the file `example.json` and is shown below. The record represents a Type II-A CRISPR-Cas operon identified from an NCBI metagenome-assembled genome containing Cas1, Cas2, Cas9, a tracrRNA, and a CRISPR array:
```json
{
  "operon_id": "GCA_947475615.1@2",
  "summary": {
   "subtype": "II-A",
   "subtype_score": "II-A",
   "operon_length": 5576,
   "n_crispr": 1,
   "n_spacers": 9,
   "n_tracr": 1,
   "n_cas": 3,
   "n_genes": 5
  },
  "metadata": {
   "source_db": "NCBI",
   "assembly_type": "MAG",
   "biosample_id": "SAMEA112175321",
   "sample_name": "Comamonadaceae bacterium",
   "taxonomy": "d__Bacteria;p__Pseudomonadota;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__;s__Comamonadaceae bacterium",
   "biome": null
  },
  "crispr": [
   {
    "crispr_repeat": "GTTCCGGCCAGAGCGCATTTCCCAATCAAATAGACT",
    "crispr_spacers": [
     "TGAAGAAATATGCGAATGTGAAAGCGAATA",
     "CAAGTGAAACCTATACAGGTTAAACAACAG",
     "AAAGCCGGTGGTGGATAGCGCCTCAAGCGC",
     "TGAACTTTCACGCCCACCTATAGGCAATCC",
     "GGTACCTTTGCGGTGGACTCCATGATGTGG",
     "TTTGCTTGCGTCTCAAAAGCTGGCGATCAA",
     "CCGATGACAGTGAGCCAAGCTGCAAATACG",
     "CTGTGCCGCCCGCTTGAATTGCGGCAAGCG",
     ""
    ]
   }
  ],
  "tracr": {
   "cm_id": "Cluster_1494",
   "evalue": 4.7e-07,
   "truncated": "00",
   "gene_overlap": "00",
   "terminator": 1,
   "confidence": "High",
   "seq_unmasked": "GATTGGGAAATGCGCTCTGGACGCTAACAAGCAGATGACTTGCAAAAGTCTGGATGCACAAAATGAAGAGGCCGCTATATGCGGCCTCTTGTCTTTTCAGA"
  },
  "cas": [
   {
    "gene_name": "Cas2",
    "hmm_name": "Cas2_5_CAS-I-II-III-IV-V-VI",
    "evalue": 5e-31,
    "score": 93.6,
    "truncated": "00",
    "length": 117,
    "protein": "MSRRAKTSLSGYRIMWMLVMFDLPVVTASERLAANQFRHSLLDMGFLRCQLSVYMRFCTSAAQVQTYCQRVEAALPNGGQVNIMQLTDKQFERVISFQGRKAQPAKKTPDQFDLFD"
   },
   {
    "gene_name": "Cas1",
    "hmm_name": "Cas1_4_CAS-I-II-III-IV-V-VI",
    "evalue": 0.0,
    "score": 229.1,
    "truncated": "00",
    "length": 309,
    "protein": "MLGRIVEVANDKRHLSMYRGFMLVQSTGEDRQEVGRVALDDMSALIANAHGLSYTNNLLVALAERGVPMVLCAANHNVVGMLWPAEGHHQQAHRMEAQIACSLPTRKRLWAAIVKSKLLNQAAVLAAAGAPAAPLQMLARQVKSGDPQNTEAQGARKYWGLLMGPLFRRDQQADGLNALLNYGYTVLRAATARAVVAAGLHPSVGLHHSHDNNAMRLVDDVMEPFRPVIDWTVWQLQSQGPCVVNADTKRALVQSLYQDLQSDAGTTPVLVAVQKLATSLAQVMLGERDKLDLPHAGVPQRYTESDDE"
   },
   {
    "gene_name": "Cas9",
    "hmm_name": "Cas9_c4",
    "evalue": 0.0,
    "score": 1215.5,
    "truncated": "00",
    "length": 1016,
    "protein": "MHMTKMRYRLALDLGSTSLGWAMLRLNVNNEPSAVIKAGVRIFSDGRNPKDGASLAVSRREARAMRRRRDRLLKRKARMMRTLLVHGFFPHDLAARKALERLEPLSLRAKGLDQTLQPAEFARALFHINQRRGFKSNRKTDKKEVDSSALKNAIGQLREAMQATGCRTVGEWLYARHQKGLPIRARYRENRSTRDDGKTKIEKSYDLYIDRAMIEAEFDALWAKQAELNPVQFHETARVEIKDCLLHQRRLKPVKPGRCTLIPEEERAPLALPSQQRFRIYQEVNNLRLIREGLTEDPLTPAQRDQLVQALETKSKVTFAQIKKVLGFSGQFNLEDDKRTELKGNATSTSLSKKEHFGSAWAGMDAAQQDSIVLQLLTEENEATLIQWLKSATGVDEITAERIANAALPEGYGSLSAKALDKILPELRREVVTFDKAVIAAGFDHHSHLSHAVTGEILPALPYYGEYLQRHVGFGSGKPEDPAEKRFGKIANPTVHIGLNQVRIVVNALIKRYGHPSEVIVEVARDLKQSQEQRKDDQKRQADNQHRNARMREQIADLLNTSPERVQTTDLHKMILWEELNRDNAADRRCPYSGAQISAAMLFSDQVEIEHILPFSQTLDDSLNNKTVALRQANRIKGNRTPWQARDDFSAQGWVIVDMLARAELMPKNKRYRFGENGYAQWLRDDKGFLARALNDTRHLSRVAREYLSLICPQNTRAIPGQMTAMLRAKFGLNNILGLNGEKNRNDHRHHAVDACVIAVTDQGMLQRFASASASAREQQLNKLVDTMPLPWESYREHVKRAVDNIWVSHKPDHGHEGAMHNDTAYGLLGKDRVHVRKVVDGQRVRKESTLKVIPFSDAKASARHGLLPDGQPRPYKGYKGDSNYCIEIVRNDKGKWEGEVISTFEAYQLVRQGGVQRLRHPTLSCSGKPLVMRLMIDDSVVILIDDVKHVLRLAYMASAGTMAFAPCNEANVDKRTRTKEMAYTFKTAGSLQKAKGRRISISPIGELRDPGFRD"
   }
  ]
 }
```

## CRISPR-Cas Prediction

To access the pipeline used to construct the CRISPR-Cas Atlas, please see this [README](docs/CAS_FINDER.md).

## License

This data and code is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** license.

**Usage Restriction**:  
This data and code is provided for *academic and non-commercial purposes only*. Any commercial use is strictly prohibited without prior permission.

You can view the full license [here](https://creativecommons.org/licenses/by-nc/4.0/).

## Cite this work

If you use the CRISPR-Cas Atlas or search tool in your research, please cite the following preprint:

```bibtex
@article{profluent2024opencrispr,
  title={Design of highly functional genome editors by modeling the universe of CRISPR-Cas sequences},
  author={Ruffolo, Jeffrey A and Nayfach, Stephen and Gallagher, Joseph and Bhatnagar, Aadyot and Beazer, Joel and Hussain, Riffat and Russ, Jordan and Yip, Jennifer and Hill, Emily and Pacesa, Martin and others},
  journal={bioRxiv},
  pages={2024--04},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

Also consider citing [these resources and tools](/docs/REFERENCES.md) utilized by the CRISPR-Cas Atlas.
