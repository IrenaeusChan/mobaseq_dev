# MOBAseq

## Installation

```bash
pip install git+https://github.com/RuiTang-2024/MOBASeq-Toolkit.git@main
```

OR

```bash
python3 -m venv .venv
source .venv/bin/activate
pip3 install -e .
```

## Usage

```bash
mobaseq --help
```

Definitions:
Metastatic Burden = (Total Cell Number in Tissue / Total Cell Number prior to Transplantation) / (Total Cell Number for sgSafe Ctrl / Total Cell Number for sgSafe Ctrl prior to Transplantation)

MOBASeq
=======
Major Changes:
    - Python CLI for ease of download and user friendly processing
    - Addressing *BUG* where R1_Barcode == R2_Barcode but R1_sgID =/= R2_sgID
    - Moba500 61 Samples: 31 Minutes to Fully Analyze using 4 Cores

Minor Changes:
    - Changed algorithm for sgID --> GeneName mapping to search through entire file rather than stopping at first "best match"
    - Allow for batch processing or individual file processing, for situations where 1 sample needs to be rerun or fails
    - TODO: Allow users ability to implement their own barcode strategy

Optimizations:
    - Introduced condition if the FASTQ size is greater than 1GB to utilize a large data reading method
    - Utilizing NumPy package to vectorize most data manipulation, decreases overall processing time
    - Implemented machine byte code version of Hamming Distance to improve computational efficiency during Hamming distance calculation
    - Utilized multithreading and parallized sample processing for amount of cores available on computer or cluster to improve processing efficiency
        - E.g. for 5 Samples (of varying sizes):
            - Single Core: 22:21:28 --> 22:55:51 Elapsed Time 34 min 53s
            - Four Cores: 05:54:21 --> 06:10:59 Elapsed Time 16 min 39s
            - Legacy Algorithm: 22:36:04 --> 03:04:35 Elapsed Time 4 Hrs 28 min 31s
    - Utilized multithreading, broadcasting strategy, and parallelized sgID matching to improve processing efficiency
        - E.g. for 25 Samples:
            - Single Core: 04:14:09 --> 04:15:50 Elapsed Time 1 min 41s
            - Four Cores: 04:08:06 --> 04:08:44 Elapsed Time 38s
            - Legacy Algorithm: 03:23:31 --> 04:14:08 Elapsed Time 50 min 37s

Questions:
    - If we dont know where LV binds to the genome, how do we use CRISPR to target genes specifically?
    - Why is blood not filtered for Cell_Num > 50 but all other Tissue types are?
    - What does CTC stand for? SuperMet
      - What's the SuperMet definition again? If detected in the blood and detected in the tissue?
    - What is unique about these sgIDs that they must be removed? 
      - sgSafe13, sgSafe1, Ccnd1 (especially this because it's not a sgSafe Ctrl)
    - Why is Peak Mode defined on a Log2 scale?
    - Why Log2 vs Log10? (Is it because assumed exponential growth?)
    - Where does this gene list come from (for Dormancy)?
      - gene_list <- c("Crebbp", "Pten", "Csnk1a1", "Gpatch8", "Zeb1", "Tsc1", "Tsc2", "Gata6", "Mga","Zmiz1", "Grhpr", "Hnf4a", "Eif5b", "Soga1", "Nup98", "Ctrl")
    - Mode2 for the normalmixEM, why for anything less than 2,000?