# MOBAseq

## Installation

```bash
pip install git+https://github.com/RuiTang-2024/MOBASeq-Toolkit.git@main

pip install git+https://github.com/IrenaeusChan/mobaseq_dev.git@main
```

OR

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -e .
```

## Usage

```bash
mobaseq --help
```

Definitions:
Metastatic Burden = (Total Cell Number in Tissue / Total Cell Number prior to Transplantation) / (Total Cell Number for sgSafe Ctrl / Total Cell Number for sgSafe Ctrl prior to Transplantation)

Metatstatic Seeding = (Total Colony Number in Tissue / Total Cell Number prior to Transplantation) / (Total Colony Number for sgSafe Ctrl / Total Cell Number for sgSafe Ctrl prior to Transplantation)

Peak Mode = How many doublings needed to reach the current Cell Number, i.e. 2^PeakMode = Cell Number

Distance/Outlier Score = Ideally if MetBurden and MetSeeding are 1, then the Distance/Outlier Score should be 0. If the Distance/Outlier Score is greater than 0, then the sample is an outlier.

If MetBurden > MetSeeding this indicates that this sgID is a Tumor Suppressor
If MetBurden < MetSeeding this indicates that this sgID is a Oncogene/Tumor Promoter

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
