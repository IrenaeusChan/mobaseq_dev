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

MOBASeq
=======
Major Changes:
    - Python CLI for ease of download and user friendly processing
    - Addressing *BUG* where R1_Barcode == R2_Barcode but R1_sgID =/= R2_sgID
    - TODO: Consideration for removing "fuzzy" matching for improved processing speed and statistical power

Minor Changes:
    - Changed algorithm for sgID --> GeneName mapping to search through entire file rather than stopping at first "best match"
    - Allow for batch processing or individual file processing, for situations where 1 sample needs to be rerun or fails
    - TODO: Allow users ability to implement their own barcode strategy

Optimizations:
    - Introduced condition if the FASTQ size is greater than 1GB to utilize a large data reading method
    - Utilizing NumPy package to vectorize most data manipulation, decreases overall processing time (e.g. 61 MOBASeq500 Samples ~2Hrs)
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
