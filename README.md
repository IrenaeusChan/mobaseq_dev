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

   Sample_ID Spike_Rsq
1      Brn11 0.9769261
2      Brn12 0.9966988
3      Brn13 0.9784124
4      Brn14 0.9619013
5      Brn15 0.9965447
6      Brn16 0.9916088
7      Brn21 0.8832450
8      Brn22 0.9290073
9      Brn23 0.9808438
10     Brn24 0.9051336
11     Brn25 0.9208960
12     Brn26 0.9259640
13     Brn27 0.9024599
14     Brn28 0.9142950
15     Brn29 0.8775403
16     Liv11 0.8535246
17     Liv12 0.9200064
18     Liv13 0.9623270
19     Liv14 0.9351993
20     Liv15 0.9318214
21     Liv16 0.9199480
22     Liv21 0.9695602
23     Liv22 0.9547380
24     Liv23 0.9721261
25     Liv24 0.9624038
26     Liv25 0.9125539
27     Liv26 0.9948438
28     Liv27 0.9769652
29     Liv28 0.9582834
30     Liv29 0.9790924
31     Lug11 0.9860086
32     Lug12 0.9874371
33     Lug13 0.9772633
34     Lug14 0.9941175
35     Lug15 0.9739853
36     Lug16 0.9613365
37     Lug21 0.9977475
38     Lug22 0.9877952
39     Lug23 0.9982292
40     Lug24 0.9464346
41     Lug25 0.9975629
42     Lug26 0.9922634
43     Lug27 0.9907117
44     Lug28 0.9930452
45     Lug29 0.9856332
46  preTran1        NA
47      WB11 0.9926839
48      WB12 0.9839951
49      WB13 0.9977735
50      WB14 0.9667134
51      WB15 0.9864507
52      WB16 0.9572261
53      WB21 0.9833962
54      WB22 0.9828065
55      WB23 0.9893883
56      WB24 0.9944391
57      WB25 0.9744131
58      WB26 0.9930776
59      WB27 0.9616109
60      WB28 0.9927409
61      WB29 0.9839260

Sample_ID	Spike_Rsq
Brn11	0.9485931257978324
Brn12	0.9926452302292328
Brn13	0.9519045783195405
Brn14	0.915119149
Brn15	0.9923019732308136
Brn16	0.9813050122641144
Brn21	0.6941788789960557
Brn22	0.8252934217745455
Brn23	0.7522827599297062
Brn24	0.7886452833582878
Brn25	0.8053323374843424
Brn26	0.8350536836040455
Brn27	0.7542822899477825
Brn28	0.8115398710187707
Brn29	0.6305405341232551
Liv11	0.4968616747668642
Liv12	0.7396728693551252
Liv13	-0.247659975
Liv14	0.8556291774105078
Liv15	0.7548136565977807
Liv16	0.80299939
Liv21	0.93218246
Liv22	0.8886144320296135
Liv23	0.9378992773815572
Liv24	0.9162385835696266
Liv25	0.4016991982220295
Liv26	0.9885123698807936
Liv27	0.9486802560132211
Liv28	0.6582210400833168
Liv29	0.9534194832380334
Lug11	0.9688283654066304
Lug12	0.97201079
Lug13	0.949344482
Lug14	0.9868943246118126
Lug15	0.9420412646980546
Lug16	0.913860882
Lug21	0.9949816665160838
Lug22	0.9728087541786647
Lug23	0.9960548862659636
Lug24	0.8806605637194967
Lug25	0.994570388
Lug26	0.982763445
Lug27	0.9793063202340196
Lug28	0.9845052621378534
Lug29	0.9679920244366176
WB11	0.9837004047130088
WB12	0.9643424371917906
WB13	0.9950394462919014
WB14	0.9258399976773944
WB15	0.969813289
WB16	0.9047032270312498
WB21	0.9630079697280616
WB22	0.9616942937928272
WB23	0.9763579866785582
WB24	0.9876107097244744
WB25	0.9429943252716854
WB26	0.9845774048135738
WB27	0.9144721749749802
WB28	0.9838273891500924
WB29	0.9641884374723289
preTran1	0