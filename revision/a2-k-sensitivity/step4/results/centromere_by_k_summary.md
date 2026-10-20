# Centromere subtype clustering at k = 21, 31, 51 -- summary

Pipeline: KaryoScope-analysis@pangenome_structure_V2 (90e2874), unchanged, gate (0.5, 5) with its own --exclude-features default (novel) and its own EXCLUDED_CHROMS-before-gating order. Inputs: the published centromere windows re-annotated with the v1-like region rebuilds HKS_human_CHM13_v1like_region_k{21,31,51} (rows v1_*; the v2 six-k family is kept as rows v2_* for the database-generation contrast only -- its vocabulary differs from v1 so it is not comparable to published). (region set, presmoothed, 100 kb bins, trimmed to centromere-class extents as the published script did). 'published' = the v1.1.0 inputs (A6 run, job 37136128) re-gated identically.

## Per k: haplotypes, minor calls, chromosomes retaining a split

| k | haplotypes | chromosomes | minor | major | chromosomes with a split | collapsed |
|---|---|---|---|---|---|---|
| v1_21 | 4551 | 18 | 395 | 4156 | 16 | chr4,chr20 |
| v1_31 | 4551 | 18 | 821 | 3730 | 15 | chr6,chr17,chr20 |
| v1_51 | 4551 | 18 | 1266 | 3285 | 17 | chr4 |
| published | 4551 | 18 | 821 | 3730 | 15 | chr6,chr17,chr20 |
| v2_21 | 4551 | 18 | 722 | 3829 | 15 | chr1,chr6,chr8 |
| v2_31 | 4551 | 18 | 1195 | 3356 | 16 | chr4,chr17 |
| v2_51 | 4551 | 18 | 1151 | 3400 | 16 | chr4,chr6 |

## Major / minor counts per chromosome

| chromosome | k=v1_21 major/minor | k=v1_31 major/minor | k=v1_51 major/minor | k=published major/minor | k2 silhouette 21/31/51/pub |
|---|---|---|---|---|---|
| chr1 | 138/1 | 138/1 | 138/1 | 138/1 | 0.196/0.156/0.274/0.156 |
| chr2 | 43/163 | 205/1 | 204/2 | 205/1 | 0.643/0.184/0.803/0.184 |
| chr3 | 176/2 | 54/124 | 102/76 | 54/124 | 0.827/0.852/0.883/0.852 |
| chr4 | 166/0 | 164/2 | 166/0 | 164/2 | 0.311/0.346/0.381/0.346 |
| chr5 | 355/3 | 354/4 | 340/18 | 354/4 | 0.486/0.417/0.626/0.417 |
| chr6 | 226/4 | 230/0 | 52/178 | 230/0 | 0.347/0.386/0.567/0.386 |
| chr7 | 262/6 | 265/3 | 266/2 | 265/3 | 0.200/0.312/0.315/0.312 |
| chr8 | 341/1 | 136/206 | 142/200 | 136/206 | 0.377/0.694/0.705/0.694 |
| chr9 | 244/2 | 244/2 | 244/2 | 244/2 | 0.411/0.316/0.369/0.316 |
| chr10 | 283/2 | 282/3 | 282/3 | 282/3 | 0.476/0.274/0.368/0.274 |
| chr11 | 388/6 | 216/178 | 118/276 | 216/178 | 0.765/0.921/0.533/0.921 |
| chr12 | 371/10 | 148/233 | 138/243 | 148/233 | 0.856/0.879/0.875/0.879 |
| chr16 | 253/2 | 252/3 | 253/2 | 252/3 | 0.229/0.756/0.821/0.756 |
| chr17 | 158/1 | 159/0 | 158/1 | 159/0 | 0.212/0.098/0.223/0.098 |
| chr18 | 35/133 | 154/14 | 129/39 | 154/14 | 0.509/0.922/0.814/0.922 |
| chr19 | 325/3 | 322/6 | 323/5 | 322/6 | 0.308/0.580/0.324/0.580 |
| chr20 | 306/0 | 306/0 | 129/177 | 306/0 | 0.451/0.449/0.583/0.449 |
| chrX | 86/56 | 101/41 | 101/41 | 101/41 | 0.806/1.000/1.000/1.000 |

## Per-haplotype major/minor concordance vs k = 31 (new)

### k = v1_21 vs k = 31(new): 3479/4551 same call = 76.44%  (haplotypes only in k=v1_21: 0; only in k=31: 0)

| chromosome | same | of | % |
|---|---|---|---|
| chr1 | 139 | 139 | 100.0 |
| chr2 | 44 | 206 | 21.4 |
| chr3 | 56 | 178 | 31.5 |
| chr4 | 164 | 166 | 98.8 |
| chr5 | 355 | 358 | 99.2 |
| chr6 | 226 | 230 | 98.3 |
| chr7 | 261 | 268 | 97.4 |
| chr8 | 137 | 342 | 40.1 |
| chr9 | 244 | 246 | 99.2 |
| chr10 | 282 | 285 | 98.9 |
| chr11 | 216 | 394 | 54.8 |
| chr12 | 148 | 381 | 38.8 |
| chr16 | 252 | 255 | 98.8 |
| chr17 | 158 | 159 | 99.4 |
| chr18 | 45 | 168 | 26.8 |
| chr19 | 319 | 328 | 97.3 |
| chr20 | 306 | 306 | 100.0 |
| chrX | 127 | 142 | 89.4 |

Changed calls (1072; full list in centromere_changes_v1_21_vs_v1_31.tsv):

- HG01074#1#CM088790.1 (chr10: major at k=v1_21 -> minor at k=31)
- HG01960#1#CM088649.1 (chr10: minor at k=v1_21 -> major at k=31)
- HG04187#2#CM099760.1 (chr10: major at k=v1_21 -> minor at k=31)
- HG00099#1#CM087322.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00099#2#CM087366.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00126#1#CM090115.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00232#1#CM089997.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00232#2#CM090028.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00235#2#CM094401.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00272#2#CM094221.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00290#1#CM089985.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00320#2#CM089970.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00321#1#CM089546.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00323#2#CM087335.1 (chr11: minor at k=v1_21 -> major at k=31)
- HG00329#2#CM094243.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00408#1#CM085959.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00423#1#CM089405.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00438#2#CM089195.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00544#1#CM089376.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00544#2#CM089391.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00558#2#CM088529.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00597#1#CM085774.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00597#2#CM085788.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00609#1#CM086354.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00639#2#CM088559.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00642#1#CM086612.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00706#1#CM100020.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00706#2#CM100065.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00733#1#CM099633.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00735#2#CM087892.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG00738#2#CM086686.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01074#2#CM088817.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01081#2#CM088588.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01099#2#CM086273.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01109#2#CM092080.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01123#2#CM089102.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01175#1#CM087930.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01192#1#CM085932.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01243#1#CM099608.1 (chr11: major at k=v1_21 -> minor at k=31)
- HG01243#2#CM099596.1 (chr11: major at k=v1_21 -> minor at k=31)
- ... 1032 more

### k = v1_51 vs k = 31(new): 3908/4551 same call = 85.87%  (haplotypes only in k=v1_51: 0; only in k=31: 0)

| chromosome | same | of | % |
|---|---|---|---|
| chr1 | 139 | 139 | 100.0 |
| chr2 | 205 | 206 | 99.5 |
| chr3 | 120 | 178 | 67.4 |
| chr4 | 164 | 166 | 98.8 |
| chr5 | 344 | 358 | 96.1 |
| chr6 | 52 | 230 | 22.6 |
| chr7 | 267 | 268 | 99.6 |
| chr8 | 336 | 342 | 98.2 |
| chr9 | 244 | 246 | 99.2 |
| chr10 | 283 | 285 | 99.3 |
| chr11 | 256 | 394 | 65.0 |
| chr12 | 361 | 381 | 94.8 |
| chr16 | 254 | 255 | 99.6 |
| chr17 | 158 | 159 | 99.4 |
| chr18 | 135 | 168 | 80.4 |
| chr19 | 319 | 328 | 97.3 |
| chr20 | 129 | 306 | 42.2 |
| chrX | 142 | 142 | 100.0 |

Changed calls (643; full list in centromere_changes_v1_51_vs_v1_31.tsv):

- HG01074#1#CM088790.1 (chr10: major at k=v1_51 -> minor at k=31)
- HG03130#1#CM087259.1 (chr10: minor at k=v1_51 -> major at k=31)
- HG00097#2#CM094082.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00128#1#CM090085.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00140#2#CM087141.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00146#1#CM090018.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00272#1#CM094194.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00290#2#CM090008.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00320#1#CM089955.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00321#2#CM089577.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00329#1#CM094231.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00344#2#CM094517.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00350#2#CM094495.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00408#2#CM085972.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00438#1#CM089177.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG005#1#CM099991.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00609#2#CM086366.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00621#2#CM087869.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00639#1#CM088545.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00658#1#CM088781.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00658#2#CM088803.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00673#2#CM089162.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00733#2#CM099646.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG00735#2#CM087892.1 (chr11: major at k=v1_51 -> minor at k=31)
- HG00738#1#CM086650.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01074#1#CM088791.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01123#1#CM089088.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01150#1#CM088603.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01192#2#CM085946.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01243#2#CM099596.1 (chr11: major at k=v1_51 -> minor at k=31)
- HG01252#2#CM088860.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01258#1#CM088116.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01358#1#CM089060.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01358#2#CM089073.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01433#1#CM086485.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01530#2#CM087563.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01784#1#CM087573.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01891#1#CM088187.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01928#1#CM088993.1 (chr11: minor at k=v1_51 -> major at k=31)
- HG01940#1#CM088636.1 (chr11: minor at k=v1_51 -> major at k=31)
- ... 603 more

### k = published vs k = 31(new): 4551/4551 same call = 100.00%  (haplotypes only in k=published: 0; only in k=31: 0)

| chromosome | same | of | % |
|---|---|---|---|
| chr1 | 139 | 139 | 100.0 |
| chr2 | 206 | 206 | 100.0 |
| chr3 | 178 | 178 | 100.0 |
| chr4 | 166 | 166 | 100.0 |
| chr5 | 358 | 358 | 100.0 |
| chr6 | 230 | 230 | 100.0 |
| chr7 | 268 | 268 | 100.0 |
| chr8 | 342 | 342 | 100.0 |
| chr9 | 246 | 246 | 100.0 |
| chr10 | 285 | 285 | 100.0 |
| chr11 | 394 | 394 | 100.0 |
| chr12 | 381 | 381 | 100.0 |
| chr16 | 255 | 255 | 100.0 |
| chr17 | 159 | 159 | 100.0 |
| chr18 | 168 | 168 | 100.0 |
| chr19 | 328 | 328 | 100.0 |
| chr20 | 306 | 306 | 100.0 |
| chrX | 142 | 142 | 100.0 |

### k = v2_31 vs k = 31(new): 3953/4551 same call = 86.86%  (haplotypes only in k=v2_31: 0; only in k=31: 0)

| chromosome | same | of | % |
|---|---|---|---|
| chr1 | 134 | 139 | 96.4 |
| chr2 | 203 | 206 | 98.5 |
| chr3 | 107 | 178 | 60.1 |
| chr4 | 164 | 166 | 98.8 |
| chr5 | 353 | 358 | 98.6 |
| chr6 | 229 | 230 | 99.6 |
| chr7 | 268 | 268 | 100.0 |
| chr8 | 318 | 342 | 93.0 |
| chr9 | 246 | 246 | 100.0 |
| chr10 | 139 | 285 | 48.8 |
| chr11 | 357 | 394 | 90.6 |
| chr12 | 360 | 381 | 94.5 |
| chr16 | 253 | 255 | 99.2 |
| chr17 | 159 | 159 | 100.0 |
| chr18 | 95 | 168 | 56.5 |
| chr19 | 321 | 328 | 97.9 |
| chr20 | 105 | 306 | 34.3 |
| chrX | 142 | 142 | 100.0 |

Changed calls (598; full list in centromere_changes_v2_31_vs_v1_31.tsv):

- HG01884#2#CM086418.1 (chr1: minor at k=v2_31 -> major at k=31)
- HG01940#2#CM088658.1 (chr1: minor at k=v2_31 -> major at k=31)
- HG02055#1#CM092205.1 (chr1: minor at k=v2_31 -> major at k=31)
- HG03225#1#CM089298.1 (chr1: minor at k=v2_31 -> major at k=31)
- NA18960#1#CM101526.1 (chr1: minor at k=v2_31 -> major at k=31)
- HG00097#1#CM094066.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00097#2#CM094081.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00099#1#CM087321.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00099#2#CM087365.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00133#2#CM090068.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00253#1#CM094463.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00253#2#CM094485.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00272#1#CM094193.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00280#1#CM087153.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00320#2#CM089969.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00329#1#CM094230.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00344#2#CM094516.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00423#1#CM089404.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00438#1#CM089176.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00609#1#CM086353.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00609#2#CM086365.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00621#1#CM087855.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00639#1#CM088544.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00658#2#CM088802.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00673#2#CM089161.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00706#2#CM100064.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00735#2#CM087891.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00738#1#CM086649.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG00741#2#CM087915.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01071#2#CM089134.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01106#1#CM087956.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01106#2#CM087968.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01109#1#CM092091.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01175#1#CM087929.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01192#2#CM085945.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01243#2#CM099595.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01252#2#CM088859.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01255#2#CM086724.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01361#2#CM089043.1 (chr10: minor at k=v2_31 -> major at k=31)
- HG01530#1#CM087539.1 (chr10: minor at k=v2_31 -> major at k=31)
- ... 558 more

## k = 21 vs k = 51 directly (the brief's instability flag)

3090/4551 same call = 67.90%; 1461 haplotypes change major/minor between k=21 and k=51.

## FISH-validated haplotypes

| chromosome | sample | note | k=v1_21 | k=v1_31 | k=v1_51 | k=published |
|---|---|---|---|---|---|---|
| chr3 | NA21144 | HSat1A deletion (manuscript: GM21144) | 1:major,2:major | 1:minor,2:minor | 1:minor,2:major | 1:minor,2:minor |
| chr5 | HG00558 | HSat3 loss | 1:minor,2:major | 1:minor,2:major | 1:minor,2:major | 1:minor,2:major |
| chr9 | HG02630 | HSat3 inversion | 1:major,2:major | 1:minor,2:major | 1:major,2:major | 1:minor,2:major |
| chr12 | HG03816 | D12Z3 variant | 1:minor,2:major | 1:minor,2:major | 1:minor,2:major | 1:minor,2:major |
| chr19 | NA18960 | D19Z3 variant (manuscript: GM18960); Major at every A6 setting | 1:major,2:major | 1:major,2:major | 1:major,2:major | 1:major,2:major |

(hap index:call; both haplotypes of the sample on that chromosome are listed when both passed QC)
