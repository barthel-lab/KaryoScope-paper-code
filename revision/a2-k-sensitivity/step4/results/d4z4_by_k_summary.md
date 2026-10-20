# D4Z4 array reconstruction at k = 21, 31, 51 -- summary

Caller: barthel-lab/HPRC-D4Z4 analyze_d4z4_arrays.py @ e09c1d7, unchanged, published constants. Inputs: the same 582 q-terminal 500 kb windows (249 chr4 / 333 chr10 T2T contigs, HPRC release 2), region set re-annotated with the v1-like rebuilds HKS_human_CHM13_v1like_region_k{21,31,51} (v1 region BED regenerated from the 2025-08 censat input with the 2025-08-25 generator; leaves and priorities identical to the installed v1), gene/subtelomeric/repeat from the v2 six-k annotation (those sets match the published inputs to <0.001%), rebuilt into the published input representation (100 bp bins, contig-frame boundaries). 'published' = the caller on the published v1.1.0 inputs.
CONTROLS: (a) the installed v1 KMC database on the same windows reproduces the published caller output on all 582 haplotypes in every column; (b) the v1-like HKS k=31 build reproduces (a) bit-for-bit (label bp identical, caller output identical). Rows k=v2_21/31/51 are the six-k (v2) family, kept as the database-generation contrast: v2 labels the whole D4Z4 array as one bSat interval, so the caller cannot resolve units there.

## Totals per k

| k | haplotypes | repeat units (published quantity: sum of repeat_count_total) | main-array repeat units | arrays | canonical | degraded |
|---|---|---|---|---|---|---|
| 21 | 582 | 13620 | 10193 | 624 | 616 | 8 |
| 31 | 582 | 28580 | 14433 | 629 | 621 | 8 |
| 51 | 582 | 33332 | 0 | 0 | 0 | 0 |
| published | 582 | 28580 | 14433 | 629 | 621 | 8 |
| v2_21 | 582 | 1286 | 54 | 15 | 6 | 9 |
| v2_31 | 582 | 1378 | 100 | 20 | 7 | 13 |
| v2_51 | 582 | 1430 | 107 | 23 | 5 | 18 |

Per chromosome:

| chromosome | k | haplotypes | repeat units | main-array units | arrays |
|---|---|---|---|---|---|
| chr4 | 21 | 249 | 8126 | 6675 | 257 |
| chr4 | 31 | 249 | 13736 | 7007 | 257 |
| chr4 | 51 | 249 | 15620 | 0 | 0 |
| chr4 | published | 249 | 13736 | 7007 | 257 |
| chr10 | 21 | 333 | 5494 | 3518 | 367 |
| chr10 | 31 | 333 | 14844 | 7426 | 372 |
| chr10 | 51 | 333 | 17712 | 0 | 0 |
| chr10 | published | 333 | 14844 | 7426 | 372 |

## Per-haplotype count concordance

| comparison | count | n | Pearson r | median abs diff | n changed by >1 unit | n changed at all | mean diff (a-b) |
|---|---|---|---|---|---|---|---|
| k=21 vs k=31 | total | 582 | 0.85747 | 23.0 | 578 | 578 | -25.704 |
| k=21 vs k=31 | main array | 582 | 0.81264 | 5.0 | 358 | 388 | -7.285 |
| k=51 vs k=31 | total | 582 | 0.97606 | 7.0 | 538 | 565 | +8.165 |
| k=51 vs k=31 | main array | 582 | nan | 22.0 | 577 | 577 | -24.799 |
| k=31 vs k=published | total | 582 | 1.00000 | 0.0 | 0 | 0 | +0.000 |
| k=31 vs k=published | main array | 582 | 1.00000 | 0.0 | 0 | 0 | +0.000 |
| k=21 vs k=51 | total | 582 | 0.84367 | 29.0 | 579 | 579 | -33.869 |
| k=21 vs k=51 | main array | 582 | nan | 14.0 | 576 | 577 | +17.514 |

## Configuration classes per k

| class | k=21 | k=31 | k=51 | k=published |
|---|---|---|---|---|
| single A (beta-sat terminal) | 424 | 424 | 0 | 424 |
| single B (no beta-sat terminal) | 125 | 125 | 0 | 125 |
| multi-array | 24 | 25 | 0 | 25 |
| degraded / no terminal | 4 | 3 | 0 | 3 |
| no array | 5 | 5 | 582 | 5 |

## Haplotypes whose configuration class changes

### k=21 vs k=31: 5 haplotype(s)

- chr10_hap1_HG01361#1#CM089027.1: single A (beta-sat terminal) (units 11) at k=21 -> degraded / no terminal (units 21) at k=31
- chr10_hap2_HG01928#2#CM089008.1: single A (beta-sat terminal) (units 14) at k=21 -> degraded / no terminal (units 32) at k=31
- chr10_hap1_HG03471#1#CM087526.1: degraded / no terminal (units 16) at k=21 -> multi-array (units 46) at k=31
- chr10_hap1_HG03540#1#CM088053.1: degraded / no terminal (units 31) at k=21 -> single A (beta-sat terminal) (units 91) at k=31
- chr4_hap1_HG03139#1#CM091853.1: degraded / no terminal (units 51) at k=21 -> single A (beta-sat terminal) (units 95) at k=31

### k=51 vs k=31: 577 haplotype(s)

- chr10_hap1_HG00097#1#CM094066.1: no array (units 80) at k=51 -> single A (beta-sat terminal) (units 69) at k=31
- chr10_hap2_HG00097#2#CM094081.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr10_hap1_HG00099#1#CM087321.1: no array (units 28) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG00099#2#CM087365.1: no array (units 129) at k=51 -> single A (beta-sat terminal) (units 114) at k=31
- chr10_hap2_HG00126#2#CM090128.1: no array (units 40) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr10_hap1_HG00128#1#CM090084.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap2_HG00128#2#CM090099.1: no array (units 47) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_HG00133#2#CM090068.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap1_HG00140#1#CM087115.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr10_hap2_HG00140#2#CM087140.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap2_HG00146#2#CM090038.1: no array (units 100) at k=51 -> multi-array (units 88) at k=31
- chr10_hap1_HG00253#1#CM094463.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr10_hap2_HG00253#2#CM094485.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr10_hap1_HG00272#1#CM094193.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 53) at k=31
- chr10_hap2_HG00272#2#CM094220.1: no array (units 34) at k=51 -> single B (no beta-sat terminal) (units 30) at k=31
- chr10_hap1_HG00280#1#CM087153.1: no array (units 96) at k=51 -> single A (beta-sat terminal) (units 76) at k=31
- chr10_hap2_HG00280#2#CM087164.1: no array (units 70) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr10_hap2_HG00320#2#CM089969.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 55) at k=31
- chr10_hap1_HG00321#1#CM089545.1: no array (units 85) at k=51 -> single A (beta-sat terminal) (units 69) at k=31
- chr10_hap2_HG00321#2#CM089576.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap1_HG00329#1#CM094230.1: no array (units 46) at k=51 -> single A (beta-sat terminal) (units 39) at k=31
- chr10_hap2_HG00344#2#CM094516.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 53) at k=31
- chr10_hap1_HG00350#1#CM094475.1: no array (units 22) at k=51 -> single A (beta-sat terminal) (units 16) at k=31
- chr10_hap2_HG00350#2#CM094494.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 55) at k=31
- chr10_hap1_HG00423#1#CM089404.1: no array (units 85) at k=51 -> single A (beta-sat terminal) (units 74) at k=31
- chr10_hap2_HG00423#2#CM089419.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_HG00438#1#CM089176.1: no array (units 28) at k=51 -> degraded / no terminal (units 30) at k=31
- chr10_hap2_HG00438#2#CM089194.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr10_hap1_HG005#1#CM099990.1: no array (units 114) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr10_hap2_HG005#2#CM099979.1: no array (units 23) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr10_hap2_HG00558#2#CM088528.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr10_hap1_HG00597#1#CM085773.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 25) at k=31
- chr10_hap2_HG00597#2#CM085787.1: no array (units 81) at k=51 -> single A (beta-sat terminal) (units 73) at k=31
- chr10_hap1_HG00609#1#CM086353.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 31) at k=31
- chr10_hap2_HG00609#2#CM086365.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap1_HG00621#1#CM087855.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 26) at k=31
- chr10_hap2_HG00621#2#CM087868.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 41) at k=31
- chr10_hap1_HG00639#1#CM088544.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_HG00639#2#CM088558.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_HG00642#2#CM086624.1: no array (units 82) at k=51 -> single B (no beta-sat terminal) (units 81) at k=31
- chr10_hap1_HG00658#1#CM088780.1: no array (units 57) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr10_hap2_HG00658#2#CM088802.1: no array (units 40) at k=51 -> single A (beta-sat terminal) (units 35) at k=31
- chr10_hap1_HG00673#1#CM089148.1: no array (units 34) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr10_hap2_HG00673#2#CM089161.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 31) at k=31
- chr10_hap2_HG00706#2#CM100064.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap2_HG00733#2#CM099645.1: no array (units 32) at k=51 -> single A (beta-sat terminal) (units 28) at k=31
- chr10_hap1_HG00735#1#CM087882.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr10_hap2_HG00735#2#CM087891.1: no array (units 21) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap1_HG00738#1#CM086649.1: no array (units 70) at k=51 -> single A (beta-sat terminal) (units 58) at k=31
- chr10_hap1_HG00741#1#CM087902.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG00741#2#CM087915.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap1_HG01071#1#CM089119.1: no array (units 85) at k=51 -> single A (beta-sat terminal) (units 68) at k=31
- chr10_hap2_HG01071#2#CM089134.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 19) at k=31
- chr10_hap1_HG01074#1#CM088790.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG01074#2#CM088816.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 41) at k=31
- chr10_hap1_HG01081#1#CM088571.1: no array (units 24) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap2_HG01081#2#CM088587.1: no array (units 55) at k=51 -> single A (beta-sat terminal) (units 45) at k=31
- chr10_hap1_HG01099#1#CM086263.1: no array (units 86) at k=51 -> single A (beta-sat terminal) (units 72) at k=31
- chr10_hap1_HG01106#1#CM087956.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap2_HG01106#2#CM087968.1: no array (units 44) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_HG01109#1#CM092091.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr10_hap2_HG01109#2#CM092079.1: no array (units 27) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap1_HG01123#1#CM089087.1: no array (units 27) at k=51 -> single A (beta-sat terminal) (units 23) at k=31
- chr10_hap1_HG01150#1#CM088602.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG01150#2#CM088618.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr10_hap2_HG01167#2#CM094546.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr10_hap1_HG01175#1#CM087929.1: no array (units 47) at k=51 -> single B (no beta-sat terminal) (units 45) at k=31
- chr10_hap2_HG01175#2#CM087943.1: no array (units 28) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap2_HG01192#2#CM085945.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_HG01243#1#CM099607.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap2_HG01243#2#CM099595.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 31) at k=31
- chr10_hap1_HG01252#1#CM088828.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr10_hap2_HG01252#2#CM088859.1: no array (units 22) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap1_HG01255#1#CM086711.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap2_HG01255#2#CM086724.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap2_HG01261#2#CM085921.1: no array (units 69) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr10_hap1_HG01358#1#CM089059.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG01358#2#CM089072.1: no array (units 19) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap1_HG01361#1#CM089027.1: no array (units 24) at k=51 -> degraded / no terminal (units 21) at k=31
- chr10_hap2_HG01361#2#CM089043.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 23) at k=31
- chr10_hap2_HG01433#2#CM086515.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG01496#2#CM086392.1: no array (units 44) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap1_HG01530#1#CM087539.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr10_hap2_HG01530#2#CM087562.1: no array (units 76) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr10_hap2_HG01784#2#CM087601.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr10_hap1_HG01786#1#CM089508.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap2_HG01786#2#CM089532.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr10_hap1_HG01884#1#CM086409.1: no array (units 34) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap1_HG01928#1#CM088992.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr10_hap2_HG01928#2#CM089008.1: no array (units 32) at k=51 -> degraded / no terminal (units 32) at k=31
- chr10_hap1_HG01934#1#CM088844.1: no array (units 23) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr10_hap2_HG01934#2#CM088885.1: no array (units 18) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap2_HG01940#2#CM088665.1: no array (units 39) at k=51 -> single B (no beta-sat terminal) (units 32) at k=31
- chr10_hap1_HG01943#1#CM089355.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr10_hap2_HG01943#2#CM089364.1: no array (units 41) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap1_HG01952#1#CM087978.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap2_HG01952#2#CM088003.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap1_HG01960#1#CM088649.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap1_HG01969#1#CM088691.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr10_hap2_HG01969#2#CM088703.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap1_HG01975#1#CM099801.1: no array (units 23) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr10_hap2_HG01975#2#CM099816.1: no array (units 78) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr10_hap1_HG01978#1#CM089250.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr10_hap2_HG01978#2#CM089276.1: no array (units 55) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr10_hap1_HG01981#1#CM086325.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 31) at k=31
- chr10_hap1_HG01993#1#CM086498.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap2_HG01993#2#CM086538.1: no array (units 46) at k=51 -> single A (beta-sat terminal) (units 34) at k=31
- chr10_hap1_HG02004#1#CM086305.1: no array (units 58) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr10_hap2_HG02004#2#CM086314.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr10_hap1_HG02027#1#CM086979.1: no array (units 41) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap2_HG02027#2#CM086991.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 29) at k=31
- chr10_hap2_HG02040#2#CM087408.1: no array (units 17) at k=51 -> single A (beta-sat terminal) (units 17) at k=31
- chr10_hap1_HG02055#1#CM092210.1: no array (units 77) at k=51 -> single A (beta-sat terminal) (units 62) at k=31
- chr10_hap2_HG02055#2#CM092198.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr10_hap1_HG02056#1#CM085874.1: no array (units 31) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr10_hap2_HG02056#2#CM085886.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap1_HG02071#1#CM088264.1: no array (units 24) at k=51 -> single A (beta-sat terminal) (units 17) at k=31
- chr10_hap2_HG02071#2#CM088287.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap1_HG02080#1#CM092065.1: no array (units 46) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap2_HG02080#2#CM092051.1: no array (units 61) at k=51 -> single A (beta-sat terminal) (units 55) at k=31
- chr10_hap1_HG02109#1#CM092165.1: no array (units 76) at k=51 -> single A (beta-sat terminal) (units 64) at k=31
- chr10_hap1_HG02129#1#CM085853.1: no array (units 121) at k=51 -> single A (beta-sat terminal) (units 102) at k=31
- chr10_hap2_HG02132#2#CM086911.1: no array (units 31) at k=51 -> single A (beta-sat terminal) (units 25) at k=31
- chr10_hap1_HG02135#1#CM088312.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr10_hap1_HG02145#1#CM099676.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 74) at k=31
- chr10_hap2_HG02145#2#CM099661.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap1_HG02148#1#CM087987.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr10_hap2_HG02148#2#CM088028.1: no array (units 68) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr10_hap1_HG02155#1#CM091915.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr10_hap2_HG02155#2#CM091927.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 25) at k=31
- chr10_hap1_HG02165#1#CM087447.1: no array (units 24) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap2_HG02165#2#CM087462.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr10_hap1_HG02178#1#CM089927.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 57) at k=31
- chr10_hap2_HG02178#2#CM089940.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr10_hap1_HG02257#1#CM089225.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap2_HG02257#2#CM089264.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr10_hap1_HG02258#1#CM085828.1: no array (units 39) at k=51 -> multi-array (units 32) at k=31
- chr10_hap2_HG02258#2#CM085840.1: no array (units 65) at k=51 -> single A (beta-sat terminal) (units 53) at k=31
- chr10_hap1_HG02273#1#CM087024.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG02273#2#CM087036.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap1_HG02280#1#CM086293.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap1_HG02293#1#CM086251.1: no array (units 67) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr10_hap2_HG02293#2#CM086240.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_HG02300#1#CM086438.1: no array (units 24) at k=51 -> single A (beta-sat terminal) (units 19) at k=31
- chr10_hap2_HG02300#2#CM086449.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap2_HG02391#2#CM087617.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap1_HG02392#1#CM099884.1: no array (units 50) at k=51 -> single B (no beta-sat terminal) (units 29) at k=31
- chr10_hap2_HG02451#2#CM086229.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr10_hap2_HG02486#2#CM088223.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 45) at k=31
- chr10_hap1_HG02514#1#CM098735.1: no array (units 40) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap2_HG02514#2#CM098751.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG02523#2#CM086942.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 53) at k=31
- chr10_hap1_HG02559#1#CM088237.1: no array (units 68) at k=51 -> single A (beta-sat terminal) (units 57) at k=31
- chr10_hap1_HG02572#1#CM088159.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap2_HG02572#2#CM088177.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr10_hap1_HG02602#1#CM099690.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_HG02602#2#CM099704.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr10_hap1_HG02615#1#CM086674.1: no array (units 32) at k=51 -> single A (beta-sat terminal) (units 23) at k=31
- chr10_hap1_HG02630#1#CM091812.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr10_hap2_HG02630#2#CM091744.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 26) at k=31
- chr10_hap1_HG02647#1#CM086562.1: no array (units 92) at k=51 -> single A (beta-sat terminal) (units 81) at k=31
- chr10_hap2_HG02647#2#CM086576.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr10_hap2_HG02668#2#CM088900.1: no array (units 70) at k=51 -> multi-array (units 59) at k=31
- chr10_hap1_HG02698#1#CM086829.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 26) at k=31
- chr10_hap1_HG02717#1#CM088962.1: no array (units 164) at k=51 -> single A (beta-sat terminal) (units 145) at k=31
- chr10_hap2_HG02717#2#CM088975.1: no array (units 104) at k=51 -> multi-array (units 82) at k=31
- chr10_hap1_HG02723#1#CM092035.1: no array (units 40) at k=51 -> multi-array (units 34) at k=31
- chr10_hap2_HG02723#2#CM092016.1: no array (units 76) at k=51 -> single A (beta-sat terminal) (units 58) at k=31
- chr10_hap1_HG02735#1#CM089334.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 28) at k=31
- chr10_hap2_HG02735#2#CM089347.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_HG02738#1#CM099923.1: no array (units 82) at k=51 -> single A (beta-sat terminal) (units 69) at k=31
- chr10_hap2_HG02738#2#CM099934.1: no array (units 46) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap1_HG02809#1#CM088348.1: no array (units 111) at k=51 -> single A (beta-sat terminal) (units 93) at k=31
- chr10_hap2_HG02809#2#CM088367.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_HG02818#2#CM091991.1: no array (units 156) at k=51 -> multi-array (units 135) at k=31
- chr10_hap2_HG02841#2#CM098675.1: no array (units 93) at k=51 -> single A (beta-sat terminal) (units 85) at k=31
- chr10_hap2_HG02886#2#CM088131.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 26) at k=31
- chr10_hap1_HG02922#1#CM087002.1: no array (units 49) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr10_hap1_HG02965#1#CM087350.1: no array (units 96) at k=51 -> single A (beta-sat terminal) (units 81) at k=31
- chr10_hap2_HG02976#2#CM087126.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap1_HG02984#1#CM098794.1: no array (units 57) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr10_hap2_HG02984#2#CM098808.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap1_HG03017#1#CM088915.1: no array (units 22) at k=51 -> single A (beta-sat terminal) (units 19) at k=31
- chr10_hap2_HG03017#2#CM088932.1: no array (units 97) at k=51 -> single A (beta-sat terminal) (units 78) at k=31
- chr10_hap1_HG03041#1#CM088718.1: no array (units 41) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_HG03050#1#CM098767.1: no array (units 44) at k=51 -> single A (beta-sat terminal) (units 41) at k=31
- chr10_hap2_HG03050#2#CM098781.1: no array (units 140) at k=51 -> single A (beta-sat terminal) (units 117) at k=31
- chr10_hap1_HG03098#1#CM092235.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap1_HG03130#1#CM087259.1: no array (units 63) at k=51 -> single A (beta-sat terminal) (units 54) at k=31
- chr10_hap2_HG03130#2#CM087273.1: no array (units 98) at k=51 -> single A (beta-sat terminal) (units 82) at k=31
- chr10_hap2_HG03139#2#CM091841.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 63) at k=31
- chr10_hap1_HG03209#1#CM087217.1: no array (units 57) at k=51 -> multi-array (units 51) at k=31
- chr10_hap2_HG03225#2#CM089321.1: no array (units 24) at k=51 -> single A (beta-sat terminal) (units 17) at k=31
- chr10_hap1_HG03239#1#CM088751.1: no array (units 78) at k=51 -> single A (beta-sat terminal) (units 64) at k=31
- chr10_hap2_HG03239#2#CM088765.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap1_HG03270#1#CM089904.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 25) at k=31
- chr10_hap2_HG03270#2#CM089916.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap1_HG03369#1#CM094387.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 63) at k=31
- chr10_hap2_HG03369#2#CM094363.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 58) at k=31
- chr10_hap1_HG03453#1#CM088015.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr10_hap2_HG03453#2#CM088039.1: no array (units 117) at k=51 -> single A (beta-sat terminal) (units 102) at k=31
- chr10_hap1_HG03470#1#CM091826.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 66) at k=31
- chr10_hap2_HG03470#2#CM091874.1: no array (units 89) at k=51 -> multi-array (units 76) at k=31
- chr10_hap1_HG03471#1#CM087526.1: no array (units 53) at k=51 -> multi-array (units 46) at k=31
- chr10_hap2_HG03471#2#CM087549.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr10_hap1_HG03486#1#CM091978.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr10_hap2_HG03486#2#CM091964.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr10_hap1_HG03516#1#CM089209.1: no array (units 69) at k=51 -> single A (beta-sat terminal) (units 62) at k=31
- chr10_hap1_HG03521#1#CM094417.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr10_hap1_HG03540#1#CM088053.1: no array (units 102) at k=51 -> single A (beta-sat terminal) (units 91) at k=31
- chr10_hap1_HG03579#1#CM088949.1: no array (units 57) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr10_hap2_HG03579#2#CM089292.1: no array (units 87) at k=51 -> single A (beta-sat terminal) (units 78) at k=31
- chr10_hap1_HG03654#1#CM086839.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 28) at k=31
- chr10_hap2_HG03654#2#CM086859.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 27) at k=31
- chr10_hap1_HG03669#1#CM086736.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr10_hap2_HG03669#2#CM086747.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap1_HG03688#1#CM086806.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr10_hap1_HG03704#1#CM086215.1: no array (units 55) at k=51 -> single B (no beta-sat terminal) (units 46) at k=31
- chr10_hap2_HG03704#2#CM086193.1: no array (units 58) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr10_hap1_HG03742#1#CM087632.1: no array (units 31) at k=51 -> single A (beta-sat terminal) (units 25) at k=31
- chr10_hap2_HG03742#2#CM087647.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr10_hap1_HG03784#1#CM094178.1: no array (units 41) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_HG03784#2#CM094205.1: no array (units 88) at k=51 -> single A (beta-sat terminal) (units 70) at k=31
- chr10_hap2_HG03804#2#CM088383.1: no array (units 65) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr10_hap1_HG03816#1#CM100031.1: no array (units 98) at k=51 -> single A (beta-sat terminal) (units 86) at k=31
- chr10_hap2_HG03816#2#CM100046.1: no array (units 48) at k=51 -> single B (no beta-sat terminal) (units 50) at k=31
- chr10_hap1_HG03831#1#CM086529.1: no array (units 67) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr10_hap1_HG03834#1#CM085800.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap1_HG03874#1#CM089436.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap2_HG03874#2#CM089478.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 28) at k=31
- chr10_hap1_HG03927#1#CM086638.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr10_hap2_HG03927#2#CM086663.1: no array (units 21) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr10_hap2_HG03942#2#CM088407.1: no array (units 35) at k=51 -> single B (no beta-sat terminal) (units 27) at k=31
- chr10_hap1_HG04115#1#CM086757.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr10_hap1_HG04157#1#CM088417.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 55) at k=31
- chr10_hap1_HG04160#1#CM088430.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap2_HG04160#2#CM088451.1: no array (units 39) at k=51 -> single B (no beta-sat terminal) (units 29) at k=31
- chr10_hap1_HG04184#1#CM086463.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap1_HG04187#1#CM099746.1: no array (units 86) at k=51 -> single A (beta-sat terminal) (units 72) at k=31
- chr10_hap2_HG04187#2#CM099760.1: no array (units 72) at k=51 -> single A (beta-sat terminal) (units 62) at k=31
- chr10_hap2_HG04199#2#CM086793.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr10_hap2_HG04204#2#CM088472.1: no array (units 78) at k=51 -> single A (beta-sat terminal) (units 65) at k=31
- chr10_hap1_HG04228#1#CM088484.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 29) at k=31
- chr10_hap2_HG04228#2#CM088512.1: no array (units 67) at k=51 -> single B (no beta-sat terminal) (units 69) at k=31
- chr10_hap1_HG06807#1#CM101404.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap2_HG06807#2#CM101438.1: no array (units 50) at k=51 -> multi-array (units 40) at k=31
- chr10_hap2_NA18522#2#CM087194.1: no array (units 57) at k=51 -> multi-array (units 49) at k=31
- chr10_hap1_NA18565#1#CM094322.1: no array (units 100) at k=51 -> multi-array (units 82) at k=31
- chr10_hap2_NA18565#2#CM094344.1: no array (units 31) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap1_NA18570#1#CM087663.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 29) at k=31
- chr10_hap2_NA18570#2#CM087681.1: no array (units 42) at k=51 -> single B (no beta-sat terminal) (units 36) at k=31
- chr10_hap1_NA18608#1#CM089843.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 39) at k=31
- chr10_hap1_NA18620#1#CM089823.1: no array (units 55) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap2_NA18620#2#CM089833.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr10_hap2_NA18747#2#CM087245.1: no array (units 18) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap1_NA18879#1#CM099718.1: no array (units 59) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr10_hap2_NA18879#2#CM099733.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr10_hap1_NA18906#1#CM091938.1: no array (units 123) at k=51 -> single A (beta-sat terminal) (units 110) at k=31
- chr10_hap1_NA18943#1#CM101663.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 35) at k=31
- chr10_hap2_NA18943#2#CM101692.1: no array (units 40) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap1_NA18944#1#CM101649.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr10_hap1_NA18945#1#CM101617.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr10_hap2_NA18945#2#CM101632.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 26) at k=31
- chr10_hap1_NA18948#1#CM101603.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap1_NA18952#1#CM089788.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap2_NA18959#2#CM101577.1: no array (units 20) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr10_hap2_NA18967#2#CM101516.1: no array (units 41) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap1_NA18970#1#CM101471.1: no array (units 93) at k=51 -> single B (no beta-sat terminal) (units 68) at k=31
- chr10_hap2_NA18970#2#CM101483.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap1_NA18971#1#CM087296.1: no array (units 26) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr10_hap2_NA18971#2#CM087284.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr10_hap1_NA18974#1#CM089776.1: no array (units 44) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr10_hap2_NA18974#2#CM089800.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 29) at k=31
- chr10_hap1_NA18976#1#CM089755.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr10_hap2_NA18976#2#CM089764.1: no array (units 65) at k=51 -> single A (beta-sat terminal) (units 53) at k=31
- chr10_hap1_NA18982#1#CM101423.1: no array (units 22) at k=51 -> single A (beta-sat terminal) (units 17) at k=31
- chr10_hap2_NA18982#2#CM101456.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 45) at k=31
- chr10_hap1_NA18983#1#CM087423.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_NA18983#2#CM087435.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap1_NA19036#1#CM089728.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr10_hap1_NA19043#1#CM087049.1: no array (units 100) at k=51 -> single A (beta-sat terminal) (units 87) at k=31
- chr10_hap2_NA19043#2#CM087075.1: no array (units 92) at k=51 -> single A (beta-sat terminal) (units 92) at k=31
- chr10_hap1_NA19087#1#CM087700.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr10_hap2_NA19087#2#CM087711.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 34) at k=31
- chr10_hap1_NA19159#1#CM099947.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 35) at k=31
- chr10_hap2_NA19159#2#CM099966.1: no array (units 55) at k=51 -> single A (beta-sat terminal) (units 54) at k=31
- chr10_hap1_NA19185#1#CM087728.1: no array (units 118) at k=51 -> multi-array (units 95) at k=31
- chr10_hap2_NA19185#2#CM087744.1: no array (units 59) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap1_NA19240#1#CM099581.1: no array (units 89) at k=51 -> multi-array (units 74) at k=31
- chr10_hap2_NA19240#2#CM099619.1: no array (units 83) at k=51 -> single A (beta-sat terminal) (units 73) at k=31
- chr10_hap1_NA19338#1#CM087762.1: no array (units 55) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr10_hap2_NA19338#2#CM087779.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 39) at k=31
- chr10_hap1_NA19443#1#CM089450.1: no array (units 130) at k=51 -> single A (beta-sat terminal) (units 111) at k=31
- chr10_hap1_NA19468#1#CM087822.1: no array (units 84) at k=51 -> multi-array (units 65) at k=31
- chr10_hap2_NA19468#2#CM087837.1: no array (units 63) at k=51 -> multi-array (units 48) at k=31
- chr10_hap1_NA19700#1#CM089625.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr10_hap2_NA19700#2#CM089642.1: no array (units 98) at k=51 -> single A (beta-sat terminal) (units 81) at k=31
- chr10_hap2_NA19776#2#CM094053.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap2_NA19835#2#CM094038.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 64) at k=31
- chr10_hap1_NA19909#1#CM094123.1: no array (units 24) at k=51 -> single A (beta-sat terminal) (units 19) at k=31
- chr10_hap2_NA19909#2#CM094136.1: no array (units 121) at k=51 -> single A (beta-sat terminal) (units 108) at k=31
- chr10_hap1_NA20129#1#CM092136.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap2_NA20282#2#CM094283.1: no array (units 49) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr10_hap1_NA20346#1#CM094253.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_NA20346#2#CM094265.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 35) at k=31
- chr10_hap1_NA20503#1#CM102463.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 41) at k=31
- chr10_hap2_NA20503#2#CM102480.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr10_hap2_NA20752#2#CM099786.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr10_hap1_NA20762#1#CM102516.1: no array (units 70) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr10_hap2_NA20799#2#CM098721.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr10_hap1_NA20805#1#CM091890.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr10_hap2_NA20805#2#CM091902.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 29) at k=31
- chr10_hap1_NA20806#1#CM102446.1: no array (units 67) at k=51 -> single A (beta-sat terminal) (units 57) at k=31
- chr10_hap1_NA20827#1#CM102548.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap2_NA20827#2#CM102568.1: no array (units 23) at k=51 -> single A (beta-sat terminal) (units 20) at k=31
- chr10_hap2_NA20850#2#CM094309.1: no array (units 105) at k=51 -> single A (beta-sat terminal) (units 89) at k=31
- chr10_hap1_NA20870#1#CM089701.1: no array (units 41) at k=51 -> single B (no beta-sat terminal) (units 39) at k=31
- chr10_hap2_NA20870#2#CM089715.1: no array (units 18) at k=51 -> single A (beta-sat terminal) (units 15) at k=31
- chr10_hap1_NA20905#1#CM087477.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 27) at k=31
- chr10_hap2_NA20905#2#CM087511.1: no array (units 30) at k=51 -> single A (beta-sat terminal) (units 25) at k=31
- chr10_hap1_NA21093#1#CM089601.1: no array (units 141) at k=51 -> single A (beta-sat terminal) (units 118) at k=31
- chr10_hap2_NA21093#2#CM089611.1: no array (units 141) at k=51 -> single A (beta-sat terminal) (units 113) at k=31
- chr10_hap1_NA21106#1#CM089560.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr10_hap2_NA21106#2#CM089591.1: no array (units 35) at k=51 -> single B (no beta-sat terminal) (units 33) at k=31
- chr10_hap1_NA21110#1#CM089669.1: no array (units 80) at k=51 -> single A (beta-sat terminal) (units 71) at k=31
- chr10_hap2_NA21110#2#CM089685.1: no array (units 55) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr10_hap1_NA21144#1#CM094098.1: no array (units 44) at k=51 -> single A (beta-sat terminal) (units 36) at k=31
- chr10_hap2_NA21144#2#CM094114.1: no array (units 84) at k=51 -> single A (beta-sat terminal) (units 67) at k=31
- chr10_hap1_NA21309#1#CM092113.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr4_hap2_HG00099#2#CM087359.1: no array (units 56) at k=51 -> single B (no beta-sat terminal) (units 49) at k=31
- chr4_hap1_HG00128#1#CM090078.1: no array (units 72) at k=51 -> single A (beta-sat terminal) (units 57) at k=31
- chr4_hap2_HG00128#2#CM090094.1: no array (units 35) at k=51 -> single B (no beta-sat terminal) (units 28) at k=31
- chr4_hap1_HG00133#1#CM090048.1: no array (units 79) at k=51 -> single A (beta-sat terminal) (units 73) at k=31
- chr4_hap2_HG00133#2#CM090063.1: no array (units 73) at k=51 -> multi-array (units 67) at k=31
- chr4_hap1_HG00140#1#CM087110.1: no array (units 84) at k=51 -> single B (no beta-sat terminal) (units 67) at k=31
- chr4_hap2_HG00140#2#CM087135.1: no array (units 82) at k=51 -> single B (no beta-sat terminal) (units 73) at k=31
- chr4_hap2_HG00146#2#CM090035.1: no array (units 71) at k=51 -> single B (no beta-sat terminal) (units 61) at k=31
- chr4_hap1_HG00272#1#CM094187.1: no array (units 102) at k=51 -> single B (no beta-sat terminal) (units 98) at k=31
- chr4_hap2_HG00272#2#CM094214.1: no array (units 49) at k=51 -> single B (no beta-sat terminal) (units 40) at k=31
- chr4_hap1_HG00280#1#CM087148.1: no array (units 58) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr4_hap2_HG00280#2#CM087159.1: no array (units 274) at k=51 -> single A (beta-sat terminal) (units 227) at k=31
- chr4_hap1_HG00320#1#CM089949.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 45) at k=31
- chr4_hap2_HG00320#2#CM089963.1: no array (units 52) at k=51 -> single B (no beta-sat terminal) (units 42) at k=31
- chr4_hap1_HG00321#1#CM089539.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr4_hap2_HG00321#2#CM089570.1: no array (units 61) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr4_hap2_HG00408#2#CM085967.1: no array (units 44) at k=51 -> single B (no beta-sat terminal) (units 38) at k=31
- chr4_hap1_HG00423#1#CM089399.1: no array (units 29) at k=51 -> single B (no beta-sat terminal) (units 24) at k=31
- chr4_hap2_HG00423#2#CM089415.1: no array (units 29) at k=51 -> single B (no beta-sat terminal) (units 24) at k=31
- chr4_hap1_HG00438#1#CM089170.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 31) at k=31
- chr4_hap2_HG00438#2#CM089189.1: no array (units 66) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr4_hap1_HG00558#1#CM088493.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr4_hap2_HG00558#2#CM088522.1: no array (units 76) at k=51 -> single A (beta-sat terminal) (units 67) at k=31
- chr4_hap2_HG00597#2#CM085782.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr4_hap1_HG00609#1#CM086349.1: no array (units 45) at k=51 -> single B (no beta-sat terminal) (units 36) at k=31
- chr4_hap2_HG00609#2#CM086359.1: no array (units 47) at k=51 -> single B (no beta-sat terminal) (units 41) at k=31
- chr4_hap1_HG00621#1#CM087849.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr4_hap2_HG00621#2#CM087863.1: no array (units 35) at k=51 -> single B (no beta-sat terminal) (units 40) at k=31
- chr4_hap1_HG00639#1#CM088538.1: no array (units 93) at k=51 -> single A (beta-sat terminal) (units 86) at k=31
- chr4_hap2_HG00639#2#CM088553.1: no array (units 61) at k=51 -> single A (beta-sat terminal) (units 62) at k=31
- chr4_hap2_HG00642#2#CM086620.1: no array (units 62) at k=51 -> single B (no beta-sat terminal) (units 39) at k=31
- chr4_hap1_HG00673#1#CM089144.1: no array (units 50) at k=51 -> single B (no beta-sat terminal) (units 50) at k=31
- chr4_hap2_HG00673#2#CM089155.1: no array (units 57) at k=51 -> single B (no beta-sat terminal) (units 54) at k=31
- chr4_hap1_HG00706#1#CM100014.1: no array (units 52) at k=51 -> single B (no beta-sat terminal) (units 46) at k=31
- chr4_hap2_HG00706#2#CM100058.1: no array (units 54) at k=51 -> single B (no beta-sat terminal) (units 46) at k=31
- chr4_hap1_HG00735#1#CM087877.1: no array (units 68) at k=51 -> single A (beta-sat terminal) (units 68) at k=31
- chr4_hap2_HG00735#2#CM087886.1: no array (units 87) at k=51 -> single B (no beta-sat terminal) (units 85) at k=31
- chr4_hap2_HG00738#2#CM086682.1: no array (units 59) at k=51 -> single B (no beta-sat terminal) (units 57) at k=31
- chr4_hap1_HG01071#1#CM089114.1: no array (units 62) at k=51 -> single A (beta-sat terminal) (units 59) at k=31
- chr4_hap2_HG01071#2#CM089130.1: no array (units 86) at k=51 -> single A (beta-sat terminal) (units 76) at k=31
- chr4_hap1_HG01109#1#CM092086.1: no array (units 71) at k=51 -> multi-array (units 57) at k=31
- chr4_hap1_HG01167#1#CM094526.1: no array (units 55) at k=51 -> single B (no beta-sat terminal) (units 37) at k=31
- chr4_hap2_HG01167#2#CM094541.1: no array (units 69) at k=51 -> single B (no beta-sat terminal) (units 54) at k=31
- chr4_hap1_HG01175#1#CM087923.1: no array (units 71) at k=51 -> single B (no beta-sat terminal) (units 68) at k=31
- chr4_hap2_HG01175#2#CM087938.1: no array (units 105) at k=51 -> single A (beta-sat terminal) (units 96) at k=31
- chr4_hap2_HG01192#2#CM085941.1: no array (units 69) at k=51 -> single A (beta-sat terminal) (units 61) at k=31
- chr4_hap1_HG01255#1#CM086705.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr4_hap1_HG01258#1#CM088111.1: no array (units 50) at k=51 -> single B (no beta-sat terminal) (units 47) at k=31
- chr4_hap2_HG01258#2#CM088140.1: no array (units 86) at k=51 -> single B (no beta-sat terminal) (units 67) at k=31
- chr4_hap1_HG01261#1#CM085908.1: no array (units 54) at k=51 -> single B (no beta-sat terminal) (units 46) at k=31
- chr4_hap1_HG01346#1#CM086586.1: no array (units 66) at k=51 -> single B (no beta-sat terminal) (units 73) at k=31
- chr4_hap2_HG01346#2#CM086599.1: no array (units 60) at k=51 -> single B (no beta-sat terminal) (units 58) at k=31
- chr4_hap1_HG01358#1#CM089053.1: no array (units 57) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr4_hap2_HG01358#2#CM089066.1: no array (units 63) at k=51 -> single B (no beta-sat terminal) (units 56) at k=31
- chr4_hap1_HG01361#1#CM089022.1: no array (units 114) at k=51 -> single B (no beta-sat terminal) (units 100) at k=31
- chr4_hap2_HG01361#2#CM089038.1: no array (units 60) at k=51 -> single B (no beta-sat terminal) (units 67) at k=31
- chr4_hap2_HG01433#2#CM086509.1: no array (units 51) at k=51 -> single B (no beta-sat terminal) (units 54) at k=31
- chr4_hap1_HG01530#1#CM087534.1: no array (units 30) at k=51 -> single B (no beta-sat terminal) (units 29) at k=31
- chr4_hap2_HG01530#2#CM087558.1: no array (units 76) at k=51 -> single A (beta-sat terminal) (units 66) at k=31
- chr4_hap1_HG01784#1#CM087568.1: no array (units 35) at k=51 -> single B (no beta-sat terminal) (units 36) at k=31
- chr4_hap1_HG01786#1#CM089502.1: no array (units 40) at k=51 -> single B (no beta-sat terminal) (units 32) at k=31
- chr4_hap2_HG01891#2#CM088207.1: no array (units 74) at k=51 -> multi-array (units 47) at k=31
- chr4_hap1_HG01928#1#CM088986.1: no array (units 60) at k=51 -> single B (no beta-sat terminal) (units 47) at k=31
- chr4_hap1_HG01940#1#CM088631.1: no array (units 52) at k=51 -> single B (no beta-sat terminal) (units 45) at k=31
- chr4_hap2_HG01940#2#CM088661.1: no array (units 44) at k=51 -> single B (no beta-sat terminal) (units 34) at k=31
- chr4_hap2_HG01952#2#CM087997.1: no array (units 66) at k=51 -> single B (no beta-sat terminal) (units 52) at k=31
- chr4_hap1_HG01981#1#CM086321.1: no array (units 50) at k=51 -> single B (no beta-sat terminal) (units 53) at k=31
- chr4_hap2_HG01981#2#CM086335.1: no array (units 107) at k=51 -> single B (no beta-sat terminal) (units 100) at k=31
- chr4_hap1_HG01993#1#CM086494.1: no array (units 64) at k=51 -> single B (no beta-sat terminal) (units 47) at k=31
- chr4_hap2_HG01993#2#CM086535.1: no array (units 64) at k=51 -> single B (no beta-sat terminal) (units 50) at k=31
- chr4_hap2_HG02004#2#CM086311.1: no array (units 95) at k=51 -> single B (no beta-sat terminal) (units 64) at k=31
- chr4_hap1_HG02015#1#CM085891.1: no array (units 38) at k=51 -> single B (no beta-sat terminal) (units 39) at k=31
- chr4_hap2_HG02015#2#CM085898.1: no array (units 78) at k=51 -> single A (beta-sat terminal) (units 68) at k=31
- chr4_hap1_HG02027#1#CM086976.1: no array (units 71) at k=51 -> single A (beta-sat terminal) (units 65) at k=31
- chr4_hap2_HG02027#2#CM086986.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr4_hap1_HG02040#1#CM087374.1: no array (units 43) at k=51 -> single B (no beta-sat terminal) (units 41) at k=31
- chr4_hap1_HG02074#1#CM088269.1: no array (units 95) at k=51 -> single B (no beta-sat terminal) (units 73) at k=31
- chr4_hap2_HG02074#2#CM088295.1: no array (units 39) at k=51 -> single B (no beta-sat terminal) (units 39) at k=31
- chr4_hap2_HG02080#2#CM092046.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 21) at k=31
- chr4_hap1_HG02083#1#CM086953.1: no array (units 41) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr4_hap1_HG02109#1#CM092159.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 52) at k=31
- chr4_hap2_HG02109#2#CM092145.1: no array (units 87) at k=51 -> single A (beta-sat terminal) (units 81) at k=31
- chr4_hap1_HG02129#1#CM085848.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr4_hap2_HG02129#2#CM085859.1: no array (units 37) at k=51 -> single B (no beta-sat terminal) (units 33) at k=31
- chr4_hap1_HG02132#1#CM086895.1: no array (units 40) at k=51 -> single B (no beta-sat terminal) (units 38) at k=31
- chr4_hap2_HG02132#2#CM086906.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr4_hap1_HG02155#1#CM091910.1: no array (units 76) at k=51 -> single A (beta-sat terminal) (units 79) at k=31
- chr4_hap2_HG02155#2#CM091922.1: no array (units 32) at k=51 -> single B (no beta-sat terminal) (units 36) at k=31
- chr4_hap1_HG02165#1#CM087442.1: no array (units 55) at k=51 -> single B (no beta-sat terminal) (units 44) at k=31
- chr4_hap2_HG02165#2#CM087458.1: no array (units 53) at k=51 -> single B (no beta-sat terminal) (units 52) at k=31
- chr4_hap1_HG02178#1#CM089923.1: no array (units 30) at k=51 -> single B (no beta-sat terminal) (units 29) at k=31
- chr4_hap2_HG02178#2#CM089936.1: no array (units 41) at k=51 -> single B (no beta-sat terminal) (units 38) at k=31
- chr4_hap1_HG02257#1#CM089219.1: no array (units 85) at k=51 -> single A (beta-sat terminal) (units 87) at k=31
- chr4_hap2_HG02257#2#CM089258.1: no array (units 75) at k=51 -> multi-array (units 62) at k=31
- chr4_hap1_HG02258#1#CM085824.1: no array (units 106) at k=51 -> single A (beta-sat terminal) (units 97) at k=31
- chr4_hap2_HG02258#2#CM085835.1: no array (units 97) at k=51 -> single A (beta-sat terminal) (units 86) at k=31
- chr4_hap2_HG02280#2#CM086280.1: no array (units 55) at k=51 -> single B (no beta-sat terminal) (units 41) at k=31
- chr4_hap2_HG02391#2#CM087611.1: no array (units 40) at k=51 -> single B (no beta-sat terminal) (units 37) at k=31
- chr4_hap1_HG02392#1#CM099880.1: no array (units 39) at k=51 -> single B (no beta-sat terminal) (units 34) at k=31
- chr4_hap2_HG02451#2#CM086224.1: no array (units 57) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr4_hap2_HG02486#2#CM088217.1: no array (units 59) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr4_hap2_HG02514#2#CM098746.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 38) at k=31
- chr4_hap2_HG02523#2#CM086936.1: no array (units 29) at k=51 -> single B (no beta-sat terminal) (units 30) at k=31
- chr4_hap1_HG02559#1#CM088233.1: no array (units 36) at k=51 -> single A (beta-sat terminal) (units 35) at k=31
- chr4_hap2_HG02559#2#CM088246.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 67) at k=31
- chr4_hap1_HG02572#1#CM088153.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr4_hap2_HG02572#2#CM088171.1: no array (units 131) at k=51 -> single A (beta-sat terminal) (units 118) at k=31
- chr4_hap1_HG02583#1#CM089488.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 62) at k=31
- chr4_hap2_HG02583#2#CM089515.1: no array (units 64) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr4_hap1_HG02602#1#CM099685.1: no array (units 42) at k=51 -> single B (no beta-sat terminal) (units 43) at k=31
- chr4_hap1_HG02615#1#CM086671.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr4_hap2_HG02615#2#CM086692.1: no array (units 127) at k=51 -> single A (beta-sat terminal) (units 110) at k=31
- chr4_hap2_HG02622#2#CM088097.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr4_hap1_HG02630#1#CM091807.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 63) at k=31
- chr4_hap2_HG02630#2#CM091738.1: no array (units 28) at k=51 -> single B (no beta-sat terminal) (units 27) at k=31
- chr4_hap1_HG02647#1#CM086557.1: no array (units 74) at k=51 -> single A (beta-sat terminal) (units 67) at k=31
- chr4_hap2_HG02668#2#CM088895.1: no array (units 68) at k=51 -> single A (beta-sat terminal) (units 68) at k=31
- chr4_hap1_HG02717#1#CM088957.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 31) at k=31
- chr4_hap2_HG02717#2#CM088970.1: no array (units 61) at k=51 -> single A (beta-sat terminal) (units 57) at k=31
- chr4_hap1_HG02723#1#CM092029.1: no array (units 187) at k=51 -> single A (beta-sat terminal) (units 141) at k=31
- chr4_hap2_HG02723#2#CM092011.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 43) at k=31
- chr4_hap1_HG02735#1#CM089329.1: no array (units 68) at k=51 -> single B (no beta-sat terminal) (units 57) at k=31
- chr4_hap2_HG02738#2#CM099931.1: no array (units 64) at k=51 -> single B (no beta-sat terminal) (units 59) at k=31
- chr4_hap1_HG02809#1#CM088343.1: no array (units 38) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr4_hap2_HG02809#2#CM088361.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr4_hap1_HG02841#1#CM098701.1: no array (units 65) at k=51 -> single A (beta-sat terminal) (units 55) at k=31
- chr4_hap1_HG02886#1#CM088084.1: no array (units 59) at k=51 -> multi-array (units 40) at k=31
- chr4_hap2_HG02886#2#CM088125.1: no array (units 39) at k=51 -> single B (no beta-sat terminal) (units 38) at k=31
- chr4_hap2_HG02922#2#CM087011.1: no array (units 105) at k=51 -> single A (beta-sat terminal) (units 80) at k=31
- chr4_hap1_HG02965#1#CM087344.1: no array (units 72) at k=51 -> single A (beta-sat terminal) (units 65) at k=31
- chr4_hap2_HG02965#2#CM087387.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 39) at k=31
- chr4_hap1_HG03017#1#CM088910.1: no array (units 96) at k=51 -> single B (no beta-sat terminal) (units 91) at k=31
- chr4_hap2_HG03017#2#CM088926.1: no array (units 58) at k=51 -> single B (no beta-sat terminal) (units 52) at k=31
- chr4_hap1_HG03041#1#CM088713.1: no array (units 42) at k=51 -> single A (beta-sat terminal) (units 29) at k=31
- chr4_hap2_HG03041#2#CM088730.1: no array (units 18) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr4_hap1_HG03050#1#CM098761.1: no array (units 157) at k=51 -> single A (beta-sat terminal) (units 113) at k=31
- chr4_hap1_HG03098#1#CM092229.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr4_hap2_HG03098#2#CM092217.1: no array (units 100) at k=51 -> single A (beta-sat terminal) (units 85) at k=31
- chr4_hap1_HG03130#1#CM087254.1: no array (units 50) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr4_hap2_HG03130#2#CM087268.1: no array (units 218) at k=51 -> single A (beta-sat terminal) (units 184) at k=31
- chr4_hap1_HG03139#1#CM091853.1: no array (units 93) at k=51 -> single A (beta-sat terminal) (units 95) at k=31
- chr4_hap2_HG03139#2#CM091838.1: no array (units 27) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr4_hap1_HG03209#1#CM087211.1: no array (units 75) at k=51 -> single A (beta-sat terminal) (units 77) at k=31
- chr4_hap2_HG03209#2#CM087229.1: no array (units 40) at k=51 -> single A (beta-sat terminal) (units 30) at k=31
- chr4_hap1_HG03225#1#CM089299.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr4_hap1_HG03239#1#CM088746.1: no array (units 45) at k=51 -> single B (no beta-sat terminal) (units 41) at k=31
- chr4_hap1_HG03453#1#CM088011.1: no array (units 81) at k=51 -> single A (beta-sat terminal) (units 64) at k=31
- chr4_hap1_HG03470#1#CM091820.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 19) at k=31
- chr4_hap2_HG03470#2#CM091868.1: no array (units 77) at k=51 -> single A (beta-sat terminal) (units 69) at k=31
- chr4_hap1_HG03471#1#CM087521.1: no array (units 78) at k=51 -> single A (beta-sat terminal) (units 67) at k=31
- chr4_hap2_HG03471#2#CM087544.1: no array (units 73) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr4_hap1_HG03486#1#CM091972.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 47) at k=31
- chr4_hap2_HG03486#2#CM091958.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr4_hap1_HG03516#1#CM089204.1: no array (units 23) at k=51 -> single B (no beta-sat terminal) (units 17) at k=31
- chr4_hap1_HG03521#1#CM094413.1: no array (units 31) at k=51 -> single A (beta-sat terminal) (units 22) at k=31
- chr4_hap2_HG03521#2#CM094435.1: no array (units 83) at k=51 -> single A (beta-sat terminal) (units 84) at k=31
- chr4_hap1_HG03540#1#CM088048.1: no array (units 29) at k=51 -> single B (no beta-sat terminal) (units 24) at k=31
- chr4_hap2_HG03540#2#CM088061.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 41) at k=31
- chr4_hap1_HG03579#1#CM088943.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 37) at k=31
- chr4_hap2_HG03579#2#CM089286.1: no array (units 105) at k=51 -> single B (no beta-sat terminal) (units 94) at k=31
- chr4_hap2_HG03654#2#CM086855.1: no array (units 48) at k=51 -> single B (no beta-sat terminal) (units 47) at k=31
- chr4_hap1_HG03669#1#CM086732.1: no array (units 35) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr4_hap1_HG03710#1#CM086867.1: no array (units 83) at k=51 -> single A (beta-sat terminal) (units 77) at k=31
- chr4_hap2_HG03710#2#CM086881.1: no array (units 45) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr4_hap1_HG03742#1#CM087626.1: no array (units 101) at k=51 -> single B (no beta-sat terminal) (units 108) at k=31
- chr4_hap2_HG03742#2#CM087643.1: no array (units 47) at k=51 -> single A (beta-sat terminal) (units 40) at k=31
- chr4_hap1_HG03804#1#CM088329.1: no array (units 37) at k=51 -> single A (beta-sat terminal) (units 34) at k=31
- chr4_hap2_HG03804#2#CM088378.1: no array (units 90) at k=51 -> single A (beta-sat terminal) (units 76) at k=31
- chr4_hap1_HG03831#1#CM086525.1: no array (units 58) at k=51 -> single B (no beta-sat terminal) (units 54) at k=31
- chr4_hap2_HG03831#2#CM086546.1: no array (units 48) at k=51 -> single B (no beta-sat terminal) (units 47) at k=31
- chr4_hap1_HG03834#1#CM085795.1: no array (units 56) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr4_hap2_HG03834#2#CM085811.1: no array (units 93) at k=51 -> single A (beta-sat terminal) (units 83) at k=31
- chr4_hap1_HG03874#1#CM089430.1: no array (units 28) at k=51 -> single B (no beta-sat terminal) (units 31) at k=31
- chr4_hap2_HG03874#2#CM089473.1: no array (units 45) at k=51 -> single B (no beta-sat terminal) (units 33) at k=31
- chr4_hap1_HG03927#1#CM086635.1: no array (units 63) at k=51 -> single B (no beta-sat terminal) (units 58) at k=31
- chr4_hap2_HG03927#2#CM086657.1: no array (units 75) at k=51 -> single B (no beta-sat terminal) (units 61) at k=31
- chr4_hap1_HG03942#1#CM088391.1: no array (units 134) at k=51 -> multi-array (units 128) at k=31
- chr4_hap2_HG03942#2#CM088402.1: no array (units 80) at k=51 -> single A (beta-sat terminal) (units 64) at k=31
- chr4_hap1_HG04160#1#CM088425.1: no array (units 128) at k=51 -> single A (beta-sat terminal) (units 109) at k=31
- chr4_hap2_HG04160#2#CM088447.1: no array (units 125) at k=51 -> single A (beta-sat terminal) (units 113) at k=31
- chr4_hap2_HG04187#2#CM099754.1: no array (units 61) at k=51 -> single A (beta-sat terminal) (units 54) at k=31
- chr4_hap1_HG04199#1#CM086778.1: no array (units 60) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr4_hap1_HG06807#1#CM101398.1: no array (units 75) at k=51 -> single A (beta-sat terminal) (units 74) at k=31
- chr4_hap2_HG06807#2#CM101432.1: no array (units 64) at k=51 -> multi-array (units 45) at k=31
- chr4_hap1_NA18505#1#CM089859.1: no array (units 39) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr4_hap2_NA18505#2#CM089869.1: no array (units 26) at k=51 -> single B (no beta-sat terminal) (units 19) at k=31
- chr4_hap1_NA18508#1#CM089648.1: no array (units 23) at k=51 -> single B (no beta-sat terminal) (units 20) at k=31
- chr4_hap1_NA18522#1#CM087174.1: no array (units 148) at k=51 -> single A (beta-sat terminal) (units 125) at k=31
- chr4_hap2_NA18522#2#CM087189.1: no array (units 63) at k=51 -> single A (beta-sat terminal) (units 60) at k=31
- chr4_hap1_NA18565#1#CM094318.1: no array (units 53) at k=51 -> single B (no beta-sat terminal) (units 54) at k=31
- chr4_hap2_NA18565#2#CM094340.1: no array (units 46) at k=51 -> single B (no beta-sat terminal) (units 40) at k=31
- chr4_hap1_NA18570#1#CM087657.1: no array (units 66) at k=51 -> multi-array (units 61) at k=31
- chr4_hap2_NA18570#2#CM087676.1: no array (units 92) at k=51 -> single A (beta-sat terminal) (units 74) at k=31
- chr4_hap2_NA18608#2#CM089851.1: no array (units 45) at k=51 -> single B (no beta-sat terminal) (units 37) at k=31
- chr4_hap1_NA18879#1#CM099714.1: no array (units 167) at k=51 -> single A (beta-sat terminal) (units 129) at k=31
- chr4_hap2_NA18879#2#CM099729.1: no array (units 46) at k=51 -> single A (beta-sat terminal) (units 41) at k=31
- chr4_hap1_NA18906#1#CM091934.1: no array (units 49) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr4_hap2_NA18943#2#CM101687.1: no array (units 58) at k=51 -> single B (no beta-sat terminal) (units 53) at k=31
- chr4_hap1_NA18948#1#CM101598.1: no array (units 52) at k=51 -> single A (beta-sat terminal) (units 42) at k=31
- chr4_hap2_NA18948#2#CM101587.1: no array (units 51) at k=51 -> single A (beta-sat terminal) (units 46) at k=31
- chr4_hap1_NA18960#1#CM101529.1: no array (units 32) at k=51 -> single A (beta-sat terminal) (units 27) at k=31
- chr4_hap2_NA18960#2#CM101557.1: no array (units 46) at k=51 -> single B (no beta-sat terminal) (units 38) at k=31
- chr4_hap1_NA18967#1#CM101496.1: no array (units 42) at k=51 -> single B (no beta-sat terminal) (units 33) at k=31
- chr4_hap1_NA18971#1#CM087293.1: no array (units 47) at k=51 -> single B (no beta-sat terminal) (units 38) at k=31
- chr4_hap2_NA18971#2#CM087280.1: no array (units 48) at k=51 -> single B (no beta-sat terminal) (units 43) at k=31
- chr4_hap1_NA18974#1#CM089772.1: no array (units 36) at k=51 -> single B (no beta-sat terminal) (units 31) at k=31
- chr4_hap2_NA18974#2#CM089797.1: no array (units 55) at k=51 -> single B (no beta-sat terminal) (units 49) at k=31
- chr4_hap1_NA18982#1#CM101419.1: no array (units 89) at k=51 -> single B (no beta-sat terminal) (units 83) at k=31
- chr4_hap2_NA18982#2#CM101452.1: no array (units 33) at k=51 -> single A (beta-sat terminal) (units 33) at k=31
- chr4_hap1_NA18983#1#CM087419.1: no array (units 62) at k=51 -> single B (no beta-sat terminal) (units 47) at k=31
- chr4_hap2_NA18983#2#CM087431.1: no array (units 38) at k=51 -> single B (no beta-sat terminal) (units 32) at k=31
- chr4_hap2_NA19043#2#CM087071.1: no array (units 168) at k=51 -> single A (beta-sat terminal) (units 128) at k=31
- chr4_hap1_NA19087#1#CM087694.1: no array (units 27) at k=51 -> single B (no beta-sat terminal) (units 27) at k=31
- chr4_hap2_NA19087#2#CM087706.1: no array (units 63) at k=51 -> single B (no beta-sat terminal) (units 51) at k=31
- chr4_hap1_NA19159#1#CM099942.1: no array (units 122) at k=51 -> single A (beta-sat terminal) (units 104) at k=31
- chr4_hap2_NA19159#2#CM099961.1: no array (units 29) at k=51 -> single B (no beta-sat terminal) (units 20) at k=31
- chr4_hap1_NA19185#1#CM087722.1: no array (units 111) at k=51 -> single A (beta-sat terminal) (units 79) at k=31
- chr4_hap1_NA19338#1#CM087756.1: no array (units 67) at k=51 -> single A (beta-sat terminal) (units 62) at k=31
- chr4_hap2_NA19338#2#CM087773.1: no array (units 108) at k=51 -> multi-array (units 93) at k=31
- chr4_hap2_NA19391#2#CM087802.1: no array (units 69) at k=51 -> single A (beta-sat terminal) (units 56) at k=31
- chr4_hap1_NA19700#1#CM089619.1: no array (units 53) at k=51 -> single A (beta-sat terminal) (units 51) at k=31
- chr4_hap2_NA19700#2#CM089636.1: no array (units 48) at k=51 -> single A (beta-sat terminal) (units 45) at k=31
- chr4_hap1_NA19776#1#CM094010.1: no array (units 51) at k=51 -> single B (no beta-sat terminal) (units 48) at k=31
- chr4_hap2_NA19776#2#CM094047.1: no array (units 51) at k=51 -> single B (no beta-sat terminal) (units 42) at k=31
- chr4_hap1_NA19835#1#CM094024.1: no array (units 46) at k=51 -> single A (beta-sat terminal) (units 32) at k=31
- chr4_hap2_NA19909#2#CM094133.1: no array (units 90) at k=51 -> single A (beta-sat terminal) (units 75) at k=31
- chr4_hap2_NA20129#2#CM092119.1: no array (units 61) at k=51 -> single A (beta-sat terminal) (units 50) at k=31
- chr4_hap1_NA20503#1#CM102457.1: no array (units 64) at k=51 -> single B (no beta-sat terminal) (units 59) at k=31
- chr4_hap1_NA20752#1#CM099768.1: no array (units 65) at k=51 -> single B (no beta-sat terminal) (units 51) at k=31
- chr4_hap2_NA20752#2#CM099781.1: no array (units 37) at k=51 -> single B (no beta-sat terminal) (units 31) at k=31
- chr4_hap1_NA20762#1#CM102510.1: no array (units 28) at k=51 -> single A (beta-sat terminal) (units 24) at k=31
- chr4_hap1_NA20799#1#CM098686.1: no array (units 29) at k=51 -> single A (beta-sat terminal) (units 27) at k=31
- chr4_hap2_NA20799#2#CM098715.1: no array (units 25) at k=51 -> single A (beta-sat terminal) (units 18) at k=31
- chr4_hap2_NA20805#2#CM091897.1: no array (units 94) at k=51 -> single B (no beta-sat terminal) (units 76) at k=31
- chr4_hap1_NA20806#1#CM102440.1: no array (units 42) at k=51 -> single B (no beta-sat terminal) (units 36) at k=31
- chr4_hap2_NA20806#2#CM102494.1: no array (units 78) at k=51 -> single A (beta-sat terminal) (units 82) at k=31
- chr4_hap1_NA20827#1#CM102542.1: no array (units 84) at k=51 -> single A (beta-sat terminal) (units 86) at k=31
- chr4_hap2_NA20827#2#CM102562.1: no array (units 23) at k=51 -> single B (no beta-sat terminal) (units 17) at k=31
- chr4_hap1_NA20850#1#CM094291.1: no array (units 94) at k=51 -> single B (no beta-sat terminal) (units 84) at k=31
- chr4_hap2_NA20850#2#CM094305.1: no array (units 45) at k=51 -> single B (no beta-sat terminal) (units 29) at k=31
- chr4_hap1_NA21102#1#CM094425.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 49) at k=31
- chr4_hap2_NA21102#2#CM094446.1: no array (units 45) at k=51 -> single B (no beta-sat terminal) (units 39) at k=31
- chr4_hap1_NA21110#1#CM089665.1: no array (units 33) at k=51 -> single B (no beta-sat terminal) (units 25) at k=31
- chr4_hap2_NA21110#2#CM089681.1: no array (units 58) at k=51 -> single B (no beta-sat terminal) (units 60) at k=31
- chr4_hap1_NA21144#1#CM094093.1: no array (units 32) at k=51 -> single B (no beta-sat terminal) (units 31) at k=31
- chr4_hap2_NA21144#2#CM094109.1: no array (units 54) at k=51 -> single A (beta-sat terminal) (units 48) at k=31
- chr4_hap1_NA21309#1#CM092107.1: no array (units 43) at k=51 -> single A (beta-sat terminal) (units 44) at k=31
- chr4_hap2_NA21309#2#CM092098.1: no array (units 81) at k=51 -> single A (beta-sat terminal) (units 71) at k=31

### k=21 vs k=51: 577 haplotype(s)

- chr10_hap1_HG00097#1#CM094066.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 80) at k=51
- chr10_hap2_HG00097#2#CM094081.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 66) at k=51
- chr10_hap1_HG00099#1#CM087321.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 28) at k=51
- chr10_hap2_HG00099#2#CM087365.1: single A (beta-sat terminal) (units 36) at k=21 -> no array (units 129) at k=51
- chr10_hap2_HG00126#2#CM090128.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 40) at k=51
- chr10_hap1_HG00128#1#CM090084.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 26) at k=51
- chr10_hap2_HG00128#2#CM090099.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 47) at k=51
- chr10_hap2_HG00133#2#CM090068.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 51) at k=51
- chr10_hap1_HG00140#1#CM087115.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 64) at k=51
- chr10_hap2_HG00140#2#CM087140.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 53) at k=51
- chr10_hap2_HG00146#2#CM090038.1: multi-array (units 26) at k=21 -> no array (units 100) at k=51
- chr10_hap1_HG00253#1#CM094463.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 53) at k=51
- chr10_hap2_HG00253#2#CM094485.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 48) at k=51
- chr10_hap1_HG00272#1#CM094193.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 66) at k=51
- chr10_hap2_HG00272#2#CM094220.1: single B (no beta-sat terminal) (units 18) at k=21 -> no array (units 34) at k=51
- chr10_hap1_HG00280#1#CM087153.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 96) at k=51
- chr10_hap2_HG00280#2#CM087164.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 70) at k=51
- chr10_hap2_HG00320#2#CM089969.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 62) at k=51
- chr10_hap1_HG00321#1#CM089545.1: single A (beta-sat terminal) (units 27) at k=21 -> no array (units 85) at k=51
- chr10_hap2_HG00321#2#CM089576.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 56) at k=51
- chr10_hap1_HG00329#1#CM094230.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 46) at k=51
- chr10_hap2_HG00344#2#CM094516.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 64) at k=51
- chr10_hap1_HG00350#1#CM094475.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 22) at k=51
- chr10_hap2_HG00350#2#CM094494.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 64) at k=51
- chr10_hap1_HG00423#1#CM089404.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 85) at k=51
- chr10_hap2_HG00423#2#CM089419.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 50) at k=51
- chr10_hap1_HG00438#1#CM089176.1: degraded / no terminal (units 21) at k=21 -> no array (units 28) at k=51
- chr10_hap2_HG00438#2#CM089194.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 48) at k=51
- chr10_hap1_HG005#1#CM099990.1: single A (beta-sat terminal) (units 37) at k=21 -> no array (units 114) at k=51
- chr10_hap2_HG005#2#CM099979.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 23) at k=51
- chr10_hap2_HG00558#2#CM088528.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 62) at k=51
- chr10_hap1_HG00597#1#CM085773.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 30) at k=51
- chr10_hap2_HG00597#2#CM085787.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 81) at k=51
- chr10_hap1_HG00609#1#CM086353.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 39) at k=51
- chr10_hap2_HG00609#2#CM086365.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 51) at k=51
- chr10_hap1_HG00621#1#CM087855.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 30) at k=51
- chr10_hap2_HG00621#2#CM087868.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 50) at k=51
- chr10_hap1_HG00639#1#CM088544.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 48) at k=51
- chr10_hap2_HG00639#2#CM088558.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 43) at k=51
- chr10_hap2_HG00642#2#CM086624.1: single B (no beta-sat terminal) (units 40) at k=21 -> no array (units 82) at k=51
- chr10_hap1_HG00658#1#CM088780.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 57) at k=51
- chr10_hap2_HG00658#2#CM088802.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 40) at k=51
- chr10_hap1_HG00673#1#CM089148.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 34) at k=51
- chr10_hap2_HG00673#2#CM089161.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 38) at k=51
- chr10_hap2_HG00706#2#CM100064.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 25) at k=51
- chr10_hap2_HG00733#2#CM099645.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 32) at k=51
- chr10_hap1_HG00735#1#CM087882.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG00735#2#CM087891.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 21) at k=51
- chr10_hap1_HG00738#1#CM086649.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 70) at k=51
- chr10_hap1_HG00741#1#CM087902.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 26) at k=51
- chr10_hap2_HG00741#2#CM087915.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 42) at k=51
- chr10_hap1_HG01071#1#CM089119.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 85) at k=51
- chr10_hap2_HG01071#2#CM089134.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 26) at k=51
- chr10_hap1_HG01074#1#CM088790.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 26) at k=51
- chr10_hap2_HG01074#2#CM088816.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 48) at k=51
- chr10_hap1_HG01081#1#CM088571.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 24) at k=51
- chr10_hap2_HG01081#2#CM088587.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 55) at k=51
- chr10_hap1_HG01099#1#CM086263.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 86) at k=51
- chr10_hap1_HG01106#1#CM087956.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 56) at k=51
- chr10_hap2_HG01106#2#CM087968.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 44) at k=51
- chr10_hap1_HG01109#1#CM092091.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 33) at k=51
- chr10_hap2_HG01109#2#CM092079.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 27) at k=51
- chr10_hap1_HG01123#1#CM089087.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 27) at k=51
- chr10_hap1_HG01150#1#CM088602.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 30) at k=51
- chr10_hap2_HG01150#2#CM088618.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 35) at k=51
- chr10_hap2_HG01167#2#CM094546.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 56) at k=51
- chr10_hap1_HG01175#1#CM087929.1: single B (no beta-sat terminal) (units 21) at k=21 -> no array (units 47) at k=51
- chr10_hap2_HG01175#2#CM087943.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 28) at k=51
- chr10_hap2_HG01192#2#CM085945.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 45) at k=51
- chr10_hap1_HG01243#1#CM099607.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG01243#2#CM099595.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 35) at k=51
- chr10_hap1_HG01252#1#CM088828.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 29) at k=51
- chr10_hap2_HG01252#2#CM088859.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 22) at k=51
- chr10_hap1_HG01255#1#CM086711.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 39) at k=51
- chr10_hap2_HG01255#2#CM086724.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 25) at k=51
- chr10_hap2_HG01261#2#CM085921.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 69) at k=51
- chr10_hap1_HG01358#1#CM089059.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 30) at k=51
- chr10_hap2_HG01358#2#CM089072.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 19) at k=51
- chr10_hap1_HG01361#1#CM089027.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 24) at k=51
- chr10_hap2_HG01361#2#CM089043.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 29) at k=51
- chr10_hap2_HG01433#2#CM086515.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 25) at k=51
- chr10_hap2_HG01496#2#CM086392.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 44) at k=51
- chr10_hap1_HG01530#1#CM087539.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG01530#2#CM087562.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 76) at k=51
- chr10_hap2_HG01784#2#CM087601.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 54) at k=51
- chr10_hap1_HG01786#1#CM089508.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 35) at k=51
- chr10_hap2_HG01786#2#CM089532.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 60) at k=51
- chr10_hap1_HG01884#1#CM086409.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 34) at k=51
- chr10_hap1_HG01928#1#CM088992.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 66) at k=51
- chr10_hap2_HG01928#2#CM089008.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 32) at k=51
- chr10_hap1_HG01934#1#CM088844.1: single A (beta-sat terminal) (units 6) at k=21 -> no array (units 23) at k=51
- chr10_hap2_HG01934#2#CM088885.1: single A (beta-sat terminal) (units 5) at k=21 -> no array (units 18) at k=51
- chr10_hap2_HG01940#2#CM088665.1: single B (no beta-sat terminal) (units 20) at k=21 -> no array (units 39) at k=51
- chr10_hap1_HG01943#1#CM089355.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 60) at k=51
- chr10_hap2_HG01943#2#CM089364.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 41) at k=51
- chr10_hap1_HG01952#1#CM087978.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 50) at k=51
- chr10_hap2_HG01952#2#CM088003.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 26) at k=51
- chr10_hap1_HG01960#1#CM088649.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 37) at k=51
- chr10_hap1_HG01969#1#CM088691.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 53) at k=51
- chr10_hap2_HG01969#2#CM088703.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 26) at k=51
- chr10_hap1_HG01975#1#CM099801.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 23) at k=51
- chr10_hap2_HG01975#2#CM099816.1: single A (beta-sat terminal) (units 27) at k=21 -> no array (units 78) at k=51
- chr10_hap1_HG01978#1#CM089250.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 56) at k=51
- chr10_hap2_HG01978#2#CM089276.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 55) at k=51
- chr10_hap1_HG01981#1#CM086325.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 33) at k=51
- chr10_hap1_HG01993#1#CM086498.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 52) at k=51
- chr10_hap2_HG01993#2#CM086538.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 46) at k=51
- chr10_hap1_HG02004#1#CM086305.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 58) at k=51
- chr10_hap2_HG02004#2#CM086314.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 56) at k=51
- chr10_hap1_HG02027#1#CM086979.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 41) at k=51
- chr10_hap2_HG02027#2#CM086991.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 37) at k=51
- chr10_hap2_HG02040#2#CM087408.1: single A (beta-sat terminal) (units 5) at k=21 -> no array (units 17) at k=51
- chr10_hap1_HG02055#1#CM092210.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 77) at k=51
- chr10_hap2_HG02055#2#CM092198.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 66) at k=51
- chr10_hap1_HG02056#1#CM085874.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 31) at k=51
- chr10_hap2_HG02056#2#CM085886.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 36) at k=51
- chr10_hap1_HG02071#1#CM088264.1: single A (beta-sat terminal) (units 6) at k=21 -> no array (units 24) at k=51
- chr10_hap2_HG02071#2#CM088287.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 26) at k=51
- chr10_hap1_HG02080#1#CM092065.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 46) at k=51
- chr10_hap2_HG02080#2#CM092051.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 61) at k=51
- chr10_hap1_HG02109#1#CM092165.1: single A (beta-sat terminal) (units 23) at k=21 -> no array (units 76) at k=51
- chr10_hap1_HG02129#1#CM085853.1: single A (beta-sat terminal) (units 31) at k=21 -> no array (units 121) at k=51
- chr10_hap2_HG02132#2#CM086911.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 31) at k=51
- chr10_hap1_HG02135#1#CM088312.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 43) at k=51
- chr10_hap1_HG02145#1#CM099676.1: single A (beta-sat terminal) (units 32) at k=21 -> no array (units 74) at k=51
- chr10_hap2_HG02145#2#CM099661.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 38) at k=51
- chr10_hap1_HG02148#1#CM087987.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 42) at k=51
- chr10_hap2_HG02148#2#CM088028.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 68) at k=51
- chr10_hap1_HG02155#1#CM091915.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 30) at k=51
- chr10_hap2_HG02155#2#CM091927.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 29) at k=51
- chr10_hap1_HG02165#1#CM087447.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 24) at k=51
- chr10_hap2_HG02165#2#CM087462.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 51) at k=51
- chr10_hap1_HG02178#1#CM089927.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 62) at k=51
- chr10_hap2_HG02178#2#CM089940.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 53) at k=51
- chr10_hap1_HG02257#1#CM089225.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 37) at k=51
- chr10_hap2_HG02257#2#CM089264.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 51) at k=51
- chr10_hap1_HG02258#1#CM085828.1: multi-array (units 15) at k=21 -> no array (units 39) at k=51
- chr10_hap2_HG02258#2#CM085840.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 65) at k=51
- chr10_hap1_HG02273#1#CM087024.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 26) at k=51
- chr10_hap2_HG02273#2#CM087036.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 25) at k=51
- chr10_hap1_HG02280#1#CM086293.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 45) at k=51
- chr10_hap1_HG02293#1#CM086251.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 67) at k=51
- chr10_hap2_HG02293#2#CM086240.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 43) at k=51
- chr10_hap1_HG02300#1#CM086438.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 24) at k=51
- chr10_hap2_HG02300#2#CM086449.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG02391#2#CM087617.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 43) at k=51
- chr10_hap1_HG02392#1#CM099884.1: single B (no beta-sat terminal) (units 25) at k=21 -> no array (units 50) at k=51
- chr10_hap2_HG02451#2#CM086229.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 60) at k=51
- chr10_hap2_HG02486#2#CM088223.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 56) at k=51
- chr10_hap1_HG02514#1#CM098735.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 40) at k=51
- chr10_hap2_HG02514#2#CM098751.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 25) at k=51
- chr10_hap2_HG02523#2#CM086942.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 62) at k=51
- chr10_hap1_HG02559#1#CM088237.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 68) at k=51
- chr10_hap1_HG02572#1#CM088159.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 50) at k=51
- chr10_hap2_HG02572#2#CM088177.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 29) at k=51
- chr10_hap1_HG02602#1#CM099690.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 43) at k=51
- chr10_hap2_HG02602#2#CM099704.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 62) at k=51
- chr10_hap1_HG02615#1#CM086674.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 32) at k=51
- chr10_hap1_HG02630#1#CM091812.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG02630#2#CM091744.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 33) at k=51
- chr10_hap1_HG02647#1#CM086562.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 92) at k=51
- chr10_hap2_HG02647#2#CM086576.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 71) at k=51
- chr10_hap2_HG02668#2#CM088900.1: multi-array (units 26) at k=21 -> no array (units 70) at k=51
- chr10_hap1_HG02698#1#CM086829.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 35) at k=51
- chr10_hap1_HG02717#1#CM088962.1: single A (beta-sat terminal) (units 41) at k=21 -> no array (units 164) at k=51
- chr10_hap2_HG02717#2#CM088975.1: multi-array (units 38) at k=21 -> no array (units 104) at k=51
- chr10_hap1_HG02723#1#CM092035.1: multi-array (units 13) at k=21 -> no array (units 40) at k=51
- chr10_hap2_HG02723#2#CM092016.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 76) at k=51
- chr10_hap1_HG02735#1#CM089334.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG02735#2#CM089347.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 45) at k=51
- chr10_hap1_HG02738#1#CM099923.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 82) at k=51
- chr10_hap2_HG02738#2#CM099934.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 46) at k=51
- chr10_hap1_HG02809#1#CM088348.1: single A (beta-sat terminal) (units 27) at k=21 -> no array (units 111) at k=51
- chr10_hap2_HG02809#2#CM088367.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 26) at k=51
- chr10_hap2_HG02818#2#CM091991.1: multi-array (units 49) at k=21 -> no array (units 156) at k=51
- chr10_hap2_HG02841#2#CM098675.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 93) at k=51
- chr10_hap2_HG02886#2#CM088131.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 30) at k=51
- chr10_hap1_HG02922#1#CM087002.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 49) at k=51
- chr10_hap1_HG02965#1#CM087350.1: single A (beta-sat terminal) (units 32) at k=21 -> no array (units 96) at k=51
- chr10_hap2_HG02976#2#CM087126.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 45) at k=51
- chr10_hap1_HG02984#1#CM098794.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 57) at k=51
- chr10_hap2_HG02984#2#CM098808.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 50) at k=51
- chr10_hap1_HG03017#1#CM088915.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 22) at k=51
- chr10_hap2_HG03017#2#CM088932.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 97) at k=51
- chr10_hap1_HG03041#1#CM088718.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 41) at k=51
- chr10_hap1_HG03050#1#CM098767.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 44) at k=51
- chr10_hap2_HG03050#2#CM098781.1: single A (beta-sat terminal) (units 42) at k=21 -> no array (units 140) at k=51
- chr10_hap1_HG03098#1#CM092235.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 42) at k=51
- chr10_hap1_HG03130#1#CM087259.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 63) at k=51
- chr10_hap2_HG03130#2#CM087273.1: single A (beta-sat terminal) (units 27) at k=21 -> no array (units 98) at k=51
- chr10_hap2_HG03139#2#CM091841.1: single A (beta-sat terminal) (units 23) at k=21 -> no array (units 71) at k=51
- chr10_hap1_HG03209#1#CM087217.1: multi-array (units 28) at k=21 -> no array (units 57) at k=51
- chr10_hap2_HG03225#2#CM089321.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 24) at k=51
- chr10_hap1_HG03239#1#CM088751.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 78) at k=51
- chr10_hap2_HG03239#2#CM088765.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 38) at k=51
- chr10_hap1_HG03270#1#CM089904.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 33) at k=51
- chr10_hap2_HG03270#2#CM089916.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 30) at k=51
- chr10_hap1_HG03369#1#CM094387.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 74) at k=51
- chr10_hap2_HG03369#2#CM094363.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 64) at k=51
- chr10_hap1_HG03453#1#CM088015.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 71) at k=51
- chr10_hap2_HG03453#2#CM088039.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 117) at k=51
- chr10_hap1_HG03470#1#CM091826.1: single A (beta-sat terminal) (units 23) at k=21 -> no array (units 71) at k=51
- chr10_hap2_HG03470#2#CM091874.1: multi-array (units 29) at k=21 -> no array (units 89) at k=51
- chr10_hap1_HG03471#1#CM087526.1: degraded / no terminal (units 16) at k=21 -> no array (units 53) at k=51
- chr10_hap2_HG03471#2#CM087549.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 71) at k=51
- chr10_hap1_HG03486#1#CM091978.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 39) at k=51
- chr10_hap2_HG03486#2#CM091964.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 62) at k=51
- chr10_hap1_HG03516#1#CM089209.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 69) at k=51
- chr10_hap1_HG03521#1#CM094417.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 71) at k=51
- chr10_hap1_HG03540#1#CM088053.1: degraded / no terminal (units 31) at k=21 -> no array (units 102) at k=51
- chr10_hap1_HG03579#1#CM088949.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 57) at k=51
- chr10_hap2_HG03579#2#CM089292.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 87) at k=51
- chr10_hap1_HG03654#1#CM086839.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 36) at k=51
- chr10_hap2_HG03654#2#CM086859.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 35) at k=51
- chr10_hap1_HG03669#1#CM086736.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 54) at k=51
- chr10_hap2_HG03669#2#CM086747.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 60) at k=51
- chr10_hap1_HG03688#1#CM086806.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 74) at k=51
- chr10_hap1_HG03704#1#CM086215.1: single B (no beta-sat terminal) (units 29) at k=21 -> no array (units 55) at k=51
- chr10_hap2_HG03704#2#CM086193.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 58) at k=51
- chr10_hap1_HG03742#1#CM087632.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 31) at k=51
- chr10_hap2_HG03742#2#CM087647.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 60) at k=51
- chr10_hap1_HG03784#1#CM094178.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 41) at k=51
- chr10_hap2_HG03784#2#CM094205.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 88) at k=51
- chr10_hap2_HG03804#2#CM088383.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 65) at k=51
- chr10_hap1_HG03816#1#CM100031.1: single A (beta-sat terminal) (units 33) at k=21 -> no array (units 98) at k=51
- chr10_hap2_HG03816#2#CM100046.1: single B (no beta-sat terminal) (units 33) at k=21 -> no array (units 48) at k=51
- chr10_hap1_HG03831#1#CM086529.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 67) at k=51
- chr10_hap1_HG03834#1#CM085800.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 52) at k=51
- chr10_hap1_HG03874#1#CM089436.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 25) at k=51
- chr10_hap2_HG03874#2#CM089478.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 38) at k=51
- chr10_hap1_HG03927#1#CM086638.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 29) at k=51
- chr10_hap2_HG03927#2#CM086663.1: single A (beta-sat terminal) (units 6) at k=21 -> no array (units 21) at k=51
- chr10_hap2_HG03942#2#CM088407.1: single B (no beta-sat terminal) (units 23) at k=21 -> no array (units 35) at k=51
- chr10_hap1_HG04115#1#CM086757.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 60) at k=51
- chr10_hap1_HG04157#1#CM088417.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 66) at k=51
- chr10_hap1_HG04160#1#CM088430.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 25) at k=51
- chr10_hap2_HG04160#2#CM088451.1: single B (no beta-sat terminal) (units 18) at k=21 -> no array (units 39) at k=51
- chr10_hap1_HG04184#1#CM086463.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 45) at k=51
- chr10_hap1_HG04187#1#CM099746.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 86) at k=51
- chr10_hap2_HG04187#2#CM099760.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 72) at k=51
- chr10_hap2_HG04199#2#CM086793.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 50) at k=51
- chr10_hap2_HG04204#2#CM088472.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 78) at k=51
- chr10_hap1_HG04228#1#CM088484.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 36) at k=51
- chr10_hap2_HG04228#2#CM088512.1: single B (no beta-sat terminal) (units 36) at k=21 -> no array (units 67) at k=51
- chr10_hap1_HG06807#1#CM101404.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 50) at k=51
- chr10_hap2_HG06807#2#CM101438.1: multi-array (units 15) at k=21 -> no array (units 50) at k=51
- chr10_hap2_NA18522#2#CM087194.1: multi-array (units 18) at k=21 -> no array (units 57) at k=51
- chr10_hap1_NA18565#1#CM094322.1: multi-array (units 29) at k=21 -> no array (units 100) at k=51
- chr10_hap2_NA18565#2#CM094344.1: single A (beta-sat terminal) (units 6) at k=21 -> no array (units 31) at k=51
- chr10_hap1_NA18570#1#CM087663.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 35) at k=51
- chr10_hap2_NA18570#2#CM087681.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 42) at k=51
- chr10_hap1_NA18608#1#CM089843.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 50) at k=51
- chr10_hap1_NA18620#1#CM089823.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 55) at k=51
- chr10_hap2_NA18620#2#CM089833.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 64) at k=51
- chr10_hap2_NA18747#2#CM087245.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 18) at k=51
- chr10_hap1_NA18879#1#CM099718.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 59) at k=51
- chr10_hap2_NA18879#2#CM099733.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 48) at k=51
- chr10_hap1_NA18906#1#CM091938.1: single A (beta-sat terminal) (units 36) at k=21 -> no array (units 123) at k=51
- chr10_hap1_NA18943#1#CM101663.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 39) at k=51
- chr10_hap2_NA18943#2#CM101692.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 40) at k=51
- chr10_hap1_NA18944#1#CM101649.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 45) at k=51
- chr10_hap1_NA18945#1#CM101617.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 60) at k=51
- chr10_hap2_NA18945#2#CM101632.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 37) at k=51
- chr10_hap1_NA18948#1#CM101603.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 53) at k=51
- chr10_hap1_NA18952#1#CM089788.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 25) at k=51
- chr10_hap2_NA18959#2#CM101577.1: single A (beta-sat terminal) (units 6) at k=21 -> no array (units 20) at k=51
- chr10_hap2_NA18967#2#CM101516.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 41) at k=51
- chr10_hap1_NA18970#1#CM101471.1: single B (no beta-sat terminal) (units 47) at k=21 -> no array (units 93) at k=51
- chr10_hap2_NA18970#2#CM101483.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 26) at k=51
- chr10_hap1_NA18971#1#CM087296.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 26) at k=51
- chr10_hap2_NA18971#2#CM087284.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 48) at k=51
- chr10_hap1_NA18974#1#CM089776.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 44) at k=51
- chr10_hap2_NA18974#2#CM089800.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 35) at k=51
- chr10_hap1_NA18976#1#CM089755.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 66) at k=51
- chr10_hap2_NA18976#2#CM089764.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 65) at k=51
- chr10_hap1_NA18982#1#CM101423.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 22) at k=51
- chr10_hap2_NA18982#2#CM101456.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 53) at k=51
- chr10_hap1_NA18983#1#CM087423.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 42) at k=51
- chr10_hap2_NA18983#2#CM087435.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 42) at k=51
- chr10_hap1_NA19036#1#CM089728.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 25) at k=51
- chr10_hap1_NA19043#1#CM087049.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 100) at k=51
- chr10_hap2_NA19043#2#CM087075.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 92) at k=51
- chr10_hap1_NA19087#1#CM087700.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 45) at k=51
- chr10_hap2_NA19087#2#CM087711.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 37) at k=51
- chr10_hap1_NA19159#1#CM099947.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 39) at k=51
- chr10_hap2_NA19159#2#CM099966.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 55) at k=51
- chr10_hap1_NA19185#1#CM087728.1: multi-array (units 54) at k=21 -> no array (units 118) at k=51
- chr10_hap2_NA19185#2#CM087744.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 59) at k=51
- chr10_hap1_NA19240#1#CM099581.1: multi-array (units 32) at k=21 -> no array (units 89) at k=51
- chr10_hap2_NA19240#2#CM099619.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 83) at k=51
- chr10_hap1_NA19338#1#CM087762.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 55) at k=51
- chr10_hap2_NA19338#2#CM087779.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 43) at k=51
- chr10_hap1_NA19443#1#CM089450.1: single A (beta-sat terminal) (units 34) at k=21 -> no array (units 130) at k=51
- chr10_hap1_NA19468#1#CM087822.1: multi-array (units 30) at k=21 -> no array (units 84) at k=51
- chr10_hap2_NA19468#2#CM087837.1: multi-array (units 24) at k=21 -> no array (units 63) at k=51
- chr10_hap1_NA19700#1#CM089625.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 71) at k=51
- chr10_hap2_NA19700#2#CM089642.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 98) at k=51
- chr10_hap2_NA19776#2#CM094053.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 52) at k=51
- chr10_hap2_NA19835#2#CM094038.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 74) at k=51
- chr10_hap1_NA19909#1#CM094123.1: single A (beta-sat terminal) (units 6) at k=21 -> no array (units 24) at k=51
- chr10_hap2_NA19909#2#CM094136.1: single A (beta-sat terminal) (units 30) at k=21 -> no array (units 121) at k=51
- chr10_hap1_NA20129#1#CM092136.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 51) at k=51
- chr10_hap2_NA20282#2#CM094283.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 49) at k=51
- chr10_hap1_NA20346#1#CM094253.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 39) at k=51
- chr10_hap2_NA20346#2#CM094265.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 36) at k=51
- chr10_hap1_NA20503#1#CM102463.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 45) at k=51
- chr10_hap2_NA20503#2#CM102480.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 39) at k=51
- chr10_hap2_NA20752#2#CM099786.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 51) at k=51
- chr10_hap1_NA20762#1#CM102516.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 70) at k=51
- chr10_hap2_NA20799#2#CM098721.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 74) at k=51
- chr10_hap1_NA20805#1#CM091890.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 56) at k=51
- chr10_hap2_NA20805#2#CM091902.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 36) at k=51
- chr10_hap1_NA20806#1#CM102446.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 67) at k=51
- chr10_hap1_NA20827#1#CM102548.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 25) at k=51
- chr10_hap2_NA20827#2#CM102568.1: single A (beta-sat terminal) (units 8) at k=21 -> no array (units 23) at k=51
- chr10_hap2_NA20850#2#CM094309.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 105) at k=51
- chr10_hap1_NA20870#1#CM089701.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 41) at k=51
- chr10_hap2_NA20870#2#CM089715.1: single A (beta-sat terminal) (units 4) at k=21 -> no array (units 18) at k=51
- chr10_hap1_NA20905#1#CM087477.1: single A (beta-sat terminal) (units 11) at k=21 -> no array (units 33) at k=51
- chr10_hap2_NA20905#2#CM087511.1: single A (beta-sat terminal) (units 9) at k=21 -> no array (units 30) at k=51
- chr10_hap1_NA21093#1#CM089601.1: single A (beta-sat terminal) (units 31) at k=21 -> no array (units 141) at k=51
- chr10_hap2_NA21093#2#CM089611.1: single A (beta-sat terminal) (units 37) at k=21 -> no array (units 141) at k=51
- chr10_hap1_NA21106#1#CM089560.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 52) at k=51
- chr10_hap2_NA21106#2#CM089591.1: single B (no beta-sat terminal) (units 17) at k=21 -> no array (units 35) at k=51
- chr10_hap1_NA21110#1#CM089669.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 80) at k=51
- chr10_hap2_NA21110#2#CM089685.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 55) at k=51
- chr10_hap1_NA21144#1#CM094098.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 44) at k=51
- chr10_hap2_NA21144#2#CM094114.1: single A (beta-sat terminal) (units 31) at k=21 -> no array (units 84) at k=51
- chr10_hap1_NA21309#1#CM092113.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 54) at k=51
- chr4_hap2_HG00099#2#CM087359.1: single B (no beta-sat terminal) (units 28) at k=21 -> no array (units 56) at k=51
- chr4_hap1_HG00128#1#CM090078.1: single A (beta-sat terminal) (units 36) at k=21 -> no array (units 72) at k=51
- chr4_hap2_HG00128#2#CM090094.1: single B (no beta-sat terminal) (units 18) at k=21 -> no array (units 35) at k=51
- chr4_hap1_HG00133#1#CM090048.1: single A (beta-sat terminal) (units 43) at k=21 -> no array (units 79) at k=51
- chr4_hap2_HG00133#2#CM090063.1: multi-array (units 41) at k=21 -> no array (units 73) at k=51
- chr4_hap1_HG00140#1#CM087110.1: single B (no beta-sat terminal) (units 40) at k=21 -> no array (units 84) at k=51
- chr4_hap2_HG00140#2#CM087135.1: single B (no beta-sat terminal) (units 41) at k=21 -> no array (units 82) at k=51
- chr4_hap2_HG00146#2#CM090035.1: single B (no beta-sat terminal) (units 38) at k=21 -> no array (units 71) at k=51
- chr4_hap1_HG00272#1#CM094187.1: single B (no beta-sat terminal) (units 51) at k=21 -> no array (units 102) at k=51
- chr4_hap2_HG00272#2#CM094214.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 49) at k=51
- chr4_hap1_HG00280#1#CM087148.1: single A (beta-sat terminal) (units 35) at k=21 -> no array (units 58) at k=51
- chr4_hap2_HG00280#2#CM087159.1: single A (beta-sat terminal) (units 140) at k=21 -> no array (units 274) at k=51
- chr4_hap1_HG00320#1#CM089949.1: single A (beta-sat terminal) (units 30) at k=21 -> no array (units 64) at k=51
- chr4_hap2_HG00320#2#CM089963.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 52) at k=51
- chr4_hap1_HG00321#1#CM089539.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 50) at k=51
- chr4_hap2_HG00321#2#CM089570.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 61) at k=51
- chr4_hap2_HG00408#2#CM085967.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 44) at k=51
- chr4_hap1_HG00423#1#CM089399.1: single B (no beta-sat terminal) (units 14) at k=21 -> no array (units 29) at k=51
- chr4_hap2_HG00423#2#CM089415.1: single B (no beta-sat terminal) (units 15) at k=21 -> no array (units 29) at k=51
- chr4_hap1_HG00438#1#CM089170.1: single A (beta-sat terminal) (units 18) at k=21 -> no array (units 35) at k=51
- chr4_hap2_HG00438#2#CM089189.1: single A (beta-sat terminal) (units 33) at k=21 -> no array (units 66) at k=51
- chr4_hap1_HG00558#1#CM088493.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 52) at k=51
- chr4_hap2_HG00558#2#CM088522.1: single A (beta-sat terminal) (units 40) at k=21 -> no array (units 76) at k=51
- chr4_hap2_HG00597#2#CM085782.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 56) at k=51
- chr4_hap1_HG00609#1#CM086349.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 45) at k=51
- chr4_hap2_HG00609#2#CM086359.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 47) at k=51
- chr4_hap1_HG00621#1#CM087849.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 37) at k=51
- chr4_hap2_HG00621#2#CM087863.1: single B (no beta-sat terminal) (units 20) at k=21 -> no array (units 35) at k=51
- chr4_hap1_HG00639#1#CM088538.1: single A (beta-sat terminal) (units 44) at k=21 -> no array (units 93) at k=51
- chr4_hap2_HG00639#2#CM088553.1: single A (beta-sat terminal) (units 37) at k=21 -> no array (units 61) at k=51
- chr4_hap2_HG00642#2#CM086620.1: single B (no beta-sat terminal) (units 30) at k=21 -> no array (units 62) at k=51
- chr4_hap1_HG00673#1#CM089144.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 50) at k=51
- chr4_hap2_HG00673#2#CM089155.1: single B (no beta-sat terminal) (units 29) at k=21 -> no array (units 57) at k=51
- chr4_hap1_HG00706#1#CM100014.1: single B (no beta-sat terminal) (units 26) at k=21 -> no array (units 52) at k=51
- chr4_hap2_HG00706#2#CM100058.1: single B (no beta-sat terminal) (units 28) at k=21 -> no array (units 54) at k=51
- chr4_hap1_HG00735#1#CM087877.1: single A (beta-sat terminal) (units 40) at k=21 -> no array (units 68) at k=51
- chr4_hap2_HG00735#2#CM087886.1: single B (no beta-sat terminal) (units 54) at k=21 -> no array (units 87) at k=51
- chr4_hap2_HG00738#2#CM086682.1: single B (no beta-sat terminal) (units 34) at k=21 -> no array (units 59) at k=51
- chr4_hap1_HG01071#1#CM089114.1: single A (beta-sat terminal) (units 34) at k=21 -> no array (units 62) at k=51
- chr4_hap2_HG01071#2#CM089130.1: single A (beta-sat terminal) (units 47) at k=21 -> no array (units 86) at k=51
- chr4_hap1_HG01109#1#CM092086.1: multi-array (units 38) at k=21 -> no array (units 71) at k=51
- chr4_hap1_HG01167#1#CM094526.1: single B (no beta-sat terminal) (units 27) at k=21 -> no array (units 55) at k=51
- chr4_hap2_HG01167#2#CM094541.1: single B (no beta-sat terminal) (units 35) at k=21 -> no array (units 69) at k=51
- chr4_hap1_HG01175#1#CM087923.1: single B (no beta-sat terminal) (units 37) at k=21 -> no array (units 71) at k=51
- chr4_hap2_HG01175#2#CM087938.1: single A (beta-sat terminal) (units 58) at k=21 -> no array (units 105) at k=51
- chr4_hap2_HG01192#2#CM085941.1: single A (beta-sat terminal) (units 38) at k=21 -> no array (units 69) at k=51
- chr4_hap1_HG01255#1#CM086705.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 45) at k=51
- chr4_hap1_HG01258#1#CM088111.1: single B (no beta-sat terminal) (units 26) at k=21 -> no array (units 50) at k=51
- chr4_hap2_HG01258#2#CM088140.1: single B (no beta-sat terminal) (units 45) at k=21 -> no array (units 86) at k=51
- chr4_hap1_HG01261#1#CM085908.1: single B (no beta-sat terminal) (units 27) at k=21 -> no array (units 54) at k=51
- chr4_hap1_HG01346#1#CM086586.1: single B (no beta-sat terminal) (units 44) at k=21 -> no array (units 66) at k=51
- chr4_hap2_HG01346#2#CM086599.1: single B (no beta-sat terminal) (units 33) at k=21 -> no array (units 60) at k=51
- chr4_hap1_HG01358#1#CM089053.1: single A (beta-sat terminal) (units 32) at k=21 -> no array (units 57) at k=51
- chr4_hap2_HG01358#2#CM089066.1: single B (no beta-sat terminal) (units 30) at k=21 -> no array (units 63) at k=51
- chr4_hap1_HG01361#1#CM089022.1: single B (no beta-sat terminal) (units 56) at k=21 -> no array (units 114) at k=51
- chr4_hap2_HG01361#2#CM089038.1: single B (no beta-sat terminal) (units 38) at k=21 -> no array (units 60) at k=51
- chr4_hap2_HG01433#2#CM086509.1: single B (no beta-sat terminal) (units 28) at k=21 -> no array (units 51) at k=51
- chr4_hap1_HG01530#1#CM087534.1: single B (no beta-sat terminal) (units 15) at k=21 -> no array (units 30) at k=51
- chr4_hap2_HG01530#2#CM087558.1: single A (beta-sat terminal) (units 46) at k=21 -> no array (units 76) at k=51
- chr4_hap1_HG01784#1#CM087568.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 35) at k=51
- chr4_hap1_HG01786#1#CM089502.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 40) at k=51
- chr4_hap2_HG01891#2#CM088207.1: multi-array (units 36) at k=21 -> no array (units 74) at k=51
- chr4_hap1_HG01928#1#CM088986.1: single B (no beta-sat terminal) (units 30) at k=21 -> no array (units 60) at k=51
- chr4_hap1_HG01940#1#CM088631.1: single B (no beta-sat terminal) (units 26) at k=21 -> no array (units 52) at k=51
- chr4_hap2_HG01940#2#CM088661.1: single B (no beta-sat terminal) (units 21) at k=21 -> no array (units 44) at k=51
- chr4_hap2_HG01952#2#CM087997.1: single B (no beta-sat terminal) (units 34) at k=21 -> no array (units 66) at k=51
- chr4_hap1_HG01981#1#CM086321.1: single B (no beta-sat terminal) (units 31) at k=21 -> no array (units 50) at k=51
- chr4_hap2_HG01981#2#CM086335.1: single B (no beta-sat terminal) (units 61) at k=21 -> no array (units 107) at k=51
- chr4_hap1_HG01993#1#CM086494.1: single B (no beta-sat terminal) (units 35) at k=21 -> no array (units 64) at k=51
- chr4_hap2_HG01993#2#CM086535.1: single B (no beta-sat terminal) (units 35) at k=21 -> no array (units 64) at k=51
- chr4_hap2_HG02004#2#CM086311.1: single B (no beta-sat terminal) (units 46) at k=21 -> no array (units 95) at k=51
- chr4_hap1_HG02015#1#CM085891.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 38) at k=51
- chr4_hap2_HG02015#2#CM085898.1: single A (beta-sat terminal) (units 42) at k=21 -> no array (units 78) at k=51
- chr4_hap1_HG02027#1#CM086976.1: single A (beta-sat terminal) (units 37) at k=21 -> no array (units 71) at k=51
- chr4_hap2_HG02027#2#CM086986.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 42) at k=51
- chr4_hap1_HG02040#1#CM087374.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 43) at k=51
- chr4_hap1_HG02074#1#CM088269.1: single B (no beta-sat terminal) (units 51) at k=21 -> no array (units 95) at k=51
- chr4_hap2_HG02074#2#CM088295.1: single B (no beta-sat terminal) (units 20) at k=21 -> no array (units 39) at k=51
- chr4_hap2_HG02080#2#CM092046.1: single A (beta-sat terminal) (units 15) at k=21 -> no array (units 33) at k=51
- chr4_hap1_HG02083#1#CM086953.1: single A (beta-sat terminal) (units 23) at k=21 -> no array (units 41) at k=51
- chr4_hap1_HG02109#1#CM092159.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 56) at k=51
- chr4_hap2_HG02109#2#CM092145.1: single A (beta-sat terminal) (units 42) at k=21 -> no array (units 87) at k=51
- chr4_hap1_HG02129#1#CM085848.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 38) at k=51
- chr4_hap2_HG02129#2#CM085859.1: single B (no beta-sat terminal) (units 17) at k=21 -> no array (units 37) at k=51
- chr4_hap1_HG02132#1#CM086895.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 40) at k=51
- chr4_hap2_HG02132#2#CM086906.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 35) at k=51
- chr4_hap1_HG02155#1#CM091910.1: single A (beta-sat terminal) (units 45) at k=21 -> no array (units 76) at k=51
- chr4_hap2_HG02155#2#CM091922.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 32) at k=51
- chr4_hap1_HG02165#1#CM087442.1: single B (no beta-sat terminal) (units 26) at k=21 -> no array (units 55) at k=51
- chr4_hap2_HG02165#2#CM087458.1: single B (no beta-sat terminal) (units 29) at k=21 -> no array (units 53) at k=51
- chr4_hap1_HG02178#1#CM089923.1: single B (no beta-sat terminal) (units 14) at k=21 -> no array (units 30) at k=51
- chr4_hap2_HG02178#2#CM089936.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 41) at k=51
- chr4_hap1_HG02257#1#CM089219.1: single A (beta-sat terminal) (units 44) at k=21 -> no array (units 85) at k=51
- chr4_hap2_HG02257#2#CM089258.1: multi-array (units 41) at k=21 -> no array (units 75) at k=51
- chr4_hap1_HG02258#1#CM085824.1: single A (beta-sat terminal) (units 59) at k=21 -> no array (units 106) at k=51
- chr4_hap2_HG02258#2#CM085835.1: single A (beta-sat terminal) (units 52) at k=21 -> no array (units 97) at k=51
- chr4_hap2_HG02280#2#CM086280.1: single B (no beta-sat terminal) (units 30) at k=21 -> no array (units 55) at k=51
- chr4_hap2_HG02391#2#CM087611.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 40) at k=51
- chr4_hap1_HG02392#1#CM099880.1: single B (no beta-sat terminal) (units 20) at k=21 -> no array (units 39) at k=51
- chr4_hap2_HG02451#2#CM086224.1: single A (beta-sat terminal) (units 39) at k=21 -> no array (units 57) at k=51
- chr4_hap2_HG02486#2#CM088217.1: single A (beta-sat terminal) (units 34) at k=21 -> no array (units 59) at k=51
- chr4_hap2_HG02514#2#CM098746.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 36) at k=51
- chr4_hap2_HG02523#2#CM086936.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 29) at k=51
- chr4_hap1_HG02559#1#CM088233.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 36) at k=51
- chr4_hap2_HG02559#2#CM088246.1: single A (beta-sat terminal) (units 37) at k=21 -> no array (units 74) at k=51
- chr4_hap1_HG02572#1#CM088153.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 53) at k=51
- chr4_hap2_HG02572#2#CM088171.1: single A (beta-sat terminal) (units 66) at k=21 -> no array (units 131) at k=51
- chr4_hap1_HG02583#1#CM089488.1: single A (beta-sat terminal) (units 37) at k=21 -> no array (units 74) at k=51
- chr4_hap2_HG02583#2#CM089515.1: single A (beta-sat terminal) (units 36) at k=21 -> no array (units 64) at k=51
- chr4_hap1_HG02602#1#CM099685.1: single B (no beta-sat terminal) (units 23) at k=21 -> no array (units 42) at k=51
- chr4_hap1_HG02615#1#CM086671.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 38) at k=51
- chr4_hap2_HG02615#2#CM086692.1: single A (beta-sat terminal) (units 66) at k=21 -> no array (units 127) at k=51
- chr4_hap2_HG02622#2#CM088097.1: single A (beta-sat terminal) (units 23) at k=21 -> no array (units 43) at k=51
- chr4_hap1_HG02630#1#CM091807.1: single A (beta-sat terminal) (units 33) at k=21 -> no array (units 60) at k=51
- chr4_hap2_HG02630#2#CM091738.1: single B (no beta-sat terminal) (units 13) at k=21 -> no array (units 28) at k=51
- chr4_hap1_HG02647#1#CM086557.1: single A (beta-sat terminal) (units 40) at k=21 -> no array (units 74) at k=51
- chr4_hap2_HG02668#2#CM088895.1: single A (beta-sat terminal) (units 23) at k=21 -> no array (units 68) at k=51
- chr4_hap1_HG02717#1#CM088957.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 45) at k=51
- chr4_hap2_HG02717#2#CM088970.1: single A (beta-sat terminal) (units 33) at k=21 -> no array (units 61) at k=51
- chr4_hap1_HG02723#1#CM092029.1: single A (beta-sat terminal) (units 91) at k=21 -> no array (units 187) at k=51
- chr4_hap2_HG02723#2#CM092011.1: single A (beta-sat terminal) (units 27) at k=21 -> no array (units 52) at k=51
- chr4_hap1_HG02735#1#CM089329.1: single B (no beta-sat terminal) (units 33) at k=21 -> no array (units 68) at k=51
- chr4_hap2_HG02738#2#CM099931.1: single B (no beta-sat terminal) (units 33) at k=21 -> no array (units 64) at k=51
- chr4_hap1_HG02809#1#CM088343.1: single A (beta-sat terminal) (units 17) at k=21 -> no array (units 38) at k=51
- chr4_hap2_HG02809#2#CM088361.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 37) at k=51
- chr4_hap1_HG02841#1#CM098701.1: single A (beta-sat terminal) (units 30) at k=21 -> no array (units 65) at k=51
- chr4_hap1_HG02886#1#CM088084.1: multi-array (units 30) at k=21 -> no array (units 59) at k=51
- chr4_hap2_HG02886#2#CM088125.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 39) at k=51
- chr4_hap2_HG02922#2#CM087011.1: single A (beta-sat terminal) (units 50) at k=21 -> no array (units 105) at k=51
- chr4_hap1_HG02965#1#CM087344.1: single A (beta-sat terminal) (units 35) at k=21 -> no array (units 72) at k=51
- chr4_hap2_HG02965#2#CM087387.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 51) at k=51
- chr4_hap1_HG03017#1#CM088910.1: single B (no beta-sat terminal) (units 55) at k=21 -> no array (units 96) at k=51
- chr4_hap2_HG03017#2#CM088926.1: single B (no beta-sat terminal) (units 31) at k=21 -> no array (units 58) at k=51
- chr4_hap1_HG03041#1#CM088713.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 42) at k=51
- chr4_hap2_HG03041#2#CM088730.1: single A (beta-sat terminal) (units 10) at k=21 -> no array (units 18) at k=51
- chr4_hap1_HG03050#1#CM098761.1: single A (beta-sat terminal) (units 79) at k=21 -> no array (units 157) at k=51
- chr4_hap1_HG03098#1#CM092229.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 43) at k=51
- chr4_hap2_HG03098#2#CM092217.1: single A (beta-sat terminal) (units 47) at k=21 -> no array (units 100) at k=51
- chr4_hap1_HG03130#1#CM087254.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 50) at k=51
- chr4_hap2_HG03130#2#CM087268.1: single A (beta-sat terminal) (units 88) at k=21 -> no array (units 218) at k=51
- chr4_hap1_HG03139#1#CM091853.1: degraded / no terminal (units 51) at k=21 -> no array (units 93) at k=51
- chr4_hap2_HG03139#2#CM091838.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 27) at k=51
- chr4_hap1_HG03209#1#CM087211.1: single A (beta-sat terminal) (units 40) at k=21 -> no array (units 75) at k=51
- chr4_hap2_HG03209#2#CM087229.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 40) at k=51
- chr4_hap1_HG03225#1#CM089299.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 48) at k=51
- chr4_hap1_HG03239#1#CM088746.1: single B (no beta-sat terminal) (units 22) at k=21 -> no array (units 45) at k=51
- chr4_hap1_HG03453#1#CM088011.1: single A (beta-sat terminal) (units 39) at k=21 -> no array (units 81) at k=51
- chr4_hap1_HG03470#1#CM091820.1: single A (beta-sat terminal) (units 12) at k=21 -> no array (units 29) at k=51
- chr4_hap2_HG03470#2#CM091868.1: single A (beta-sat terminal) (units 40) at k=21 -> no array (units 77) at k=51
- chr4_hap1_HG03471#1#CM087521.1: single A (beta-sat terminal) (units 41) at k=21 -> no array (units 78) at k=51
- chr4_hap2_HG03471#2#CM087544.1: single A (beta-sat terminal) (units 36) at k=21 -> no array (units 73) at k=51
- chr4_hap1_HG03486#1#CM091972.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 48) at k=51
- chr4_hap2_HG03486#2#CM091958.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 51) at k=51
- chr4_hap1_HG03516#1#CM089204.1: single B (no beta-sat terminal) (units 11) at k=21 -> no array (units 23) at k=51
- chr4_hap1_HG03521#1#CM094413.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 31) at k=51
- chr4_hap2_HG03521#2#CM094435.1: single A (beta-sat terminal) (units 45) at k=21 -> no array (units 83) at k=51
- chr4_hap1_HG03540#1#CM088048.1: single B (no beta-sat terminal) (units 14) at k=21 -> no array (units 29) at k=51
- chr4_hap2_HG03540#2#CM088061.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 43) at k=51
- chr4_hap1_HG03579#1#CM088943.1: single A (beta-sat terminal) (units 21) at k=21 -> no array (units 45) at k=51
- chr4_hap2_HG03579#2#CM089286.1: single B (no beta-sat terminal) (units 54) at k=21 -> no array (units 105) at k=51
- chr4_hap2_HG03654#2#CM086855.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 48) at k=51
- chr4_hap1_HG03669#1#CM086732.1: single A (beta-sat terminal) (units 51) at k=21 -> no array (units 35) at k=51
- chr4_hap1_HG03710#1#CM086867.1: single A (beta-sat terminal) (units 40) at k=21 -> no array (units 83) at k=51
- chr4_hap2_HG03710#2#CM086881.1: single A (beta-sat terminal) (units 32) at k=21 -> no array (units 45) at k=51
- chr4_hap1_HG03742#1#CM087626.1: single B (no beta-sat terminal) (units 65) at k=21 -> no array (units 101) at k=51
- chr4_hap2_HG03742#2#CM087643.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 47) at k=51
- chr4_hap1_HG03804#1#CM088329.1: single A (beta-sat terminal) (units 19) at k=21 -> no array (units 37) at k=51
- chr4_hap2_HG03804#2#CM088378.1: single A (beta-sat terminal) (units 46) at k=21 -> no array (units 90) at k=51
- chr4_hap1_HG03831#1#CM086525.1: single B (no beta-sat terminal) (units 30) at k=21 -> no array (units 58) at k=51
- chr4_hap2_HG03831#2#CM086546.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 48) at k=51
- chr4_hap1_HG03834#1#CM085795.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 56) at k=51
- chr4_hap2_HG03834#2#CM085811.1: single A (beta-sat terminal) (units 50) at k=21 -> no array (units 93) at k=51
- chr4_hap1_HG03874#1#CM089430.1: single B (no beta-sat terminal) (units 23) at k=21 -> no array (units 28) at k=51
- chr4_hap2_HG03874#2#CM089473.1: single B (no beta-sat terminal) (units 21) at k=21 -> no array (units 45) at k=51
- chr4_hap1_HG03927#1#CM086635.1: single B (no beta-sat terminal) (units 32) at k=21 -> no array (units 63) at k=51
- chr4_hap2_HG03927#2#CM086657.1: single B (no beta-sat terminal) (units 36) at k=21 -> no array (units 75) at k=51
- chr4_hap1_HG03942#1#CM088391.1: multi-array (units 78) at k=21 -> no array (units 134) at k=51
- chr4_hap2_HG03942#2#CM088402.1: single A (beta-sat terminal) (units 38) at k=21 -> no array (units 80) at k=51
- chr4_hap1_HG04160#1#CM088425.1: single A (beta-sat terminal) (units 66) at k=21 -> no array (units 128) at k=51
- chr4_hap2_HG04160#2#CM088447.1: single A (beta-sat terminal) (units 67) at k=21 -> no array (units 125) at k=51
- chr4_hap2_HG04187#2#CM099754.1: single A (beta-sat terminal) (units 32) at k=21 -> no array (units 61) at k=51
- chr4_hap1_HG04199#1#CM086778.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 60) at k=51
- chr4_hap1_HG06807#1#CM101398.1: single A (beta-sat terminal) (units 38) at k=21 -> no array (units 75) at k=51
- chr4_hap2_HG06807#2#CM101432.1: multi-array (units 33) at k=21 -> no array (units 64) at k=51
- chr4_hap1_NA18505#1#CM089859.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 39) at k=51
- chr4_hap2_NA18505#2#CM089869.1: single B (no beta-sat terminal) (units 12) at k=21 -> no array (units 26) at k=51
- chr4_hap1_NA18508#1#CM089648.1: single B (no beta-sat terminal) (units 12) at k=21 -> no array (units 23) at k=51
- chr4_hap1_NA18522#1#CM087174.1: single A (beta-sat terminal) (units 70) at k=21 -> no array (units 148) at k=51
- chr4_hap2_NA18522#2#CM087189.1: single A (beta-sat terminal) (units 33) at k=21 -> no array (units 63) at k=51
- chr4_hap1_NA18565#1#CM094318.1: single B (no beta-sat terminal) (units 27) at k=21 -> no array (units 53) at k=51
- chr4_hap2_NA18565#2#CM094340.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 46) at k=51
- chr4_hap1_NA18570#1#CM087657.1: multi-array (units 34) at k=21 -> no array (units 66) at k=51
- chr4_hap2_NA18570#2#CM087676.1: single A (beta-sat terminal) (units 50) at k=21 -> no array (units 92) at k=51
- chr4_hap2_NA18608#2#CM089851.1: single B (no beta-sat terminal) (units 25) at k=21 -> no array (units 45) at k=51
- chr4_hap1_NA18879#1#CM099714.1: single A (beta-sat terminal) (units 84) at k=21 -> no array (units 167) at k=51
- chr4_hap2_NA18879#2#CM099729.1: single A (beta-sat terminal) (units 27) at k=21 -> no array (units 46) at k=51
- chr4_hap1_NA18906#1#CM091934.1: single A (beta-sat terminal) (units 25) at k=21 -> no array (units 49) at k=51
- chr4_hap2_NA18943#2#CM101687.1: single B (no beta-sat terminal) (units 29) at k=21 -> no array (units 58) at k=51
- chr4_hap1_NA18948#1#CM101598.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 52) at k=51
- chr4_hap2_NA18948#2#CM101587.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 51) at k=51
- chr4_hap1_NA18960#1#CM101529.1: single A (beta-sat terminal) (units 16) at k=21 -> no array (units 32) at k=51
- chr4_hap2_NA18960#2#CM101557.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 46) at k=51
- chr4_hap1_NA18967#1#CM101496.1: single B (no beta-sat terminal) (units 19) at k=21 -> no array (units 42) at k=51
- chr4_hap1_NA18971#1#CM087293.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 47) at k=51
- chr4_hap2_NA18971#2#CM087280.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 48) at k=51
- chr4_hap1_NA18974#1#CM089772.1: single B (no beta-sat terminal) (units 18) at k=21 -> no array (units 36) at k=51
- chr4_hap2_NA18974#2#CM089797.1: single B (no beta-sat terminal) (units 28) at k=21 -> no array (units 55) at k=51
- chr4_hap1_NA18982#1#CM101419.1: single B (no beta-sat terminal) (units 47) at k=21 -> no array (units 89) at k=51
- chr4_hap2_NA18982#2#CM101452.1: single A (beta-sat terminal) (units 22) at k=21 -> no array (units 33) at k=51
- chr4_hap1_NA18983#1#CM087419.1: single B (no beta-sat terminal) (units 31) at k=21 -> no array (units 62) at k=51
- chr4_hap2_NA18983#2#CM087431.1: single B (no beta-sat terminal) (units 21) at k=21 -> no array (units 38) at k=51
- chr4_hap2_NA19043#2#CM087071.1: single A (beta-sat terminal) (units 77) at k=21 -> no array (units 168) at k=51
- chr4_hap1_NA19087#1#CM087694.1: single B (no beta-sat terminal) (units 16) at k=21 -> no array (units 27) at k=51
- chr4_hap2_NA19087#2#CM087706.1: single B (no beta-sat terminal) (units 35) at k=21 -> no array (units 63) at k=51
- chr4_hap1_NA19159#1#CM099942.1: single A (beta-sat terminal) (units 68) at k=21 -> no array (units 122) at k=51
- chr4_hap2_NA19159#2#CM099961.1: single B (no beta-sat terminal) (units 13) at k=21 -> no array (units 29) at k=51
- chr4_hap1_NA19185#1#CM087722.1: single A (beta-sat terminal) (units 55) at k=21 -> no array (units 111) at k=51
- chr4_hap1_NA19338#1#CM087756.1: single A (beta-sat terminal) (units 32) at k=21 -> no array (units 67) at k=51
- chr4_hap2_NA19338#2#CM087773.1: multi-array (units 53) at k=21 -> no array (units 108) at k=51
- chr4_hap2_NA19391#2#CM087802.1: single A (beta-sat terminal) (units 28) at k=21 -> no array (units 69) at k=51
- chr4_hap1_NA19700#1#CM089619.1: single A (beta-sat terminal) (units 30) at k=21 -> no array (units 53) at k=51
- chr4_hap2_NA19700#2#CM089636.1: single A (beta-sat terminal) (units 20) at k=21 -> no array (units 48) at k=51
- chr4_hap1_NA19776#1#CM094010.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 51) at k=51
- chr4_hap2_NA19776#2#CM094047.1: single B (no beta-sat terminal) (units 24) at k=21 -> no array (units 51) at k=51
- chr4_hap1_NA19835#1#CM094024.1: single A (beta-sat terminal) (units 24) at k=21 -> no array (units 46) at k=51
- chr4_hap2_NA19909#2#CM094133.1: single A (beta-sat terminal) (units 47) at k=21 -> no array (units 90) at k=51
- chr4_hap2_NA20129#2#CM092119.1: single A (beta-sat terminal) (units 30) at k=21 -> no array (units 61) at k=51
- chr4_hap1_NA20503#1#CM102457.1: single B (no beta-sat terminal) (units 33) at k=21 -> no array (units 64) at k=51
- chr4_hap1_NA20752#1#CM099768.1: single B (no beta-sat terminal) (units 31) at k=21 -> no array (units 65) at k=51
- chr4_hap2_NA20752#2#CM099781.1: single B (no beta-sat terminal) (units 20) at k=21 -> no array (units 37) at k=51
- chr4_hap1_NA20762#1#CM102510.1: single A (beta-sat terminal) (units 13) at k=21 -> no array (units 28) at k=51
- chr4_hap1_NA20799#1#CM098686.1: single A (beta-sat terminal) (units 14) at k=21 -> no array (units 29) at k=51
- chr4_hap2_NA20799#2#CM098715.1: single A (beta-sat terminal) (units 7) at k=21 -> no array (units 25) at k=51
- chr4_hap2_NA20805#2#CM091897.1: single B (no beta-sat terminal) (units 51) at k=21 -> no array (units 94) at k=51
- chr4_hap1_NA20806#1#CM102440.1: single B (no beta-sat terminal) (units 21) at k=21 -> no array (units 42) at k=51
- chr4_hap2_NA20806#2#CM102494.1: single A (beta-sat terminal) (units 50) at k=21 -> no array (units 78) at k=51
- chr4_hap1_NA20827#1#CM102542.1: single A (beta-sat terminal) (units 45) at k=21 -> no array (units 84) at k=51
- chr4_hap2_NA20827#2#CM102562.1: single B (no beta-sat terminal) (units 10) at k=21 -> no array (units 23) at k=51
- chr4_hap1_NA20850#1#CM094291.1: single B (no beta-sat terminal) (units 50) at k=21 -> no array (units 94) at k=51
- chr4_hap2_NA20850#2#CM094305.1: single B (no beta-sat terminal) (units 21) at k=21 -> no array (units 45) at k=51
- chr4_hap1_NA21102#1#CM094425.1: single A (beta-sat terminal) (units 30) at k=21 -> no array (units 54) at k=51
- chr4_hap2_NA21102#2#CM094446.1: single B (no beta-sat terminal) (units 23) at k=21 -> no array (units 45) at k=51
- chr4_hap1_NA21110#1#CM089665.1: single B (no beta-sat terminal) (units 15) at k=21 -> no array (units 33) at k=51
- chr4_hap2_NA21110#2#CM089681.1: single B (no beta-sat terminal) (units 34) at k=21 -> no array (units 58) at k=51
- chr4_hap1_NA21144#1#CM094093.1: single B (no beta-sat terminal) (units 17) at k=21 -> no array (units 32) at k=51
- chr4_hap2_NA21144#2#CM094109.1: single A (beta-sat terminal) (units 29) at k=21 -> no array (units 54) at k=51
- chr4_hap1_NA21309#1#CM092107.1: single A (beta-sat terminal) (units 26) at k=21 -> no array (units 43) at k=51
- chr4_hap2_NA21309#2#CM092098.1: single A (beta-sat terminal) (units 34) at k=21 -> no array (units 81) at k=51

### k=31 vs k=published: 0 haplotype(s)


## The 25 published multi-array haplotypes at each k

Every haplotype here is one contig window, so both arrays sit on one contig by construction whenever num_arrays >= 2.

| haplotype | k=21 arrays/units | k=31 arrays/units | k=51 arrays/units | k=published arrays/units |
|---|---|---|---|---|
| chr10_hap2_HG00146#2#CM090038.1 | 2/26 | 2/88 | 0/100 | 2/88 |
| chr10_hap1_HG02258#1#CM085828.1 | 2/15 | 2/32 | 0/39 | 2/32 |
| chr10_hap2_HG02668#2#CM088900.1 | 3/26 | 3/59 | 0/70 | 3/59 |
| chr10_hap2_HG02717#2#CM088975.1 | 8/38 | 10/82 | 0/104 | 10/82 |
| chr10_hap1_HG02723#1#CM092035.1 | 2/13 | 2/34 | 0/40 | 2/34 |
| chr10_hap2_HG02818#2#CM091991.1 | 4/49 | 5/135 | 0/156 | 5/135 |
| chr10_hap1_HG03209#1#CM087217.1 | 2/28 | 2/51 | 0/57 | 2/51 |
| chr10_hap2_HG03470#2#CM091874.1 | 2/29 | 2/76 | 0/89 | 2/76 |
| chr10_hap1_HG03471#1#CM087526.1 | 1/16 | 2/46 | 0/53 | 2/46 |
| chr10_hap2_HG06807#2#CM101438.1 | 3/15 | 3/40 | 0/50 | 3/40 |
| chr10_hap2_NA18522#2#CM087194.1 | 2/18 | 2/49 | 0/57 | 2/49 |
| chr10_hap1_NA18565#1#CM094322.1 | 2/29 | 2/82 | 0/100 | 2/82 |
| chr10_hap1_NA19185#1#CM087728.1 | 2/54 | 2/95 | 0/118 | 2/95 |
| chr10_hap1_NA19240#1#CM099581.1 | 6/32 | 6/74 | 0/89 | 6/74 |
| chr10_hap1_NA19468#1#CM087822.1 | 6/30 | 7/65 | 0/84 | 7/65 |
| chr10_hap2_NA19468#2#CM087837.1 | 7/24 | 7/48 | 0/63 | 7/48 |
| chr4_hap2_HG00133#2#CM090063.1 | 2/41 | 2/67 | 0/73 | 2/67 |
| chr4_hap1_HG01109#1#CM092086.1 | 2/38 | 2/57 | 0/71 | 2/57 |
| chr4_hap2_HG01891#2#CM088207.1 | 2/36 | 2/47 | 0/74 | 2/47 |
| chr4_hap2_HG02257#2#CM089258.1 | 2/41 | 2/62 | 0/75 | 2/62 |
| chr4_hap1_HG02886#1#CM088084.1 | 2/30 | 2/40 | 0/59 | 2/40 |
| chr4_hap1_HG03942#1#CM088391.1 | 2/78 | 2/128 | 0/134 | 2/128 |
| chr4_hap2_HG06807#2#CM101432.1 | 2/33 | 2/45 | 0/64 | 2/45 |
| chr4_hap1_NA18570#1#CM087657.1 | 2/34 | 2/61 | 0/66 | 2/61 |
| chr4_hap2_NA19338#2#CM087773.1 | 2/53 | 2/93 | 0/108 | 2/93 |

Still multi-array at every k in {21,31,51}: 0/25.

## Mechanism: bsat fragmentation inside the published array spans, per k

Within each haplotype's published array span [bsat_start, bsat_end] (from the published caller output), the share of bases the region set labels bsat_specific and the number of bsat_specific intervals. The caller merges bsat fragments closer than BSR_MERGE_THRESHOLD = 400 bp into one BSR unit, so unit counts fall when fragments fuse.

| k | array-span bp (sum) | bsat_specific bp in spans | bsat share | bsat intervals in spans | mean bsat interval (bp) | mean gap between bsat intervals (bp) |
|---|---|---|---|---|---|---|
| 21 | 66660700 | 31756600 | 0.476 | 80387 | 395 | 415 |
| 31 | 66660700 | 22158200 | 0.332 | 89681 | 247 | 499 |
| 51 | 66660700 | 6698800 | 0.100 | 48555 | 138 | 1226 |
| published | 66660700 | 22158200 | 0.332 | 89681 | 247 | 499 |
| v2_21 | 66660700 | 50373600 | 0.756 | 2805 | 17959 | 5021 |
| v2_31 | 66660700 | 50211700 | 0.753 | 3305 | 15193 | 4163 |
| v2_51 | 66660700 | 48546100 | 0.728 | 17825 | 2723 | 754 |

## Why k = 51 returns zero arrays -- the size gate, not the exon gate

An array element must satisfy BOTH `fractions['exon'] >= EXON_DENSE_THRESHOLD` (0.5) and `size >= MIN_ARRAY_ELEMENT_SIZE` (1500) at analyze_d4z4_arrays.py:169. Measured over the first 60 chr4 haplotypes, on merged BSR clusters produced by the caller's own get_bsr_units():

| k | BSR clusters | median cluster size (bp) | share >= 1500 bp | median exon fraction | share exon >= 0.5 | array elements (both gates) |
|---|---|---|---|---|---|---|
| 21 | 1,848 | 2,700 | 84.5% | 1.00 | 88.9% | 1,506 |
| 31 | 3,174 | 800 | 47.2% | 1.00 | 84.5% | 1,466 |
| 51 | 3,550 | 200 | **0.8%** | 1.00 | 84.6% | **0** |

The exon gate is insensitive to k (84-89% of clusters pass at every value). The size gate is what fails: longer k-mers fragment the beta-satellite signal into ever smaller pieces (median merged cluster 2,700 -> 800 -> 200 bp), and at k = 51 essentially nothing reaches the hard-coded 1,500 bp floor, so no cluster is ever promoted to an array element and every haplotype reports zero arrays. This is the same constant the A6 sweep identified as the most sensitive of the caller's fourteen (invariant only over 901-1,600 bp with the published input): the k-sensitivity of the D4Z4 census runs through MIN_ARRAY_ELEMENT_SIZE, which was tuned against k = 31 data.

## HMM validation (Fig. S8B): Pearson r of main-array counts vs HMM counts, single-array chr4 haplotypes

HMM counts: HPRC-D4Z4 `Hailey-notes/HPRC_chr4_counts.xlsx` (534 contigs; supplied 2026-09-03), joined on the contig identifier exactly as plot_hmm_correlation.py does (chr4, num_arrays == 1, prefix stripped, inner join). Published: r = 0.998, n = 231.

| k | eligible (chr4, 1 array) | matched n | Pearson r | median abs diff | n with abs diff > 1 | mean (KS - HMM) |
|---|---|---|---|---|---|---|
| 21 | 239 | 231 | 0.9424 | 1.0 | 53 | -2.225 |
| 31 | 239 | 231 | 0.9977 | 1.0 | 29 | -0.797 |
| 51 | 0 | 0 | nan | nan | 0 | +nan |
| published | 239 | 231 | 0.9977 | 1.0 | 29 | -0.797 |

Cross-check of the supplied HMM table against an independent re-run of the Methods recipe on the same windows (job3_d4z4/hmm/run_hmm.sh): n = 239, Pearson r = 0.9991, identical counts for 97 of 239 contigs.

