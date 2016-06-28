# _map_align_
_map_align_ takes two contact maps and returns an alignment that attempts to maximize the number of overlapping contacts while minimizing the number of gaps.


![example image](https://raw.githubusercontent.com/sokrypton/map_align/master/map_align_fig.png)

### Installation
```sh
$ g++ -O3 -std=c++0x -o map_align main.cpp
```

### Usage
```
-------------------------------------------------------------------
MAP_ALIGN
-------------------------------------------------------------------
-a             contact map A               [REQUIRED]
-b             contact map B               [REQUIRED]
-gap_o         gap opening penalty         [Default=-1]
-gap_e         gap extension penalty       [Default=-0.01]
-sep_cut       seq seperation cutoff       [Default=3]
-iter          number of iterations        [Default=20]
-silent
-------------------------------------------------------------------
Experimental features
-------------------------------------------------------------------
-use_gap_ss    penalize gaps at secondary structure elements(SSE)
-gap_ss_w      gap penality weight at SSE  [Default=2]
-use_prf       use sequence profile
-prf_w         profile weight              [Default=1]
-------------------------------------------------------------------
```

### Experimental features
- We are still testing these features and there might be bugs!
- "gap_ss" uses the SS info provided in the "PRF" line to increasing gap penalities within secondary structure elements, favoring gaps in loop or regions of missing density. 
- "prf" uses the frequencies provided in the "PRF" line to assist in alignment. This option was intented to help align regions void of contact information.  The average frequencies of input profiles (from both maps) is used to compute the background frequencies. WARNING: This option may hurt finding the optimial alignment when aligning non-homologous proteins that share the same fold due to convergent evolution.

```sh
$ map_align -a A.map -b B.map
```

### contact map format
- ```LEN 440``` - [len]gth
- ```CON 0 4 1.0```  - [con]tact, i, j and weight.
- ```PRF 0 A H 0.01 ... 0.01``` (optional) profile, i, amino acid (AA), secondary structure (SS) and profile frequencies (20 values). The order of the frequencies does not matter, as long as they match between the two contact maps being compared. H = Helix; E = Sheet; all other characters treated equally. 
