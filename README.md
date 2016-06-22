# _map_align_
_map_align_ takes two contact maps and returns an alignment that tries to maximize the number of overlapping contacts while minimizing the number of gaps.

### Installation
```sh
$ g++ -O3 -std=c++0x -o map_align main.cpp
```

### Usage
```
-----------------------------------------------------
MAP_ALIGN                      
-----------------------------------------------------
-a      contact map A             [REQUIRED]
-b      contact map B             [REQUIRED]
-prf    use sequence profile      [Default=0]
-gap_o  gap opening penalty       [Default=-1]
-gap_e  gap extension penalty     [Default=-0.01]
-ss_cut seq seperation cutoff     [Default=3]
-----------------------------------------------------
```
```sh
$ map_align -a A.map -b B.map
```

### contact map format
- ```LEN 440``` - [len]gth
- ```CON 0 4 1.0```  - [con]tact, i, j and weight.
- ```PRF 0 A H 0.01 ... 0.01``` (optional) profile, i, amino acid (AA), secondary structure (SS) and profile frequencies (20 values). The order of the frequencies does not matter, as long as they match between the two contact maps being compared. The SS is not currently used in map_align and can be "X".