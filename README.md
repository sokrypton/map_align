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
Advanced options
-------------------------------------------------------------------
-range_a       trim map A to specified range(s) (eg. 0-20 50-100)
-range_b       trim map B to specified range(s)
-------------------------------------------------------------------
Experimental features
-------------------------------------------------------------------
-use_gap_ss    penalize gaps at secondary structure elements(SSE)
-gap_ss_w      gap penality weight at SSE  [Default=2]
-use_prf       use sequence profile
-prf_w         profile weight              [Default=1]
-------------------------------------------------------------------
```
```sh
$ map_align -a A.map -b B.map
```

### contact map format
- ```LEN 440``` - [len]gth
- ```CON 0 4 1.0```  - [con]tact, i, j and weight.
- ```PRF 0 A H 0.01 ... 0.01``` (optional) profile, i, amino acid (AA), secondary structure (SS) and profile frequencies (20 values). The order of the frequencies does not matter, as long as they match between the two contact maps being compared. H = Helix; E = Sheet; all other characters treated equally.

### parsing output
*  the output will be a single line (if -silent is used).
   * ```MAX params map_a map_b contact_sco gap_sco total_sco 0:0 1:1 2:2 ...```
   * the alignment is provided as ```0:0``` with index of first and second map.
   * if "-use_prf" flag is used, the output will include an extra profile_sco column:
      * ```MAX params map_a map_b contact_sco gap_sco profile_sco total_sco 0:0 1:1 2:2 ...```
      
### Experimental features
- WARNING: these are experimental features and may not work correctly!
- "gap_ss" uses the SS info provided in the "PRF" line to increasing gap penalities within secondary structure elements, favoring gaps in loop or regions of missing density. 
- "prf" uses the frequencies provided in the "PRF" line to assist in alignment. This option was intented to help align regions void of contact information.  The average frequencies of input profiles (from both maps) is used to compute the background frequencies. WARNING: This option may hurt finding the optimial alignment when aligning non-homologous proteins that share the same fold due to convergent evolution.

### Convert GREMLIN/CCMPRED results to .map files
* To do this you'll need the following files (see "mk_map" directory).
* aln file: containing the alignment used as input to GREMLIN/CCMPRED (after gap removal).
* cut file: One line containing the full sequence, the second line containing the trimmed sequence (with "-" to indicated positions removed). This file is used to determine the mapping from the matrix file to the full length sequence.
* mtx file: symmetric matrix containing length x length values of coupling results.
* chk file: for profile generation, binary file from [csbuild](https://github.com/cangermueller/csblast). To generate this file you'll need the a3m/fas file before gap removal. This alignment is then given as input to csbuild. [csbuild](https://github.com/cangermueller/csblast) is used to add context-specific pseudocounts to the profiles.
   * ```csbuild -i TMP.a3m -I a3m -o TMP.chk -O chk -D csblast-2.2.3/data/K4000.crf```
* A perl script is provided to convert the data to a contact map file: 
   * ```perl mk_map.pl -aln TMP.aln -cut TMP.cut -mtx TMP.mtx -chk TMP.chk -map TMP.map -do_apc``` 
   * the "-do_apc" flag is REQUIRED if no APC (Average Product Correction) was performed to input mtx.


### (work in progress) Convert PDB to .map files 
* pdb file: single chain, renumbered to match
* chk file: (see above for more information)
   * ```hhblits -o /dev/null -d DATABASE -cpu ??? -id 90 -cov 75 -e 1E-10 -n 2 -i TMP.fasta -oa3m TMP.a3m```
   * ```csbuild -i TMP.a3m -I a3m -o TMP.chk -O chk -D csblast-2.2.3/data/K4000.crf```
* ss file: output from [stride](http://webclu.bio.wzw.tum.de/stride/) (secondary structure)
   * ```stride TMP.pdb > TMP.ss```
* A perl script is provided to convert the data to a contact map file: 
   * ```perl pdb2map.pl -stride TMP.ss -pdb TMP.pdb -chk TMP.chk -map TMP.map```
