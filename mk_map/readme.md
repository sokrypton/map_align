### Convert GREMLIN/CCMPRED results to .map files
* see "cmd" file for full example!
* aln file: containing the alignment used as input to GREMLIN/CCMPRED (after gap removal).
* cut file: One line containing the full sequence, the second line containing the trimmed sequence (with "-" to indicated positions removed). This file is used to determine the mapping from the matrix file to the full length sequence.
* mtx file: symmetric matrix containing length x length values of coupling results.
* chk file: (NOT required) for profile generation, binary file from [csbuild](https://github.com/cangermueller/csblast). To generate this file you'll need the a3m/fas file before gap removal. This alignment is then given as input to csbuild. [csbuild](https://github.com/cangermueller/csblast) is used to add context-specific pseudocounts to the profiles.
   * ```csbuild -i TMP.a3m -I a3m -o TMP.chk -O chk -D csblast-2.2.3/data/K4000.crf```
* A perl script is provided to convert the data to a contact map file: 
   * ```perl mk_map.pl -aln TMP.aln -cut TMP.cut -mtx TMP.mtx -chk TMP.chk -map TMP.map -do_apc``` 
   * the "-do_apc" flag is REQUIRED if no APC (Average Product Correction) was performed to input mtx.
