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

### Convert PDB to .map files 
* see "cmd_pdb2map" file for full example!
* pdb file: single chain, single model, no HETATM
* chk file: (see above for more information)
   * WARNING: the PDB numbering should match the the input FASTA file
   * ```hhblits -o /dev/null -d DATABASE -cpu ??? -id 90 -cov 75 -e 1E-10 -n 2 -i TMP.fasta -oa3m TMP.a3m```
   * ```csbuild -i TMP.a3m -I a3m -o TMP.chk -O chk -D csblast-2.2.3/data/K4000.crf```
* ss file: output from [stride](http://webclu.bio.wzw.tum.de/stride/) (secondary structure)
   * ```stride TMP.pdb > TMP.ss```
   * WARNING: stride only reads ATOM records!
* A perl script is provided to convert the data to a contact map file: 
   * ```perl pdb2map.pl -pdb TMP.pdb -map TMP.map -chk TMP.chk -ss TMP.ss ```
   * ss and chk files are NOT required

