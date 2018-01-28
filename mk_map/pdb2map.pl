#!/usr/bin/perl

#######################################
my %AA = (
ALA => 'A', ARG => 'R', ASN => 'N', ASP => 'D', CYS => 'C',
GLU => 'E', GLN => 'Q', GLY => 'G', HIS => 'H', ILE => 'I',
LEU => 'L', LYS => 'K', MET => 'M', PHE => 'F', PRO => 'P',
SER => 'S', THR => 'T', TRP => 'W', TYR => 'Y', VAL => 'V',
UNK => 'X');

my %aa=('A' => 'ALA', 'R' => 'ARG', 'N' => 'ASN', 'D' => 'ASP', 'C' => 'CYS',
        'E' => 'GLU', 'Q' => 'GLN', 'G' => 'GLY', 'H' => 'HIS', 'I' => 'ILE',
        'L' => 'LEU', 'K' => 'LYS', 'M' => 'MET', 'F' => 'PHE', 'P' => 'PRO',
        'S' => 'SER', 'T' => 'THR', 'W' => 'TRP', 'Y' => 'TYR', 'V' => 'VAL',
        'X' => 'UNK');

my %BB = ("N"=>1,"C"=>1,"O"=>1,"CA"=>1);

my %AT;
%{$AT{G}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1);
%{$AT{X}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1);
%{$AT{A}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1);
%{$AT{V}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG1"=>1,"CG2"=>1);
%{$AT{I}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG1"=>1,"CG2"=>1,"CD1"=>1);
%{$AT{L}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD1"=>1,"CD2"=>1);
%{$AT{M}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"SD"=>1,"CE"=>1);
%{$AT{F}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD1"=>1,"CD2"=>1,"CE1"=>1,"CE2"=>1,"CZ"=>1);
%{$AT{Y}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD1"=>1,"CD2"=>1,"CE1"=>1,"CE2"=>1,"CZ"=>1,"OH"=>1);
%{$AT{W}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD1"=>1,"CD2"=>1,"NE1"=>1,"CE2"=>1,"CE3"=>1,"CZ2"=>1,"CZ3"=>1,"CH2"=>1);
%{$AT{R}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD"=>1,"NE"=>1,"CZ"=>1,"NH1"=>1,"NH2"=>1);
%{$AT{H}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"ND1"=>1,"CD2"=>1,"CE1"=>1,"NE2"=>1);
%{$AT{K}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD"=>1,"CE"=>1,"NZ"=>1);
%{$AT{D}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"OD1"=>1,"OD2"=>1);
%{$AT{E}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD"=>1,"OE1"=>1,"OE2"=>1);
%{$AT{N}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"OD1"=>1,"ND2"=>1);
%{$AT{Q}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD"=>1,"OE1"=>1,"NE2"=>1);
%{$AT{C}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"SG"=>1);
%{$AT{S}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"OG"=>1);
%{$AT{T}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"OG1"=>1,"CG2"=>1);
%{$AT{P}} = ("OXT"=>1,"N"=>1,"C"=>1,"O"=>1,"CA"=>1,"CB"=>1,"CG"=>1,"CD"=>1);
#######################################


my ($pdb,$chk,$ss,$map,$new_pdb);
while ($arg = shift())
{
	if ($arg =~ m/^-/)
	{
		$arg = lc($arg);
		if ($arg eq "-pdb") {$pdb = shift();next;}
		if ($arg eq "-chk") {$chk = shift();next;}
		if ($arg eq "-ss" ) {$ss  = shift();next;}
		if ($arg eq "-map") {$map = shift();next;}
		if ($arg eq "-new_pdb") {$new_pdb = shift();next;}
	}
}

unless(-e $pdb){die("-pdb '$pdb' not found");}
unless(defined $map){die("-map '$map' not defined");}
unless(defined $new_pdb){die("-new_pdb '$new_pdb' not defined");}

my @MAP;
my @NEW_PDB;
my $errors = 0;

my ($start,$pdb_seq,$ref_lines,$ref_CON) = get_contacts($pdb);
my %lines = %$ref_lines;
my @CON = @$ref_CON;

if(!-e $chk){
	my $len = length($pdb_seq);
	push(@MAP,"LEN $len");
	my $n = 0;
	while(exists $CON[$n]){
		my $i = $CON[$n][2] - $start;
		my $j = $CON[$n][3] - $start;
		push(@MAP,"CON $i $j 1");
		$n++;
	}
	my $resn = $start;
	for(my $i = 0; $i < length($pdb_seq); $i++){
		if(substr($pdb_seq,$i,1) ne "-"){
			for my $line (@{$lines{$resn}}){
				push(@NEW_PDB,substr($line,0,22).sprintf("%4s",$resn-$start).substr($line,26));
			}
		}
		$resn++;
	}
}
else{
	my ($chk_seq,@CHK_cs) = parse_checkpoint_file($chk);
	my $pdb_seq_tmp = $pdb_seq;
	$pdb_seq_tmp =~ s/-//g;
	if($pdb_seq_tmp eq $chk_seq){
		my %SS = get_stride($ss,"RENUM");
		my $len = length($pdb_seq_tmp);
		push(@MAP,"LEN $len");
		my $n = 0;
		while(exists $CON[$n]){
			my $i = $CON[$n][0];
			my $j = $CON[$n][1];
			push(@MAP,"CON $i $j 1");
			$n++;
		}
		for(my $n = 0; $n < $len; $n++){
			my @fix_cs;for my $i (0..19){$fix_cs[$i] = sprintf("%.6f",$CHK_cs[$n][$i]);}

			my $aa = substr($chk_seq,$n,1);
			my $ss = "X";
			if(exists $SS{$n+1}){
				if($SS{$n+1}{"AA"} ne $aa){
					print STDERR "ERROR: chk[$n]='$aa' != ss[$n]='".$SS{$n+1}{"AA"}."'\n";
					$errors++;
				}				
				$ss = $SS{$n+1}{"SS"};
			}

			push(@MAP,"PRF $n $aa $ss @fix_cs");
		}
		my $resn = $start;
		for(my $i = 0; $i < length($pdb_seq); $i++){
			if(substr($pdb_seq,$i,1) ne "-"){
				for my $line (@{$lines{$resn}}){
					push(@NEW_PDB,substr($line,0,22).sprintf("%4s",$i).substr($line,26));
				}
			}
			$resn++;
		}
	}
	else{
		# if pdb sequence and chk sequence do NOT match
		# check if they match when ignoring "missing density" positions
		my $ok = 1;
		my $j = 0;
		for(my $i = ($start-1); $i < length($chk_seq); $i++){
			if($j < length($pdb_seq)){
				my $c = substr($chk_seq,$i,1);
				my $p = substr($pdb_seq,$j,1);
				if($p ne "-" and $p ne $c){
					$ok = 0;
					print STDERR "ERROR: chk[$i]='$c' != pdb[$j]='$p'\n";
					$errors++;
				}
			}
			$j++;
		}
		if($ok == 1){
			my %SS = get_stride($ss);

			my $len = length($chk_seq);
			push(@MAP,"LEN $len");
			my $n = 0;
			while(exists $CON[$n]){
				my $i = $CON[$n][2] - 1;
				my $j = $CON[$n][3] - 1;
				push(@MAP,"CON $i $j 1");
				$n++;
			}
			for(my $n = 0; $n < $len; $n++){
				my @fix_cs;for my $i (0..19){$fix_cs[$i] = sprintf("%.6f",$CHK_cs[$n][$i]);}
				my $aa = substr($chk_seq,$n,1);
				my $ss = "X";
				if(exists $SS{$n+1}){
					if($SS{$n+1}{"AA"} ne $aa){
						print STDERR "ERROR: chk[$n]='$aa' != ss[$n]='".$SS{$n+1}{"AA"}."'\n";
						$errors++;
					}				
					$ss = $SS{$n+1}{"SS"};
				}
				push(@MAP,"PRF $n $aa $ss @fix_cs");
			}
			my $resn = $start;
			for(my $i = 0; $i < length($pdb_seq); $i++){
				if(substr($pdb_seq,$i,1) ne "-"){
					for my $line (@{$lines{$resn}}){
						push(@NEW_PDB,substr($line,0,22).sprintf("%4s",$resn - 1).substr($line,26));
					}
				}
				$resn++;
			}
		}
	}
}
if($errors == 0){
	open(MAP,">$map");
	for my $line (@MAP){
		print MAP "$line\n";
	}
	close(MAP);
	open(NEW_PDB,">$new_pdb");
	for my $line (@NEW_PDB){
		print NEW_PDB "$line\n";
	}
	close(NEW_PDB);
}
sub parse_checkpoint_file
{
	my $filename = shift;
	my $buf;
	my @aa_order = split( //, 'ACDEFGHIKLMNPQRSTVWY' );
	my @altschul_mapping = ( 0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18 );
	my @output;
	open(INPUT,$filename) or die("Couldn't open $filename for reading.\n");
	read(INPUT,$buf,4) or die("Couldn't read $filename!\n");
	my $seqlen = unpack("i",$buf);
	read(INPUT,$buf,$seqlen) or die("Premature end: $filename.\n");
	my $seqstr = unpack("a$seqlen",$buf);
	for (my $i = 0 ; $i < $seqlen ; ++$i ) {
	        read(INPUT,$buf,160) or die("Premature end: $filename, line: $i\n");
        	my @w = unpack("d20",$buf);
		for (my $j = 0 ; $j < 20 ; ++$j ) {
			$output[$i][$j] = $w[$altschul_mapping[$j]];
		}
	}
	close(INPUT);
	return($seqstr,@output);
}
sub distance{return sqrt((($_[0]-$_[3])**2)+(($_[1]-$_[4])**2)+(($_[2]-$_[5])**2));}
sub get_contacts
{

	my ($pdb) = @_;

	my $seqsep_cutoff = 3;
	my $ca_cutoff = 20;
	my $min_cutoff = 5;

	my @CA;
	my %XYZ;
	my %SEQ;
	my %PDB_LINES;
	open(PDB,$pdb);
	while(my $line = <PDB>){
		chomp($line);
		if(substr($line,0,6) eq "HETATM" and substr($line,17,3) eq "MSE"){
			$line =~ s/^HETATM/ATOM  /;
			substr($line,17,3,"MET");
			substr($line,12,2,"SD") if substr($line,12,2) eq "SE";
		}
		if(substr($line,0,4) eq "ATOM"){
			my $aa = substr($line,17,3);
			if(exists $AA{$aa}){
				my $resn = int(substr($line,22,5));
				my $resi = $AA{$aa};
				$SEQ{$resn} = $resi;
				my $atom = substr($line,12,4);
				$atom =~ s/ //g;
				if($atom eq "CA"){push(@CA,$resn);}
				if(exists $AT{$resi}{$atom}){
					my $x = substr($line,30,8)+0;
					my $y = substr($line,38,8)+0;
					my $z = substr($line,46,8)+0;
					@{$XYZ{$resn}{$atom}} = ($x,$y,$z);
				}
				if(exists $AT{"X"}{$atom}){
					push(@{$PDB_LINES{$resn}},$line);
				}
			}
		}
	}
	close(PDB);
	my @CON;
	my $i = 0;
	while(exists $CA[$i]){
		my $j = $i+1;
		while(exists $CA[$j]){
			my $sep = abs($CA[$i]-$CA[$j]);
			if($sep >= $seqsep_cutoff){
				if(distance(@{$XYZ{$CA[$i]}{"CA"}},@{$XYZ{$CA[$j]}{"CA"}}) <= $ca_cutoff){
					my $min_dist;
					for my $ai (keys %{$XYZ{$CA[$i]}}){
						for my $aj (keys %{$XYZ{$CA[$j]}}){
							my $dist = sprintf("%.1f",distance(@{$XYZ{$CA[$i]}{$ai}},@{$XYZ{$CA[$j]}{$aj}}));
							if(!defined $min_dist or $dist < $min_dist){$min_dist = $dist;}
						}
					}
					if($min_dist <= $min_cutoff){
						push(@CON,([$i,$j,$CA[$i],$CA[$j]]));
					}
				}
			}
			$j++;
		}
		$i++;
	}
	my $seqstr;
	for my $i ($CA[0]..$CA[-1]){
		if(exists $SEQ{$i}){$seqstr .= $SEQ{$i};}
		else{$seqstr .= "-";}
	}
	return($CA[0],$seqstr,\%PDB_LINES,\@CON);
}
sub get_stride
{
	my ($pdb,$mode) = @_;
	my %stride_convert = ('H' => 'H', 'G' => 'H', 'I' => 'H','E' => 'E', 'B' => 'E', 'b' => 'E', 'T' => 'C', 'C' => 'C');
	#	   H	    Alpha helix
	#	   G	    3-10 helix
	#	   I	    PI-helix
	#	   E	    Extended conformation
	#	   B or	b   Isolated bridge
	#	   T	    Turn
	#	   C	    Coil (none of the above)

	my $seqstr;
	my %DATA;
	open(SS,$pdb);
	while(my $line = <SS>){
		chomp($line);
		if(substr($line,0,3) eq "ASG"){
			#REM  |---Residue---|    |--Structure--|   |-Phi-|   |-Psi-|  |-Area-|      ~~~~
			#ASG  SER A    1    1    C          Coil    360.00    175.11      78.8      ~~~~
			my ($tag,$resi,$chain,$resn,$num,$SS,$SS_verbose,$phi,$psi,$area) = split(/\s+/,$line);

			if($mode eq "RENUM"){$resn = $num;}

			$DATA{$resn}{"AA"} = $AA{$resi};
			$seqstr .= $AA{$resi};
			$DATA{$resn}{"SS"} = $stride_convert{$SS};
			$DATA{$resn}{"phi"} = $phi;
			$DATA{$resn}{"psi"} = $psi;
			$DATA{$resn}{"area"} = $area;
		}
	}
	close(SS);
	return(%DATA);
}
##########################
