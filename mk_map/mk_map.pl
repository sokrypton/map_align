#!/usr/bin/perl

my ($aln,$cut,$mtx,$chk,$map);
my $do_apc = 0;
while ($arg = shift())
{
	if ($arg =~ m/^-/)
	{
		$arg = lc($arg);
		if ($arg eq "-aln") {$aln = shift();next; }
		if ($arg eq "-cut") {$cut = shift();next; }
		if ($arg eq "-mtx") {$mtx = shift();next; }
		if ($arg eq "-chk") {$chk = shift();next; }
		if ($arg eq "-map") {$map = shift();next; }
		if ($arg eq "-do_apc") {$do_apc = 1; next; }
	}
}
unless(-e $aln){die("-aln '$aln' not found");}
unless(-e $cut){die("-cut '$cut' not found");}
unless(-e $mtx){die("-mtx '$mtx' not found");}
unless(defined $map){die("-map '$map' not defined");}

# params for prob calculation (24Aug2015)
# see: https://elifesciences.org/content/4/e09248/figure16
my $N_v_Bx = -0.53;
my $N_v_By = 5.46;
my $N_v_Sx = 0.50;
my $N_v_Sy = 0.58;
my $N_v_H = 1;
my $N_v_Y = 0;

# get sequence
open(CUT,$cut);
my $seq = <CUT>;chomp($seq);
close(CUT);
my $len = length($seq);

# for P(contact) calculation
my $nf = get_nf($aln);
my $N_v_B = $N_v_By * pow($nf,$N_v_Bx);
my $N_v_S = $N_v_Sy * pow($nf,$N_v_Sx);

open(MAP,">$map");
print MAP "LEN $len\n";
# get contacts
my @CST = parse_mtx($cut, $mtx, $do_apc, 3, 1.5);
my $c = 0;
while(exists $CST[$c]){
	my ($i,$j,$sc) = @{$CST[$c]};
	my $prob = sprintf("%.3f",$N_v_H/(1+exp(-$N_v_S * ($sc - $N_v_B))) + $N_v_Y);	
	print MAP "CON $i $j $prob\n";
	$c++;
}
if(-e $chk){
	# get profile
	my @CHK_cs = parse_checkpoint_file($chk);
	my @SEQ = split(//,$seq);
	my $n = 0;
	while(exists $SEQ[$n]){
		my @fix_cs;for my $i (0..19){$fix_cs[$i] = sprintf("%.6f",$CHK_cs[$n][$i]);}
		print MAP "PRF $n ".$SEQ[$n]." X @fix_cs\n";
		$n++;
	}
}
close(MAP);

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
	return(@output);
}
sub mean{my($data) = @_;if(not @$data){return 0;}else{my $total = 0;foreach (@$data){$total += $_;}my $average = $total / @$data;return $average;}}
sub pow{return($_[0] ** $_[1]);}
sub parse_mtx
{
	my ($cut, $mtx, $do_apc, $cst_sep, $cst_num) = @_;
	# read cut file
	my @FULL;my @CUT;
	open(CUT,$cut);my $n = 0;
	while(my $line = <CUT>){
		chomp($line);
		if($n == 0){@FULL = split (//,$line);}
		if($n == 1){@CUT = split (//,$line);}
		$n++;
	}
	close(CUT);
	if(scalar(@FULL) != scalar(@CUT)){die("Error with -cut '$cut' file");}
	# save cut mapping
	my @f2c;my @c2f;my @f2a;
	my $f = 0;my $c = 0;
	while(exists $FULL[$f]){
		$f2a[$f] = $FULL[$f];
		if($CUT[$f] ne "-"){$f2c[$f] = $c;$c2f[$c] = $f;$c++;}
		$f++;
	}
	# read mtx file
	my @vals; my $v = 0;my @MTX;my $n = 0;
	open(MAT,$mtx);
	while(my $line = <MAT>){
		chomp($line);
		if(substr($line,0,1) ne "#"){
			my @tmp = split(/\s+/,$line);my $nn = 0;
			foreach $t (@tmp){if($nn > $n){push(@vals,$t);}$MTX[$n][$nn] = $t;$nn++;}
			$n++;
		}
	}
	close(MAT);
	my @SCO;my @vSCO;
	if($do_apc == 1)
	{
		# save mean
		my @row_mean;
		my $i = 0;
		while(exists $MTX[$i]){
			my @rows;
			my $j = 0;my $r = 0;
			while(exists $MTX[$i][$j]){
				if($i != $j){push(@rows,$MTX[$i][$j]);}
				$j++;
			}
			$row_mean[$i] = mean(\@rows);
			$i++;
		}
		my $all_mean = mean(\@vals);
		# do APC, get sub-matrix with seq_sep >=3
		my $i = 0;
		while(exists $MTX[$i]){
			my $j = $i + 1;
			while(exists $MTX[$i][$j]){
				my $apc = ($row_mean[$i]*$row_mean[$j])/$all_mean;
				my $sco = $MTX[$i][$j] - $apc;
				push(@SCO,([$sco,$c2f[$i],$c2f[$j]]));
				if(abs($c2f[$j]-$c2f[$i]) >= $cst_sep){push(@vSCO,$sco);}
				$j++;
			}
			$i++;
		}
	}
	else
	{
		my $i = 0;
		while(exists $MTX[$i]){
			my $j = $i + 1;
			while(exists $MTX[$i][$j]){
				my $sco = $MTX[$i][$j];
				push(@SCO,([$sco,$c2f[$i],$c2f[$j]]));
				if(abs($c2f[$j]-$c2f[$i]) >= $cst_sep){push(@vSCO,$sco);}
				$j++;
			}
			$i++;
		}
	}
	@vSCO = sort {$b <=> $a} @vSCO;

	# get top 1.5L preds
	my $preds = sprintf("%.0f", ($c * $cst_num));
	my $sco_tot = 0;
	for my $n (0..($preds-1)){$sco_tot += $vSCO[$n];}

	my $sco_avg = $sco_tot/$preds;
	my $sco_min = $vSCO[($preds-1)];

	@SCO = sort {$b->[0] <=> $a->[0]} @SCO;
	my $p = 0;
	my $p_ss = 0;
	my @CST;
	while(exists $SCO[$p] and $p_ss < $preds) # going through stats array
	{
		my ($sco,$i,$j) = @{$SCO[$p]};
		my $ss = abs($j-$i);
		if($ss >= $cst_sep)
		{
			my $s_sco = (0.5 * ($sco - $sco_min))/($sco_avg - $sco_min) + 0.5;
			push(@CST,([$i,$j,$s_sco]));
			$p_ss++;
		}
		$p++;
	}
	return(@CST);
}
sub get_nf
{
	my $msa = $_[0];
	my @SEQ;
	my $n = 0;
	open(MSA,$msa);
	while(my $line = <MSA>)
	{
		chomp($line);
		my @tm = split(/\s+/,$line);
		if(exists $tm[1]){push(@SEQ,$tm[1]);}
		else{push(@SEQ,$tm[0]);}
		$n++;
	}
	close(MSA);
	my $len = length($SEQ[0]);
	my $chk = $len * 0.8;
	my @N;
	my $n_80 = 0;
	for my $i (0..($n-1))
	{
		my $w_80 = $N[$i];
		for my $j (($i+1)..($n-1))
		{
			my $ham = (($SEQ[$i] ^ $SEQ[$j]) =~ tr/\0//);
			if($ham > $chk)
			{
				$N[$j] += 1;
				$w_80 += 1;
			}
			$j++;
		}
		$n_80 += 1/($w_80+1);
		$i++;
	}
	return($n_80/sqrt($len));
}
