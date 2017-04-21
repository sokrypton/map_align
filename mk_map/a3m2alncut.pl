#!/usr/bin/perl

# this script takes an a3m or fas alignment file and remove positions > 50% gaps and sequences < 75% coverage of query sequence
# it returns aln "aln" file with one sequence per line
# it also returns an "cut" file that reports the regions that were cut

if(scalar(@ARGV) < 3){
	die("perl a3m2alncut.pl a3m/fas aln cut");
}
cut(@ARGV);
sub cut
{
        my ($a3m,$aln,$cut) = @_;
        my $ratio = 0.5;
        my @col;
        my $n = 0;
        my @seq;my $s = 0;
        open(A3M,$a3m);
        while(my $line = <A3M>)
        {
                chomp($line);
                if(length($line) > 0)
                {
                        if(substr($line,0,1) eq ">")
                        {
                                if(exists $seq[$s]){$s++;}
                                $seq[$s][0] = $line;
                        }
                        else
                        {
                                $line =~ s/[a-z]//g; #make fasta
                                $seq[$s][1] .= $line;
                        }
                }
        }
        close(A3M);
        my $rlen = length($seq[0][1]);
        my $n = 0;
        my @col_score;
        my $col_score_tot = 0;
        while(exists $seq[$n])
        {
                my $gap = int($seq[$n][1] =~ tr/-//);
                my $non_gap = ($rlen-$gap)/$rlen;

                my $m = 0;
                for my $p (split(//,$seq[$n][1]))
                {
                        if($p ne "-"){$col_score[$m] += $non_gap;}
                        $col[$m] .= $p;
                        $m++;
                }
                $col_score_tot += $non_gap;
                $n++;
        }
        my $m = 0;
        my $row_score_tot = 0;
        while(exists $col[$m])
        {
                my $clen = length($col[$m]);
                my $gap = int($col[$m] =~ tr/-//);
                my $non_gap = ($clen-$gap)/$clen;

                $row_score_tot += $non_gap;
                $m++;
        }
        my $cut_seq;
        my $cut_len = 0;
        my @col_check;
        my $m = 0;
        while(exists $col[$m])
        {
                my $clen = length($col[$m]);
                my $gap = int($col[$m] =~ tr/-//);
                my $non_gap = ($clen-$gap);

                my ($val,$cutoff);
                my $val = $non_gap/$col_score_tot;
                my $cutoff = $ratio;

                if($val >= $cutoff){$col_check[$m] = 1;}else{$col_check[$m] = 0;}

                if($col_check[$m] == 0){$cut_seq .= "-";}
                else{$cut_seq .= substr($seq[0][1],$m,1); $cut_len++;}
                $m++;
        }
        open(CUT,">$cut");
        print CUT $seq[0][1]."\n";
        print CUT $cut_seq."\n";
        close(CUT);
        open(ALN,">$aln");
        my $n = 0;
        my $N = 0;
        while(exists $seq[$n])
        {
                my $m = 0;
                my $new_seq;
                while($m < $rlen)
                {
                        if($col_check[$m] == 1)
                        {
                                my $p = substr($seq[$n][1],$m,1);
                                $new_seq .= $p;
                        }
                        $m++;
                }
                my $gap = int($new_seq =~ tr/-//);
                if(($gap/$cut_len) < 0.25)
                {
                        print ALN $new_seq."\n";
                        $N++;
                }
                $n++;
        }
        close(ALN);
}
