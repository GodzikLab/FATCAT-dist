#!/usr/bin/perl -w 
use strict;

# Zhanwen Li @ 12/10/2008
# Similar to FATCATDB but no calculation in this script

my $me = "$0";
my $script = $0; $script =~ s/^.*\///;
my $usage = 
    "$script [-a alignment-file]/[-r report-file] [-f function] [-o output-file]\n".
    "  Note: both alignment-file or report-file are allowed as input\n".
    "        alignment-file: the detailed alignment result of the database searching\n".
    "        report-file: the one line one pair format of the searching result\n".
    "  Possible function:\n".
    "        -f report      (report the alignment results in a one line one pair format in text, default)\n".
    "        -f extalign    (extract the alignment results for given protein pairs or by alignment features such as gap, rmsd, pval etc)\n".
    "        -f extreport   (extract reports for given protein pairs or by alignment features such as gap, rmsd, pval etc)\n".
    "        -f id          (extract search hit protein IDs)\n".
    "  The flowing options are for -f report or -f extalign or -f extreport:\n".
    "        -g gap       (gap threshold for output list)\n".
    "        -t twist     (twist threshold for output list)\n".
    "        -alilen len  (align-length threshold for output list)\n".
    "        -equlen len  (threshold for length of equivalent positions)\n".
    "        -rmsd num    (rmsd threshold for output list)\n".
    "        -pval num    (Pvalue threshold for output list)\n".
    "        -n num       (the first n records of the given input report file or alignment file)\n".
    "  -pair              (a list of protein pairs to extract, works with -f extreport/extalign)\n".
    "  -h yes/no          (print header or not, default yes)\n";         

die "$usage" unless @ARGV>0;

my %para_h = @ARGV;
my $ali_f; if (!($ali_f = delete $para_h{"-a"})) { $ali_f = ""; }
my $report_f; if (!($report_f = delete $para_h{"-r"})) { $report_f = ""; }
my $func; if (!($func = delete $para_h{"-f"})) { die "Please enter a function option\n\n $usage"; }
my $gap; if (!($gap = delete $para_h{"-g"})) { $gap = ""; }
my $twist; if (!($twist = delete $para_h{"-t"})) { $twist = ""; }
my $alilen; if (!($alilen = delete $para_h{"-alilen"})) { $alilen = ""; }
my $equlen; if (!($equlen = delete $para_h{"-equlen"})) { $equlen = ""; }
my $rmsd; if (!($rmsd = delete $para_h{"-rmsd"})) { $rmsd = ""; }
my $pval; if (!($pval = delete $para_h{"-pval"})) { $pval = ""; }
my $num; if (exists($para_h{"-num"})) { $num = $para_h{"-num"}; }
my $pair_f; if (!($pair_f = delete $para_h{"-pair"})) { $pair_f = ""; }
my $out_f; if (!($out_f = delete $para_h{"-o"})) { die "Please enter output file\n\n $usage"; }
my $header; if (!($header = delete $para_h{"-h"})) { $header = "yes"; }

if (!$ali_f  &&  !$report_f) { die "Please enter input file: alignment file or report file\n"; }
if ($ali_f && $report_f) { print "Warning: only alignment file is used for parsing\n"; }

if ($func =~ /^report$/i) {
    printReport($ali_f, $gap, $twist, $alilen, $equlen, $rmsd, $pval, $out_f);
}  elsif ($func =~ /^extalign$/i) {
    extactAlign($ali_f, $pair_f,  $gap, $twist, $alilen, $equlen, $rmsd, $pval, $out_f);
} elsif ($func =~ /^extreport$/i) {
    extractReport($ali_f, $report_f, $pair_f, $gap, $twist, $alilen, $equlen, $rmsd, $pval, $out_f);
} elsif ($func =~ /^id/i) {
    extHitId($ali_f, $report_f,  $gap, $twist, $alilen, $equlen, $rmsd, $pval, $num, $out_f);
}


sub extHitId {
    my ($ali_f, $report_f,  $gap, $twist, $alilen, $equlen, $rmsd, $pval, $num, $out_f) = @_;

    if (!$report_f || !(-e "$report_f") || (-z "$report_f")) {
	if ($ali_f && (-e "$ali_f") && !(-z "$ali_f")) {
	    $report_f = "___temp_report___";
	    printReport($ali_f, $gap, $twist, $alilen, $equlen, $rmsd, $pval, $report_f);
	} 
    }
    if (!open(REP, "<$report_f") ) { return; }
    if (!open(OUT, ">$out_f")) { close(REP); return; }
    my $c = 0;
    while (<REP>) {
	if (!/\S/ || /^\#/) { next; }
	$c++;
	if (defined($num) && $num > 0 && $c >= $num) { last; }
	my ($id1, $id2) = split;
	print OUT "$id2\n";
    }
    close(OUT);
    close(REP);
}
	

sub extractReport {
    my ($ali_f, $report_f, $pair_f, $gap, $twist, $alilen, $equlen, $rmsd, $pval, $out_f) = @_;

    if ($ali_f && (-e "$ali_f") && !(-z "$ali_f")) {
	my $temp_f = "__aln4pairs__";
	extactAlign($ali_f, $pair_f,  $gap, $twist, $alilen, $equlen, $rmsd, $pval, $temp_f);
	system("$me -a $temp_f -f report -o $out_f");
	system("rm $temp_f");
    } elsif ($report_f) {

	my %pairs_h = readPair($pair_f);

	open(REP, "<$report_f") || die "can't open $report_f\n"; 
	open(OUT, ">$out_f") || die "Can't open file $out_f\n";
	if ($header =~ /YES/i) {
	    print OUT "# pdb1 pdb1 len1 len2 twist iniLen iniRmsd optLen optRmsd chainRms score alnLen gap pval AfpNum Identity(%) Similarity(%)\n";
	}
	
	while (<REP>) {
	    my $l = $_;

	    my ($protein_1, $len1, $protein_2, $len2,$twist_i, $iniLen, $iniRmsd, $optEqu,
		$optRmsd, $chainRmsd, $score, $aliLen_i, $gap_i, $pval_i, $afpNum, $ident, $sim) = split(/\s+/, $l);

	    if (length($gap)>0 && $gap_i>$gap) { next; }
	    if (length($twist)>0 && $twist_i>$twist) { next; }
	    if (length($alilen)>0 && $aliLen_i<$alilen) { next; }
	    if (length($equlen)>0 && $optEqu<$equlen) { next; }
	    if (length($rmsd)>0 && $optRmsd>$rmsd) { next; }
	    if (length($pval)>0 && $pval_i>$pval) { next; }
	    if (keys(%pairs_h) && !exists($pairs_h{"$protein_1:$protein_2"}) && !exists($pairs_h{"$protein_2:$protein_1"})) { next; }
	    print OUT "$l";
	}

	close(REP);
	close(OUT);
    }
}



sub readPair {
    my $pair_f = shift;
    my %pairs_h;

    if (open(PAIR, "<$pair_f")) {
	while (<PAIR>) {
	    s/^\s+//;
	    my ($p1, $p2) = split;
	    $p1 =~ s/\.pdb$//;
	    $p2 =~ s/\.pdb$//;
	    $pairs_h{"$p1:$p2"} = 1;
	    $pairs_h{"$p2:$p1"} = 1; 
	}
	close(PAIR);
    }
    return %pairs_h;
}


    
sub extactAlign {
    my ($ali_f, $pair_f,  $gap, $twist, $alilen, $equlen, $rmsd, $pval, $out_f) = @_;

    my %pairs_h = readPair($pair_f);
   
    my $aln = "";
    my @columns = ();
    open(ALN, "<$ali_f") || die "Can't open file $ali_f\n";
    open(OUT, ">$out_f") || die "Can't open file $out_f\n";
    while (<ALN>) {
	my $l = $_;
	if ($l =~ /^Align/) {
	    @columns = parseOnePair($aln);
	    my $pass = checkOnePair(\@columns, $gap, $twist, $alilen, $equlen, $rmsd, $pval, \%pairs_h);
	    if ($pass) { print OUT "$aln"; }
	    $aln = $l; 
	} else {
	    $aln .= $l;
	}
    }
    # print the last alignment
    @columns = parseOnePair($aln);
    my $pass = checkOnePair(\@columns, $gap, $twist, $alilen, $equlen, $rmsd, $pval, \%pairs_h);
    if ($pass) { print OUT "$aln"; }
    
    close(ALN);
    close(OUT);

}


sub checkOnePair {
    my ($columns_p, $gap, $twist, $alilen, $equlen, $rmsd, $pval, $pairs_hp) = @_;

    my ($protein_1, $len1, $protein_2, $len2, $factor, $ali) = @$columns_p;
    if (!$factor || !$ali) { return ""; } # no alignment
    
    my ($twist_i, $iniLen, $iniRmsd, $optEqu, $optRmsd, $chainRmsd, $score, $aliLen_i, $gap_i, $pval_i, $afpNum, $ident, $sim) = split(/\t+/, $factor);
    if (length($gap)>0 && $gap_i>$gap) { return ""; }
    if (length($twist)>0 && $twist_i>$twist) { return ""; }
    if (length($alilen)>0 && $aliLen_i<$alilen) { return ""; }
    if (length($equlen)>0 && $optEqu<$equlen) { return ""; }
    if (length($rmsd)>0 && $optRmsd>$rmsd) { return ""; }
    if (length($pval)>0 && $pval_i>$pval) { return ""; }
    if (keys(%$pairs_hp) && !exists($pairs_hp->{"$protein_1:$protein_2"}) && !exists($pairs_hp->{"$protein_2:$protein_1"})) { return ""; }
    return "1";
}
    
    

sub printReport {
    my ($ali_f, $gap, $twist, $alilen, $equlen, $rmsd, $pval, $out_f) = @_;
    
    
    open(ALN, "<$ali_f") || die "Can't open file $ali_f\n";
    open(OUT, ">$out_f") || die "Can't open file $out_f\n";
    
    if ($header =~ /YES/i) {
	print OUT "# pdb1 pdb1 len1 len2 twist iniLen iniRmsd optLen optRmsd chainRms score alnLen gap pval AfpNum Identity(%) Similarity(%)\n";
    }

    my $aln = "";
    my @columns = ();
    while (<ALN>) {
	my $l = $_;
	if ($l =~ /^Align/) {
	    @columns = parseOnePair($aln);
	    my $report_line = printOnePair(\@columns, $gap, $twist, $alilen, $equlen, $rmsd, $pval);
	    print OUT "$report_line";
	    $aln = $l; 
	} else {
	    $aln .= $l;
	}
    }
    # print the last alignment
    @columns = parseOnePair($aln);
    my $report_line = printOnePair(\@columns, $gap, $twist, $alilen, $equlen, $rmsd, $pval);
    print OUT "$report_line";
    
    close(ALN);
    close(OUT);
}


sub parseOnePair{
    my ($aln) = @_;
    
    my $ALN_START = 14; # the alignment lines start at position 14, ie 15th char

    my $factor = '';
    my ($protein_1, $protein_2) = ('',, '');
    my ($len1, $len2) = ('','');
    my ($chain_1, $blockIndex, $chain_2) = ('','','');
    
    my @lines = split(/\n/, $aln);
    
    my $i = 0;

    while ($i <= $#lines) {
	my $line = $lines[$i];
	
	if ($line =~ /^Align\s+(\S+)\s+(\d+)\s+with\s+(\S+)\s+(\d+)/) {
	    ($protein_1, $len1, $protein_2, $len2) = ($1,$2,$3,$4);
	    $protein_1 =~ s/\.pdb$//;
	    $protein_2 =~ s/\.pdb$//;
	} elsif ($line =~ /^Twists/) {
	    my @words = split(/\s+/, $line);
	    for (my $j = 1; $j < $#words; $j += 2) {
		$factor .= "$words[$j]\t";
	    }
	} elsif ($line =~ /^P\-value/) {
	    my @words = split(/\s+/, $line);
	    chop($words[5]);
	    chop($words[7]);
	    $factor .= "$words[1]\t$words[3]\t$words[5]\t$words[7]";
	} elsif ($line =~ /^Chain 1/) {
	    
	    # parse chain 1 line
	    $chain_1 .= substr($line, $ALN_START);   

	    # parse the twist blocks line
	    $blockIndex .= substr($lines[++$i], $ALN_START);

	    # parse the chain 2 line
	    $chain_2 .= substr($lines[++$i], $ALN_START);
	    
	}
	$i++;
    }
    if (!$chain_1 || !$chain_2 || !$blockIndex) { return (); }  # no alignment;
    return ($protein_1, $len1, $protein_2, $len2, $factor, "$chain_1:$blockIndex:$chain_2");
}
	    
	

sub  printOnePair {
    my ($columns_p, $gap, $twist, $alilen, $equlen, $rmsd, $pval) = @_;

    my ($protein_1, $len1, $protein_2, $len2, $factor, $ali) = @$columns_p;
    if (!$factor || !$ali) { return ""; } # no alignment
    
    my ($twist_i, $iniLen, $iniRmsd, $optEqu, $optRmsd, $chainRmsd, $score, $aliLen_i, $gap_i, $pval_i, $afpNum, $ident, $sim) = split(/\t+/, $factor);
    if (length($gap)>0 && $gap_i>$gap) { return ""; }
    if (length($twist)>0 && $twist_i>$twist) { return ""; }
    if (length($alilen)>0 && $aliLen_i<$alilen) { return ""; }
    if (length($equlen)>0 && $optEqu<$equlen) { return ""; }
    if (length($rmsd)>0 && $optRmsd>$rmsd) { return ""; }
    if (length($pval)>0 && $pval_i>$pval) { return ""; }
    return "$protein_1\t$len1\t$protein_2\t$len2\t$factor\n";
}
    
