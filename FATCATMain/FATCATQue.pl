#!/usr/bin/perl -w

#program for runing FATCAT on a list of structure pairs
#Y.Y, latest update 09/06/03
#note: the input file lists the code (not pdb file name), pdb file is in code.pdb format

#it requires FATCAT environment variate

my $fatcat;
if(!($fatcat = $ENV{'FATCAT'}))  { die "Stop: FATCAT environment variate is not found\n"; }

my $prog = "$fatcat/FATCATMain/FATCAT";

my $begtime = time();

die "FATCATQue.pl logfile listfile parameters..\n" unless (@ARGV >= 1);
my $logfile = shift;
my $listfile = shift;
my @para = @ARGV;

open(LOG, ">$logfile") || die "can't open file $logfile\n";

my $curr = `pwd`;
open($IN, $listfile) || die "open $listfile error in FATCATQue, curr-dir $curr\n";

my ($code1, $code2, $pdb1, $pdb2);

while(<$IN>)	{
    if(/^\#/)	{ next; }
    ($code1, $code2) = split;
       
    $command = $prog;
    if($code1 =~ s/^([^:]+):(\d+)\+(\d+)/$1/)	{
	$command = "$command -s1 $2 -l1 $3";
    } # define beginning and ending positions for protein 1

    $pdb1 = "$code1.pdb";
    $command = "$command -p1 $pdb1";
    
    if($code2 =~ s/^([^:]+):(\d+)\+(\d+)/$1/)	{
	$command = "$command -s2 $2 -l2 $3";
    } # define beginning and ending positions for protein 2
    
    $pdb2 = "$code2.pdb";
    $command = "$command -p2 $pdb2";
    
    $command = "$command -o $code1.$code2 @para";
    
    #print "command: $command\n";
    system($command);
}
close($IN);

my $endtime = time();
my $timeuse = ($endtime - $begtime);
my $unit = "mins";
if($timeuse > 3600)	{ 
    $timeuse /= 3600; $unit = "hours"; 
} else	{
    $timeuse /= 60; $unit = "mins";
}
my $timeused = sprintf("#Time used %.2f $unit\n", $timeuse);
print LOG $timeused;
close(LOG);
