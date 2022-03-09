#!/usr/bin/perl -w

#runing FATCAT for a query against a list of structure pairs
#note: when it is used for database-searching on cluster, be sure to use option "-q" in the parameters
#Y.Y, latest update 09/06/03
#Y.Y, update, feb 25,04, add the database-split function 

#it requires FATCAT environment variate

my $fatcat;
if(!($fatcat = $ENV{'FATCAT'}))  { die "Stop: FATCAT environment variate is not found\n"; }

my $prog = "$fatcat/FATCATMain/FATCAT";

my $usage = "FATCATSearch.pl query target-list FATCAT-parameters(refer FATCAT usage)\n".
	    "  when only part of the items from the target-list are used (e.g divide job in cluster), use target-list:beg:end\n".
	    "  for instance, /home/yye/data/pdb/scop163_95.descript:0:99, i.e search query against items 0-99 (100 in total)\n";
die "FATCATSearch.pl query target-list parameters..\n" unless (@ARGV >= 1);
my $code1 = shift;
my $listfile = shift;
my @para = @ARGV;
my $begitem = 0;
my $enditem = 1000000; # a very large number
for(my $i = 0; $i < @para; $i ++)	{
	if($para[$i] eq "-beg" && @para > $i + 1)	{ $begitem = $para[$i + 1]; }
	elsif($para[$i] eq "-end" && @para > $i + 1)	{ $enditem = $para[$i + 1]; }
}

my $command0 = "$prog -p1 $code1";
$code1 =~ s/\.pdb$//;
my $command;

my $item = 0;
my $IN;
open($IN, $listfile) || die "open $listfile error\n";
while(<$IN>)	{
	if(/^#/)	{ next; }
	if($item < $begitem || $item > $enditem)	{ $item ++; next; }
	$command = $command0;
	($code2) = split;
	if($code2 =~ s/^([^:]+):(\d+)\+(\d+)/$1/)	{
		$command = "$command -p2 $code2.pdb -s2 $2 -l2 $3";
	} # specificy starting residue and length
	else	{ $command = "$command -p2 $code2.pdb"; }
	$command = "$command -o $code1.$code2 @para";
	#print "command: $command\n";
	system($command);
	$item ++;
}
close($IN);
