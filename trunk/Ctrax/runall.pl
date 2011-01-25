#!/usr/bin/perl

use strict;

my $outdir = "/localhome/bransonk/tracking_data/olympiad/FlyBowl/CtraxTest20101118";

my $filename = $ARGV[0];

open(FILE,"<$filename");
my @childs = ();

my $count = 1;
while(my $line = <FILE>){
    chomp $line;
    if(length($line) == 0){
	next;
    }
    if($line =~ /^#/){
	next;
    }
    $line =~ /(^.+\/)?([^\/]+$)/;
    my $rootdir = $1;
    my $expdir = $2;
    my $cmd = "./ctraxcmd20110111.sh $expdir";
    my $pid = fork();
    if($pid){
	push(@childs, $pid);
	$count++;
    }
    else{
	sleep($count);
	print "$cmd\n";
	`$cmd`;
	exit(0);
    }
	
}

close(FILE);

foreach (@childs) {
    waitpid($_, 0);
    print "child $_ finished\n";
}
