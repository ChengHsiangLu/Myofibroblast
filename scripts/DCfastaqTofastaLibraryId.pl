#!/usr/bin/perl -w
use strict;

############external parameters:

use vars qw($opt_L);
use Getopt::Std;
getopts("L:");

my @temptag=();
my $ln;
my $i=0;
my $flag = 0;

while($ln = <>) {
	chomp($ln);
	if ($ln =~ /^@/o) { 
		$flag = 1;
	}
	else {
		if ($flag == 1 && $ln =~ /^[ACGTN]/o) {
			$i+=1;
			printf(STDOUT "%s\n%s\n",">".$opt_L."_".$i,$ln);
			$flag = 0;
		}
	}
}

