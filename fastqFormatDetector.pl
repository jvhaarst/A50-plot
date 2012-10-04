#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(max min);

# initiate
my $sequence;
my $skipline;
my $quality;
my @characters;
my @numbers;
my $minimum=255;
my $maximum=0;
my $DEBUG=0;

# Go through the file
while(<STDIN>){					# Skip first line with identifier
    $sequence = <STDIN>;			# Skip second line with sequence
    $skipline = <STDIN>;			# Skip third line with '+'
    $quality =  <STDIN>; chomp $quality;	# Fourth line contains quality
    @characters = split(//,$quality);		# Divide into characters
    # Convert the chracters on the array into their equivalent ASCII codes
    # And thus into the quality values
    @numbers = map(ord, @characters);
    if ($minimum > min(@numbers)){		# A new minimum was found
    	$minimum = min(@numbers) ;
    	print $quality."\n" if $DEBUG;
    	print "minimum changed to $minimum(".chr($minimum).")\n" if $DEBUG;
    }
    if ($maximum < max(@numbers)){		# A new maximum was found
    	$maximum = max(@numbers);
    	print $quality."\n" if $DEBUG;
    	print "maximum changed to $maximum(".chr($maximum).")\n" if $DEBUG;
    }
}
print "minimum is: $minimum(".chr($minimum).")\n";
print "maximum is: $maximum(".chr($maximum).")\n";

if ($minimum == 33 && $maximum ==  73){ print "Sanger\n"; }
if ($minimum == 59 && $maximum == 104){ print "Solexa\n"; }
if ($minimum == 64 && $maximum == 104){ print "Illumina 1.3\n"; }
if ($minimum == 66 && $maximum == 104){ print "Illumina 1.5\n"; }
if ($minimum == 33 && $maximum ==  74){ print "Illumina 1.8\n"; }

