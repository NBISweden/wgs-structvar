#!/usr/bin/perl
# SVvcf2bed.pl     pall@pallVB     2015/06/15 10:59:10

#use Modern::Perl 2011;
#use warnings;
#use strict; # No need with modern perl...
$|=1;
use Data::Dumper;
use POSIX qw(floor ceil);

my $DEBUG = 1;

chomp(my $inputfile = shift);
open(IN, "< $inputfile") or die; 

###
# This should be tuned for Fermikit and Manta as we want to use those
# two at present
#




###
# bed is zero based, start should be start -1, end = end
#


while (<IN>){

  next if (/^\#/);
  
  my @cols = split(/\t/);
  my $svtype = "";
  my $samples = "";

  if ($cols[7] =~ m/SVTYPE\=(\w+)?/){
    $svtype = $1;
    if ($cols[7] =~ m/SF\=([\d\,]+)/){$samples = $1}

      if ($svtype eq 'DEL'){
	if ($cols[7] =~ m/END\=(\d+)?\;/){
	  print $cols[0], "\t", $cols[1]-1, "\t", $1, "\t", $svtype, "\t", $samples, "\n";
	}

      } elsif ($svtype eq 'INS'){

	###
	# manta does not always have SVLEN!
	# 
	#
	# fermikit has end-start = |svlen-qgap|, order END, SVLEN, QGAP
	#

	# first case, fermikit, print start/end directly:
	if ($cols[7] =~ m/END\=(\d+)?\;.*SVLEN\=(\d+)\;.*QGAP\=([\-\d]+)\;/){
	  print $cols[0], "\t", $cols[1]-1, "\t", $1, "\t", $svtype, "\t", $samples, "\n";
	}
	
	# now manta, first case, when SVLEN is given:
	elsif ($cols[7] =~ m/END\=(\d+)?\;.*SVLEN\=(\d+)/){
	  print $cols[0], "\t", $cols[1]-1, "\t", $cols[1] + $2, "\t", $svtype, "\t", $samples, "\n";
 	}
	# no SVLEN, check the lenght of LEFT/RIGHT_SVINSSEQ:
	elsif ($cols[7] =~ m/END\=(\d+)?\;.*LEFT_SVINSSEQ\=([ACGT]+)\;RIGHT_SVINSSEQ\=([ACGT]+)/){
	  my $svlen = length($2) + length($3);
	  print $cols[0], "\t", $cols[1]-1, "\t", $cols[1] + $svlen, "\t", $svtype, "\t", $samples, "\n";

	}
	
	else {
	  warn "Non standard INS found...\n" if ($DEBUG > 0);
	}

      } elsif ($svtype eq 'COMPLEX'){
	my ($start,$stop) = coord2kbint($cols[1]);
	print $cols[0], "\t", $start, "\t", $stop, "\t", $svtype, "\t", $samples, "\n";
	
      } else {
	
	if ($cols[7] =~ m/END\=(\d+)?\;/){
	  print $cols[0], "\t", $cols[1]-1, "\t", $1, "\t", $svtype, "\t", $samples, "\n";
	}
	#die "unknown SV type $svtype...exit\n";

      }
  }

}

close IN;

sub coord2kbint {
  my $coord = shift;
  return(floor($coord/1000)*1000, floor($coord/1000)*1000+999);
}


=head1 NAME

SVvcf2bed.pl

=head1 SYNOPSIS

SVvcf2bed.pl <file.vcf>

=head1 DESCRIPTION

Convert fermi/MANTA structural variant vcf calls to bed.
Start/End coords for COMPLEX events are not trivial to deduce, so I simply report
the kb window they occur in. 


=head1 AUTHOR

Pall Isolfur Olason

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Pall Isolfur Olason

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
