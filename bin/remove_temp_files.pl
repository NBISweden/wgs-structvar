#!/usr/bin/env perl
# Resolving the following issue will make this script obsolete:
# https://github.com/nextflow-io/nextflow/issues/165

use strict;
use warnings;
use feature qw( say );

use File::Spec;
use File::Path qw( remove_tree );
use File::Glob qw( bsd_glob );


my $log_file = shift // '.nextflow.log';

## This regex is very brittle, need to debug carefully, what happens with
## spaces for example?
## Example line (seen null instead of lustre once)              (......................)
# Aug-24 12:25:33.343 [main] DEBUG nextflow.Session - Work-dir: /home/rajohvik/glob/work [lustre]
my $work_dir_rx = qr{
    Work-dir: \s*
    (.*?)
    \s \[.*?\]
}x;

## I want the sha checksum.
## Example line                                                 (.........)
# Aug-24 12:25:37.224 [Actor Thread 7] INFO  nextflow.Session - [be/91a2c1] Submitted process > download_masks (1)
my $session_id_rx = qr{
    INFO \s+ nextflow\.Session \s+ - \s+
    \[ ([^]]+) \]
    \s+ Submitted
}x;

my $work_dir;
my @sessions;

open my $LOG, '<', $log_file or die "Can't open logfile <$log_file>: $!\n";
while (my $line = <$LOG>) {
    if ( $line =~ $work_dir_rx ) {
        $work_dir = $1;
        next;
    }
    if ( $line =~ $session_id_rx ) {
        push @sessions, $1;
    }
}

for my $session ( @sessions ) {
    my $dir = File::Spec->catfile($work_dir, $session);
    my @dir_candidates = bsd_glob("${dir}*"); # Treat space in pattern properly

    if ( @dir_candidates > 1 ) {
        say STDERR "Found several candidate dirs for $session: @dir_candidates";
        next;
    }
    if ( @dir_candidates == 0 ) {
        say STDERR "Can't find directory for $session ($dir)";
        next;
    }
    if ( ! -d $dir_candidates[0] ) {
        say STDERR "Can't find directory for $session ($dir)";
        next;
    }

    remove_tree( $dir_candidates[0], { safe => 1 });
    say "Cleaned $session (@dir_candidates)";
}
