#!/usr/bin/env genome-perl

use strict;
use warnings;
use Data::Dumper;
use IO::Dir;
use IO::File;

my $scaf_merge_directory = '/gscmnt/temp113/finishing/scratch/chimp_merge/analysis/scaf_join/join_out';

my $dir = IO::Dir->new($scaf_merge_directory);
my $fail;
my %merges;

die unless $dir;
my $counter;
my $total;
while (my $file = $dir->read){ #and ++$counter < 5){
    next if $file =~ /^\.+$/;
    $total++;
    scoop_from_file("$scaf_merge_directory/$file");
}
my @scaffold_nums;
my %danger_scafs;
foreach my $key(keys %merges){
    my @merge = sort {$a->{position} <=> $b->{position}} @{$merges{$key}};
    push @scaffold_nums,$_->{scaffold_num} foreach @merge;
}
while (my $num = shift @scaffold_nums){
    my $times = grep {$num eq $_} @scaffold_nums;
    print "$num is merged ".++$times." times!\n" and $danger_scafs{$num} = $times if $times and !$danger_scafs{$num};
}
foreach my $key(keys %merges){
    my @merge = sort {$a->{position} <=> $b->{position}} @{$merges{$key}};
    my $print;
    my $repeat_num;
    foreach my $num (keys %danger_scafs){
        if (grep { $num eq $_->{scaffold_num} } @merge){
            $print = 1;
            $repeat_num = $num;
        }
    }
    print "$key : $repeat_num\n" if $print;
}



sub scoop_from_file{
    my $file = shift;
    my $io = IO::File->new("< $file");
    my $storable_string;
    my $counter;
    my @datums;
    while (my $line = $io->getline){
        next unless $line =~ /VAR1.*bless/;
        @datums = parse_scaffold($io);
    }
    $merges{$file} = \@datums;
}

sub parse_scaffold{
    my $io = shift;
    my @datums;
    while (my $line = $io->getline){
        last if $line =~ /\);/;
        if ($line =~/bless/){
            $line = $io->getline;
            my ($orientation) = $line =~ /orientation.*=>.*'([\+-])'/;
            $line = $io->getline;
            my ($position) = $line =~ /position.*=>.*(\d+)/;
            $line = $io->getline;
            my ($scaffold_num) = $line =~ /Contig(\d+)/;
            my $datum = {orientation => $orientation, position => $position, scaffold_num => $scaffold_num};
            push @datums, $datum;

        }
    }
    return @datums;
}

=pod

=head1 NAME
ScriptTemplate - template for new perl script

=head1 SYNOPSIS

=head1 DESCRIPTION 

=cut

#$HeadURL$
#$Id$


