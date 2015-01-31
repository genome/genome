package Genome::Model::Tools::Vcf::VcfCompare;

#this package is for parsing of joinx's vcf-compare output

use strict;
use warnings;

use IO::File;

#File looks like the below:
# 0	1	partial_match	369
# 1	0	partial_match	316
# 1	1	partial_match	16816
# 0	1	exact_match	375
# 1	0	exact_match	317
# 1	1	exact_match	16810
# 0	1	partial_miss	0
# 1	0	partial_miss	5
# 1	1	partial_miss	0
# 0	1	complete_miss	369
# 1	0	complete_miss	311
# 1	1	complete_miss	0


sub new {
    $DB::single = 1;
    my ($class, $file) = @_;

    my $self;
    my $fh = IO::File->new($file) or die "Unable to open $file to parse vcf-compare results\n";

    my $header_line = $fh->getline;
    chomp $header_line;
    my ($file1_name, $file2_name, $type, @samples) = split "\t", $header_line;
    $self->{_file_names} = [$file1_name, $file2_name]; 

    while(my $line = $fh->getline) {
        chomp $line;
        my ($file1, $file2, $type, @counts) = split "\t", $line;
        @{$self->{_melty_table}{$file1}{$file2}{$type}}{@samples} = @counts;
    }
    bless $self, $class;
}

sub unique_count {
    my ($self, $file, $type, $sample) = @_;
    my $count_pointer = $self->{_melty_table};
    for my $name (@{$self->{_file_names}}) {
        if($name eq $file) {
            $count_pointer = $count_pointer->{1};
        }
        else {
            $count_pointer = $count_pointer->{0};
        }
        unless(defined $count_pointer) {
            die "$file not present in joinx output\n";
        }
    }
    if(exists($count_pointer->{$type})) {
        if(exists($count_pointer->{$type}{$sample})) {
            return $count_pointer->{$type}{$sample};
        }
        else {
            die "$sample not present in joinx output\n";
        }
    }
    else {
        die "$type not valid in joinx output\n";
    }
}

sub joint_count {
    my ($self, $type, $sample) = @_;
    if(exists($self->{_melty_table}{1}{1}{$type})) {
        if(exists($self->{_melty_table}{1}{1}{$type}{$sample})) {
            return $self->{_melty_table}{1}{1}{$type}{$sample};
        }
        else {
            die "$sample not present in joinx output\n";
        }
    }
    else {
        die "$type not valid in joinx output\n";
    }
}

1;
