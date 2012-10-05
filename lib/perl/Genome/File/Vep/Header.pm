package Genome::File::Vep::Header;

use Carp qw/confess/;
use Genome;

use strict;
use warnings;

my @expected_header = qw/
    Uploaded_variation
    Location
    Allele
    Gene
    Feature
    Feature_type
    Consequence
    cDNA_position
    CDS_position
    Protein_position
    Amino_acids
    Codons
    Existing_variation
    Extra
    /;

sub new {
    my ($class, $lines) = @_;
    
    my $expected = "#".join("\t", @expected_header);
    if (!$lines) {
        return $class->new([$expected]);
    }

    my @meta;
    my @fields;

    my $idx = 0;
    chomp @$lines;
    for (; $idx <= $#$lines; ++$idx) {
        my $line = $lines->[$idx];

        if ($line =~ /^##/) {
            push(@meta, $line);
        } elsif ($line =~ /^#/) {
            confess "Invalid Vep header:\n$line\nExpected:\n$expected" unless $expected eq $line;
            $line =~ s/^#//g;
            @fields = split("\t", $line);
            last;
        } else {
            confess "Unexpected data in vep header: $line";
        }
    }

    if (@fields == 0 || $idx != $#$lines) {
        confess "Invalid Vep header:\n\t" . join("\n\t", @$lines);
    }

    return bless { _meta => \@meta, _fields => \@fields }, $class;
}

sub to_string {
    my $self = shift;
    return join("\n", @{$self->{_meta}}, "#".join("\t", @{$self->{_fields}}));
}

1;
