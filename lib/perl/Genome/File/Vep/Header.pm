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
    my ($class, $txt) = @_;
    if (defined $txt) {
        chomp $txt;
        my $expected = join("\t", @expected_header);
        confess "Invalid Vep header:\n$txt\nExpected:\n$expected" unless $expected eq $txt;
    }
    return bless { fields => \@expected_header }, $class;
}

sub to_string {
    return join("\t", @expected_header);
}

1;
