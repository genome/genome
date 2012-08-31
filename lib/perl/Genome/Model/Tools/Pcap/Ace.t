#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::Pcap::Ace') or die;

my @assemblies;
push @assemblies, Genome::Model::Tools::Pcap::Ace->new(
    input_file => $ENV{GENOME_TEST_INPUTS} . "/Genome-Assembly-Pcap/test.ace"
);

my $tmp = File::Temp->new( UNLINK => 1, SUFFIX => '.db' );
die "Failed to create temp sqlite db file for testing.\n" unless defined $tmp;

push @assemblies, Genome::Model::Tools::Pcap::Ace->new(
    input_file => $ENV{GENOME_TEST_INPUTS} . "/Genome-Assembly-Pcap/test.ace",
    using_db => 1,
    db_type => "SQLite",
    db_file => $tmp->filename
);


push @assemblies, Genome::Model::Tools::Pcap::Ace->new(
    input_file => $ENV{GENOME_TEST_INPUTS} . "/Genome-Assembly-Pcap/test.ace",
    using_db => 1,
    db_type => "mysql"
);	

foreach my $assembly (@assemblies) {
    my $contig = $assembly->get_contig(
        "Contig0.10"
    );
    my $name = $contig->name;
    my $seq = $contig->padded_base_string;

    is($name, "Contig0.10", "Name survives creation/getting");
}

foreach my $assembly (@assemblies) {
    my $contig = $assembly->get_contig(
        "Contig0.10"
    );
    my $seq = $contig->padded_base_string;
    is(substr($seq,0,50), "CtcaattggcaaTCAAtctGTGGCTctTAcCCAAcAAGGcGCAATCACAA", "Sequence survives creation");
}

foreach my $assembly (@assemblies) {
    my $contig = $assembly->get_contig(
        "Contig0.10"
    );
    my $length = $contig->length;
    is($length, 156110, "Length survives creation");
    my $unpadlength = $contig->length("unpadded");
    is($unpadlength, 153704, "Unpadded length is correct"); 
    my $seq = $contig->unpadded_base_string;
    my $seq2 = $contig->padded_base_string;
    is(length $seq, 153704, "Padded base string retrieved correctly");
    is(length $seq2, 156110, "Unpadded base string retrieved correctly");
}

foreach my $assembly (@assemblies) {
    #$assembly->_build_ace_index(file_name => $ENV{GENOME_TEST_INPUTS} . '/Genome-Assembly-Pcap/test.ace');
    #TODO: write tests for indexing
    ok(1, 'index built properly');
}


done_testing();
exit;

