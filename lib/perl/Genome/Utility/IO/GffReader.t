#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 3;

use_ok('Genome::Utility::IO::GffReader');
my $gff_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-IO-GffReader/alignment_summary.gff';

my $reader = Genome::Utility::IO::GffReader->create(
   input => $gff_file,
);
isa_ok($reader,'Genome::Utility::IO::GffReader');
my $data = $reader->next;
my @headers = keys %{$data};
my $expected_headers = Genome::Utility::IO::GffReader->headers;
ok( (scalar(@headers) == scalar(@{$expected_headers}) ), 'Found expected header count.' );

$data = $reader->next_with_attributes_hash_ref;
my $attributes = $data->{attributes_hash_ref};
#print Data::Dumper::Dumper($attributes) ."\n";

my $gtf_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-IO-GffReader/annotation.gtf';
my $gtf_reader = Genome::Utility::IO::GffReader->create(
   input => $gtf_file,
);
while (my $gtf_data = $gtf_reader->next_with_attributes_hash_ref) {
    my $attributes = $gtf_data->{attributes_hash_ref};
    #print Data::Dumper::Dumper($attributes) ."\n";
}
