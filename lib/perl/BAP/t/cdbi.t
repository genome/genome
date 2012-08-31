#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 5042;

use BAP::DB::Organism;
use BAP::DB::SequenceSet;


$BAP::DB::DBI::db_env = 'dev';

my $organism = BAP::DB::Organism->insert({
                                          organism_name    => 'automated_test',
                                          ncbi_taxonomy_id => 31337,
                                          gram_stain       => '+',
                                          locus            => 'ATMTD_TST',
                                         });

my $sequence_set = BAP::DB::SequenceSet->insert({
                                                 sequence_set_name => 'automated_test',
                                                 organism_id       => $organism->organism_id(),
                                                 software_version  => 0.0,
                                                 data_version      => '01-01-1970'
                                                });

foreach my $i (1..10) {

    my $sequence = BAP::DB::Sequence->insert({
                                              sequence_name   => 'automated_test',
                                              sequence_set_id => $sequence_set->sequence_set_id(),
                                              sequence_string => 'GATTACAAA' x 100,
                                             });

    for my $j (1..100) {

        my $start = 1 + (9 * ($j - 1));
        my $end = $start + 8;
        
        my $coding_gene = BAP::DB::CodingGene->insert({
                                                       gene_name       => join('.', 'automated_test', $i, $j),
                                                       sequence_id     => $sequence->sequence_id(),
                                                       sequence_string => 'GATTACAAA',
                                                       start           => $start,
                                                       end             => $end,
                                                       strand          => 1,
                                                       source          => 'automated_test',
                                                      });

        my $protein = BAP::DB::Protein->insert({
                                                protein_name    => $coding_gene->gene_name(),
                                                gene_id         => $coding_gene->gene_id(),
                                                sequence_string => 'DYK',
                                               });

    }

}

BAP::DB::DBI->dbi_commit();

ok(defined($organism));
ok(defined($organism->organism_id()));
ok($organism->organism_id() =~ /^\d+$/);
ok($organism->organism_name() eq 'automated_test');
ok($organism->ncbi_taxonomy_id() == 31337);
ok($organism->gram_stain() eq '+');
ok($organism->locus() eq 'ATMTD_TST');

ok(defined($sequence_set));
ok(defined($sequence_set->sequence_set_id()));
ok($sequence_set->sequence_set_id() =~ /^\d+$/);
ok($sequence_set->sequence_set_name() eq 'automated_test');

my @sequences = $sequence_set->sequences();

ok(scalar(@sequences) == 10);

foreach my $sequence (@sequences) {

    ok(length($sequence->sequence_string) == 900);
    ok($sequence->sequence_string eq ('GATTACAAA' x 100));

    my @coding_genes = $sequence->coding_genes();

    ok(scalar(@coding_genes) == 100);

    foreach my $coding_gene (@coding_genes) {

        ok($coding_gene->gene_id =~ /^\d+$/);
        ok($coding_gene->sequence_string eq 'GATTACAAA');
        
        my @proteins = $coding_gene->protein();

        ok(scalar(@proteins) == 1);
        ok($coding_gene->gene_name() eq $proteins[0]->protein_name());
        ok($proteins[0]->sequence_string() eq 'DYK');

    }

}

$sequence_set->delete();
$organism->delete();

BAP::DB::DBI->dbi_commit();
