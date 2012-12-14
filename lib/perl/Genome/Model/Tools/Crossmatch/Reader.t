#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper;
use Test::More;

use_ok('Genome::Model::Tools::Crossmatch::Reader') or die;

my $file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Crossmatch/reader.crossmatch';
my $reader = Genome::Model::Tools::Crossmatch::Reader->create(input => $file);
ok($reader, 'create reader');
my $aln_cnt = 0;
my @alns = _alns();
while ( my $aln = $reader->next ) {
    is_deeply($aln, $alns[$aln_cnt], 'alignment matches');
    $aln_cnt++;
}
is($aln_cnt, @alns, 'got correct number of alignments');

done_testing();


###

sub _alns {
    return (
          {
            'per_sub' => '0.57',
            'discrepancies' => [
                                 {
                                   'base' => 'A',
                                   'mutation' => 'S',
                                   'sequence' => 'ttctctAtttgct',
                                   'number' => 1,
                                   'query_pos' => '111',
                                   'subject_pos' => '681'
                                 },
                                 {
                                   'base' => 'C',
                                   'mutation' => 'S',
                                   'sequence' => 'ctccacCaggcgc',
                                   'number' => 1,
                                   'query_pos' => '172',
                                   'subject_pos' => '742'
                                 },
                                 {
                                   'base' => 'C',
                                   'mutation' => 'S',
                                   'sequence' => 'tcatatCgccgct',
                                   'number' => 1,
                                   'query_pos' => '575',
                                   'subject_pos' => '1145'
                                 },
                                 {
                                   'base' => 'A',
                                   'mutation' => 'S',
                                   'sequence' => 'aaggaaAaggcag',
                                   'number' => 1,
                                   'query_pos' => '660',
                                   'subject_pos' => '1230'
                                 },
                                 {
                                   'base' => 'T',
                                   'mutation' => 'S',
                                   'sequence' => 'atcttcTtcagag',
                                   'number' => 1,
                                   'query_pos' => '810',
                                   'subject_pos' => '1380'
                                 }
                               ],
            'sw_score' => '821',
            'per_ins' => '0.00',
            'query_name' => 'Contig-1.3',
            'bases_before' => '280',
            'subject_stop' => '1453',
            'query_start' => '1',
            'subject_start' => '571',
            'per_del' => '0.00',
            'bases_after' => '0',
            'subject_name' => '/gscuser/kchen/sata114/kchen/Hs_build36/all_fragments/Homo_sapiens.NCBI36.45.dna.chromosome.16.fa',
            'query_stop' => '883'
          },
          {
            'per_sub' => '0.00',
            'sw_score' => '38',
            'per_ins' => '0.00',
            'query_name' => 'Contig1.4.1.-2',
            'subject_stop' => '572',
            'bases_before' => '1122',
            'query_start' => '1',
            'subject_start' => '611',
            'per_del' => '0.00',
            'bases_after' => '453',
            'subject_name' => '/gscuser/kchen/sata114/kchen/Hs_build36/all_fragments/Homo_sapiens.NCBI36.45.dna.chromosome.16.fa',
            'query_stop' => '40'
          },
          {
            'per_sub' => '1.77',
            'discrepancies' => [
                                 {
                                   'base' => 'A',
                                   'mutation' => 'S',
                                   'sequence' => 'cggacgAgcggct',
                                   'number' => 1,
                                   'query_pos' => '103',
                                   'subject_pos' => '549'
                                 },
                                 {
                                   'base' => 'G',
                                   'mutation' => 'S',
                                   'sequence' => 'tgttctGttttac',
                                   'number' => 1,
                                   'query_pos' => '141',
                                   'subject_pos' => '511'
                                 },
                                 {
                                   'base' => 'G',
                                   'mutation' => 'S',
                                   'sequence' => 'ggaaggGgcatca',
                                   'number' => 1,
                                   'query_pos' => '194',
                                   'subject_pos' => '458'
                                 },
                                 {
                                   'base' => 'A',
                                   'mutation' => 'S',
                                   'sequence' => 'gttgccAtatctc',
                                   'number' => 1,
                                   'query_pos' => '247',
                                   'subject_pos' => '405'
                                 },
                                 {
                                   'base' => 'G',
                                   'mutation' => 'S',
                                   'sequence' => 'ctttatGaacaat',
                                   'number' => 1,
                                   'query_pos' => '277',
                                   'subject_pos' => '375'
                                 },
                                 {
                                   'base' => 'C',
                                   'mutation' => 'S',
                                   'sequence' => 'ccaccaCgcgtgg',
                                   'number' => 1,
                                   'query_pos' => '431',
                                   'subject_pos' => '221'
                                 },
                                 {
                                   'base' => 'G',
                                   'mutation' => 'S',
                                   'sequence' => 'ccacgcGtggcta',
                                   'number' => 1,
                                   'query_pos' => '434',
                                   'subject_pos' => '218'
                                 },
                                 {
                                   'base' => 'A',
                                   'mutation' => 'D',
                                   'sequence' => 'tggctaAtttttt',
                                   'number' => 1,
                                   'query_pos' => '441',
                                   'subject_pos' => '210'
                                 },
                                 {
                                   'base' => 'T',
                                   'mutation' => 'S',
                                   'sequence' => 'tctgtgTtgccca',
                                   'number' => 1,
                                   'query_pos' => '477',
                                   'subject_pos' => '174'
                                 }
                               ],
            'sw_score' => '354',
            'per_ins' => '0.00',
            'query_name' => 'Contig1.4.1.-2',
            'subject_stop' => '158',
            'bases_before' => '1122',
            'query_start' => '41',
            'subject_start' => '611',
            'per_del' => '0.22',
            'bases_after' => '0',
            'subject_name' => '/gscuser/kchen/sata114/kchen/Hs_build36/all_fragments/Homo_sapiens.NCBI36.45.dna.chromosome.16.fa',
            'query_stop' => '493'
          }
    );
}
