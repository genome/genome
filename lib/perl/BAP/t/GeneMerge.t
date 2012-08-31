use strict;
use warnings;

use Bio::SeqFeature::Generic;

use Data::Dumper;
use Test::More tests => 28;

BEGIN {
    use_ok('BAP::GeneMerge');
}

my @features = (
                create_feature( 2112, 2191,  1, 'TEST01' ),
                create_feature( 2191, 2271,  1, 'TEST02' ),
                create_feature( 2030, 2111,  1, 'TEST03' ),
                create_feature( 2100, 2200, -1, 'TEST04' ),
               );

my $graph = BAP::GeneMerge::graph_overlaps(1, 1, \@features);

ok($graph->has_edge('TEST01', 'TEST02'), 'TEST01 overlaps TEST02');
ok($graph->has_edge('TEST01', 'TEST04'), 'TEST01 overlaps TEST04');
ok($graph->has_edge('TEST02', 'TEST04'), 'TEST02 overlaps TEST04');
ok($graph->has_edge('TEST03', 'TEST04'), 'TEST03 overlaps TEST04');

ok(!$graph->has_edge('TEST01', 'TEST03'), 'TEST01 does not overlap TEST03');
ok(!$graph->has_edge('TEST02', 'TEST03'), 'TEST02 does not overlap TEST03');

my $pfam_feature   = Bio::SeqFeature::Generic->new(-tag => { 'pfam_evidence'   => 1 });

my $short_feature = Bio::SeqFeature::Generic->new(
                                                  -start => 1,
                                                  -end   => 100,
                                                 );

my $long_feature = Bio::SeqFeature::Generic->new(
                                                 -start => 1,
                                                 -end   => 200,
                                                );
                                                
my $short_blastp_feature = Bio::SeqFeature::Generic->new(
                                                         -start => 1,
                                                         -end   => 100,
                                                         -tag => { 'blastp_evidence' => 1 },
                                                        );

my $long_blastp_feature = Bio::SeqFeature::Generic->new(
                                                        -start => 1,
                                                        -end   => 200,
                                                        -tag   => { 'blastp_evidence' => 1 },
                                                       );
                                                       
ok(!BAP::GeneMerge::better_gene($pfam_feature, $pfam_feature), 'pfam v pfam');
ok(!BAP::GeneMerge::better_gene($pfam_feature, $short_blastp_feature), 'pfam v blastp');
ok(BAP::GeneMerge::better_gene($short_blastp_feature, $pfam_feature), 'blastp v pfam');
ok(BAP::GeneMerge::better_gene($short_blastp_feature, $long_blastp_feature), 'blastp v blastp');
ok(BAP::GeneMerge::better_gene($short_feature, $long_feature), 'no evidence v no evidence');
ok(!BAP::GeneMerge::better_gene($short_blastp_feature, $short_feature), 'blastp v no evidence');
ok(BAP::GeneMerge::better_gene($short_feature, $short_blastp_feature), 'no evidence v blastp');
        
@features = (
             create_feature( 1,     901,  1, 'TEST05' ),
             create_feature( 600,  1500,  1, 'TEST06' ),
             create_feature( 200,  2500,  1, 'TEST07' ),
             create_feature( 2200, 3100,  1, 'TEST08' ),
            );

$features[2]->add_tag_value('blastp_evidence');
$features[3]->add_tag_value('pfam_evidence');

$graph = BAP::GeneMerge::graph_overlaps(300, 40, \@features);
BAP::GeneMerge::fancy_tag_overlapping($graph, \@features);

ok(!$features[0]->has_tag('delete_overlap'), 'TEST05 is not deleted');
ok(!$features[1]->has_tag('delete_overlap'), 'TEST06 is not deleted');
ok($features[2]->has_tag('delete_overlap'), 'TEST07 is deleted');
ok(!$features[3]->has_tag('delete_overlap'), 'TEST08 is not deleted');
ok($features[2]->has_tag('delete_for'), 'TEST07 has delete_for tag');

my ($delete_for) = $features[2]->each_tag_value('delete_for');

is($delete_for, 'TEST08', 'TEST07 deleted in favor of TEST08');


@features = (
             create_feature( 1,   100, -1, 'rfam_1',    'rfam'    ),
             create_feature( 100, 200,  1, 'rnammer_1', 'rnammer' ),
             create_feature( 300, 400,  1, 'rfam_2',    'rfam'    ),
             create_feature( 401, 501,  1, 'rnammer_2', 'rnammer' ),
            );

BAP::GeneMerge::tag_redundant_rfam(\@features);

ok($features[0]->has_tag('redundant'),  'rfam_1 is redundant');
ok(!$features[1]->has_tag('redundant'), 'rnammer_1 is not redundant');
ok(!$features[2]->has_tag('redundant'), 'rfam_2 is not redundant');
ok(!$features[3]->has_tag('redundant'), 'rnammer_2 is not redundant');

@features = (
             create_feature( 1,    1000, 1, 'glimmer3_1', 'glimmer3' ),
             create_feature( 1,    100,  1, 'rfam_1',     'rfam'     ),
            );

$features[0]->add_tag_value('type'       => 'coding');
$features[1]->add_tag_value('type'       => 'rRNA'  );

BAP::GeneMerge::tag_rna_overlap(\@features);

ok($features[0]->has_tag('delete_rrna_overlap'), 'glimmer3_1 is deleted');

@features = (
             create_feature( 1,    1000, 1, 'glimmer3_1', 'glimmer3' ),
             create_feature( 991,  1090, 1, 'rfam_2',     'rfam'     ),
             create_feature( 951,  1050, 1, 'rfam_3',     'rfam'     ),
            );

$features[0]->add_tag_value('type'       => 'coding');
$features[1]->add_tag_value('type'       => 'rRNA'  );
$features[2]->add_tag_value('type'       => 'rRNA'  );
$features[2]->add_tag_value('overlap_50' => 1       );

BAP::GeneMerge::tag_rna_overlap(\@features);

ok(!$features[0]->has_tag('delete_rrna_overlap'), 'glimmer3_1 is not deleted');

@features = (
             create_feature( 1,    1000, 1, 'glimmer3_1', 'glimmer3' ),
             create_feature( 991,  1090, 1, 'rfam_2',     'rfam'     ),
             create_feature( 951,  1050, 1, 'rfam_3',     'rfam'     ),
            );

$features[0]->add_tag_value('type'       => 'coding');
$features[1]->add_tag_value('type'       => 'rRNA'  );
$features[2]->add_tag_value('type'       => 'rRNA'  );

BAP::GeneMerge::tag_rna_overlap(\@features);

ok($features[0]->has_tag('delete_rrna_overlap'), 'glimmer3_1 is deleted');

@features = (
             create_feature( 1,    1000, 1, 'genemark_1', 'genemark' ),
             create_feature( 1000, 1099, 1, 'rnammer_1',  'rnammer'  ),
            );

$features[0]->add_tag_value('type'       => 'coding');
$features[1]->add_tag_value('type'       => 'rRNA'  );

BAP::GeneMerge::tag_rna_overlap(\@features);

ok($features[0]->has_tag('delete_rrna_overlap'), 'genemark_1 is deleted');


sub create_feature {

    my ($start, $end, $strand, $name, $source) = @_;


    unless (defined($source)) { $source = 'unknown'; }

    return Bio::SeqFeature::Generic->new(
                                         -seq_id       => 'TEST',
                                         -display_name => $name,
                                         -start        => $start,
                                         -end          => $end,
                                         -strand       => $strand,
                                         -source       => $source,
                                        );

}
