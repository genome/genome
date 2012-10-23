#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;

#use Test::More tests => 5;
use Test::More skip_all => "test data needs to be regenerated in dwdev";

BEGIN {
        use_ok('Genome::Model::Tools::Hgmi::SendToPap');
}


my $s = Genome::Model::Tools::Hgmi::SendToPap->create(
          'locus_tag' => "ROSINTL182DFT_test",
          'sequence_set_id' => 71,
	  'sequence_name'   => "Roseburia_intestinalis_ROSINTL182DFT_test_1.0.1.newb",
	  'organism_name'   => "Roseburia_intestinalis",
          'taxon_id' => 166486,
          'dev' => 1,
);

my $sequence = "ATGGGGGCCTGACCCGA";
my $result = $s->make_into_pep("-1", $sequence, "TEST.1");
is($result->seq, "SGQAP", 'make sequence into peptide');

ok($s->get_gene_peps(), 'run get_gene_peps()');
my $s1 = Genome::Model::Tools::Hgmi::SendToPap->create(
          'locus_tag' => "ROSINTL182DFT_test",
          'sequence_set_id' => 71,
	  'sequence_name'   => "Roseburia_intestinalis_ROSINTL182DFT_test_1.0.1.newb",
	  'organism_name'   => "Roseburia_intestinalis",
          'taxon_id' => 166486,
          'dev' => 1,
          'keep_pep' => 0,
          'gram_stain' => 'negative',
);

$s1->get_gene_peps();
my $file = $s->pep_file();
ok( -f $file, 'keep-pep option works');
my $size = -s $file;
ok($size gt 0, 'file is not empty');
#my $dir = $file =~ s/\/pap-{[a-zA-Z0-9.]}+.fa$//;
#diag( $dir);
diag( $file);
diag( $size);

remove \1, qq{ $file };

#$s->mgap_to_biosql(1);

# todo:


# $Id$
