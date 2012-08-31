#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper;
use File::Compare;
use File::Temp;
use Storable;
use Test::More;

use_ok('Genome::Model::Tools::Consed::AceReader') or die;
use_ok('Genome::Model::Tools::Consed::AceWriter') or die;

my $ace_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Consed/M_BB0392D19.ace';
my $reader = Genome::Model::Tools::Consed::AceReader->create(
    file => $ace_file,
);
ok($reader, 'create ace reader');

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $ace_out = $tmpdir.'/ace';
my $writer = Genome::Model::Tools::Consed::AceWriter->create(
    file => $ace_out,
);
ok($writer, 'create ace ewriter');
while ( my $ctg = $reader->next_contig ) {
    #print Dumper({ map { $_ => $ctg->{$_} } (qw/ reads /)});
    ok($writer->add_contig(contig => $ctg), 'add ctg '.$ctg->{name});
}
ok($writer->add_contig_tags($reader->contig_tags), 'add contig tags');
ok($writer->add_assembly_tags($reader->assembly_tags), 'add assembly tags');
ok($writer->close, 'close the ewriter');
is(File::Compare::compare($ace_out, $ace_file), 0, 'ace file matches');

#print "gvimdiff $ace_out $ace_file\n"; <STDIN>;
done_testing();
exit;

