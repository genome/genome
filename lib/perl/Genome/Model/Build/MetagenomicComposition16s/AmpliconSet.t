#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Test::More;

use_ok('Genome::Model::Build::MetagenomicComposition16s::AmpliconSet') or die;

my $tempdir = File::Temp::tempdir(CLEANUP => 1);
mkdir $tempdir.'/fasta';
mkdir $tempdir.'/classification';
my $file_base_name = 'H_GV-933124G-S.MOCK';
my $amplicon_set = Genome::Model::Build::MetagenomicComposition16s::AmpliconSet->create(
    name => 'V1_V3',
    primers => [qw/ c b a/],
    file_base_name => $file_base_name,
    directory => $tempdir,
    classifier => 'rdp2-1',
);
ok($amplicon_set, 'Created amplicon set');
is($amplicon_set->name, 'V1_V3', 'name');
is_deeply([$amplicon_set->primers], [qw/ c b a/], 'Primers');

my $writer = $amplicon_set->seq_writer_for('processed');
ok($writer, 'seq writer');
my $seq = { id => 'FZ0V7MM01A01AQ', seq => 'ATGC', qual => 'AAAA', desc => 'blah', };
$writer->write($seq);
$writer->close;
my $reader = $amplicon_set->seq_reader_for('processed');
ok($reader, 'seq reader');
my $new_seq = $reader->read;
is_deeply($new_seq, $seq, 'seq matches');

is($amplicon_set->file_base_name, 'H_GV-933124G-S.MOCK', 'file base name');
is($amplicon_set->fasta_dir, $tempdir.'/fasta', 'fasta dir base name');
is($amplicon_set->processed_fasta_file, $tempdir.'/fasta/'.$file_base_name.'.V1_V3.processed.fasta', 'processed fasta file');
is($amplicon_set->processed_qual_file, $tempdir.'/fasta/'.$file_base_name.'.V1_V3.processed.fasta.qual', 'processed qual file');
is($amplicon_set->chimera_free_fasta_file, $tempdir.'/fasta/'.$file_base_name.'.V1_V3.chimera_free.fasta', 'chimera free fasta file');
is($amplicon_set->chimera_free_qual_file, $tempdir.'/fasta/'.$file_base_name.'.V1_V3.chimera_free.fasta.qual', 'chimera free qual file');
is($amplicon_set->oriented_fasta_file, $tempdir.'/fasta/'.$file_base_name.'.V1_V3.oriented.fasta', 'oriented fasta file');
is($amplicon_set->oriented_qual_file, $tempdir.'/fasta/'.$file_base_name.'.V1_V3.oriented.fasta.qual', 'oriented qual file');

my %input = $amplicon_set->amplicon_iterator_input_fasta_and_qual;
is_deeply(
    {file => $amplicon_set->processed_fasta_file, qual_file => $amplicon_set->processed_qual_file,},
    \%input,
    'amplicon iterator input fasta is processed',
);
ok(
    Genome::Sys->copy_file($amplicon_set->processed_fasta_file, $amplicon_set->chimera_free_fasta_file),
    'copy processed fasta to chimera free',
);
ok(
    Genome::Sys->copy_file($amplicon_set->processed_qual_file, $amplicon_set->chimera_free_qual_file),
    'copy processed qual to chimera free',
);
%input = $amplicon_set->amplicon_iterator_input_fasta_and_qual;
is_deeply(
    {file => $amplicon_set->chimera_free_fasta_file, qual_file => $amplicon_set->chimera_free_qual_file,},
    \%input,
    'amplicon iterator input fasta is chimera free',
);

is($amplicon_set->chimera_dir, $tempdir.'/chimera', 'chimera dir');
is($amplicon_set->chimera_file, $tempdir.'/chimera/'.$file_base_name.'.V1_V3.chimera', 'chimera file');

my $classification_file = $amplicon_set->classification_file;
is($classification_file, $tempdir.'/classification/H_GV-933124G-S.MOCK.V1_V3.rdp2-1', 'classification file');
my $fh = Genome::Sys->open_file_for_writing($classification_file);
$fh->print("FZ0V7MM01A01AQ;-;Root:1.0;Bacteria:1.0;Fusobacteria:1.0;Fusobacteria:1.0;Fusobacteriales:1.0;Fusobacteriaceae:1.0;Fusobacterium:1.0;\n");
$fh->close;

my $amplicon = $amplicon_set->next_amplicon;
ok($amplicon, 'amplicon');
is($amplicon->{name}, 'FZ0V7MM01A01AQ', 'amplicon name');
ok($amplicon->{seq}, 'amplicon seq');
ok($amplicon->{classification}, 'amplicon classification');

#print "$tempdir\n"; <STDIN>;
done_testing();
exit;

